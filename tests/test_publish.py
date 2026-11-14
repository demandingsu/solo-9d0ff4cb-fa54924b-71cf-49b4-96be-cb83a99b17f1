"""Regression tests for paired-result publishing, rollback and stop semantics.

Covers the failure mode where the R2 gzip publish failed after R1 had
already landed: the worker used to die leaving the job stuck in
``running`` (or ``stopping``) with a downloadable R1 half-result.
"""

from __future__ import annotations

import gzip
import io
import threading
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app.main import create_app
from app.qc import (
    JobStopped,
    PublishFaults,
    publish_outputs,
)


def fastq_record(name: str, seq: str, qual: str) -> str:
    return f"@{name}\n{seq}\n+\n{qual}\n"


def gz_bytes(text: str) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        handle.write(text.encode("ascii"))
    return buffer.getvalue()


def gz_pair(r1: str, r2: str) -> tuple[bytes, bytes]:
    return gz_bytes(r1), gz_bytes(r2)


@pytest.fixture()
def client(tmp_path: Path):
    app = create_app(data_dir=tmp_path / "data")
    with TestClient(app) as test_client:
        yield test_client


def submit(
    client: TestClient,
    r1: bytes,
    r2: bytes,
    threshold: int = 20,
    **extra,
) -> dict:
    return client.post(
        "/jobs",
        files={
            "r1": ("R1.fastq.gz", r1, "application/gzip"),
            "r2": ("R2.fastq.gz", r2, "application/gzip"),
        },
        data={"min_quality": str(threshold), **extra},
    )


def wait_for(client: TestClient, job_id: str, statuses, timeout=20.0) -> dict:
    deadline = time.monotonic() + timeout
    payload = {}
    while time.monotonic() < deadline:
        payload = client.get(f"/jobs/{job_id}").json()
        if payload["status"] in statuses:
            return payload
        time.sleep(0.01)
    raise AssertionError(f"job {job_id} did not reach {statuses}: {payload}")


def read_gz(data: bytes) -> str:
    return gzip.decompress(data).decode("ascii")


HIGH = "I" * 4
PAIR_R1 = fastq_record("a/1", "ACGT", HIGH)
PAIR_R2 = fastq_record("a/2", "ACGT", HIGH)


def job_dir(client: TestClient, job_id: str) -> Path:
    return client.app.state.manager.jobs_dir / job_id


# ---------------------------------------------------------------------------
# Engine-level publish unit tests
# ---------------------------------------------------------------------------


def test_publish_rolls_back_when_r2_fails(tmp_path: Path):
    part1, part2 = tmp_path / "r1.part", tmp_path / "r2.part"
    out1, out2 = tmp_path / "r1.gz", tmp_path / "r2.gz"
    part1.write_bytes(b"r1-data")
    part2.write_bytes(b"r2-data")

    event = threading.Event()
    with pytest.raises(Exception) as excinfo:
        publish_outputs(
            part1, part2, out1, out2, event,
            faults=PublishFaults(fail_r2=True),
        )
    assert excinfo.value.code == "publish_failed"

    # No half-published result and no temporary file survives.
    for path in (part1, part2, out1, out2):
        assert not path.exists(), f"{path} must be removed after failed publish"


def test_publish_stop_before_any_move_cleans_up(tmp_path: Path):
    part1, part2 = tmp_path / "r1.part", tmp_path / "r2.part"
    out1, out2 = tmp_path / "r1.gz", tmp_path / "r2.gz"
    part1.write_bytes(b"r1-data")
    part2.write_bytes(b"r2-data")

    event = threading.Event()
    event.set()
    with pytest.raises(JobStopped):
        publish_outputs(part1, part2, out1, out2, event)
    for path in (part1, part2, out1, out2):
        assert not path.exists()


def test_publish_stop_between_moves_removes_staged_r1(tmp_path: Path):
    """A stop landing after the R1 move still leaves nothing published."""
    part1, part2 = tmp_path / "r1.part", tmp_path / "r2.part"
    out1, out2 = tmp_path / "r1.gz", tmp_path / "r2.gz"
    part1.write_bytes(b"r1-data")
    part2.write_bytes(b"r2-data")

    event = threading.Event()

    real_replace = __import__("os").replace

    def stop_after_first_replace(src, dst):
        result = real_replace(src, dst)
        if Path(dst) == out1:
            event.set()
        return result

    import os

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.setattr(os, "replace", stop_after_first_replace)
        with pytest.raises(JobStopped):
            publish_outputs(part1, part2, out1, out2, event)

    for path in (part1, part2, out1, out2):
        assert not path.exists()


def test_publish_slow_window_is_interruptible(tmp_path: Path):
    part1, part2 = tmp_path / "r1.part", tmp_path / "r2.part"
    out1, out2 = tmp_path / "r1.gz", tmp_path / "r2.gz"
    part1.write_bytes(b"r1-data")
    part2.write_bytes(b"r2-data")

    event = threading.Event()

    def stop_soon():
        time.sleep(0.2)
        event.set()

    thread = threading.Thread(target=stop_soon)
    thread.start()
    started = time.monotonic()
    with pytest.raises(JobStopped):
        publish_outputs(
            part1, part2, out1, out2, event,
            faults=PublishFaults(slow_seconds=30.0),
        )
    thread.join()
    # The 30 s sleep must be interrupted promptly, not run to completion.
    assert time.monotonic() - started < 2.0
    for path in (part1, part2, out1, out2):
        assert not path.exists()


# ---------------------------------------------------------------------------
# API-level paired publish behaviour
# ---------------------------------------------------------------------------


def test_r2_publish_failure_fails_job_and_removes_everything(client):
    response = submit(
        client, *gz_pair(PAIR_R1, PAIR_R2), fail_r2_publish="true"
    )
    assert response.status_code == 202, response.text
    job_id = response.json()["job_id"]

    done = wait_for(client, job_id, {"succeeded", "failed", "stopped"})
    assert done["status"] == "failed", done
    assert done["error"]["code"] == "publish_failed", done
    assert done["kept_pairs"] is None
    assert done["downloads"] is None
    assert done["phase"] is None

    # No leftover R1 half-result or .part files anywhere on disk.
    directory = job_dir(client, job_id)
    assert not directory.exists(), f"leftover files: {list(directory.parent.glob('*'))}"

    for side in ("r1", "r2"):
        resp = client.get(f"/jobs/{job_id}/download/{side}")
        assert resp.status_code == 409

    # A failed job cannot be stopped afterwards.
    assert client.post(f"/jobs/{job_id}/stop").status_code == 409


def test_failed_publish_job_is_not_stuck_running(client):
    """The worker must always settle the job into a terminal state."""
    job_id = submit(
        client, *gz_pair(PAIR_R1, PAIR_R2), fail_r2_publish="1"
    ).json()["job_id"]
    # Give generous time; the job must reach a terminal state and stay there.
    done = wait_for(client, job_id, {"succeeded", "failed", "stopped"})
    assert done["status"] == "failed"
    time.sleep(0.1)
    again = client.get(f"/jobs/{job_id}").json()
    assert again["status"] == "failed"


def test_stop_during_slow_publish_waits_then_reports_stopped(client):
    job_id = submit(
        client, *gz_pair(PAIR_R1, PAIR_R2), slow_publish_seconds="30"
    ).json()["job_id"]

    # Wait until the worker has entered the publish phase.
    deadline = time.monotonic() + 10
    phase = None
    while time.monotonic() < deadline:
        payload = client.get(f"/jobs/{job_id}").json()
        phase = payload.get("phase")
        if payload["status"] == "running" and phase == "publishing":
            break
        time.sleep(0.01)
    assert phase == "publishing", payload

    stop_resp = client.post(f"/jobs/{job_id}/stop")
    assert stop_resp.status_code == 200
    assert stop_resp.json()["status"] == "stopping"

    # The stop must be served quickly even though the publish was asked to
    # sleep for 30 s, and the final state must be stopped (never succeeded).
    started = time.monotonic()
    done = wait_for(client, job_id, {"stopped", "succeeded", "failed"}, timeout=10)
    assert time.monotonic() - started < 5.0
    assert done["status"] == "stopped", done
    assert done["downloads"] is None
    assert not job_dir(client, job_id).exists()
    for side in ("r1", "r2"):
        assert client.get(f"/jobs/{job_id}/download/{side}").status_code == 409


def test_status_and_stop_of_other_jobs_responsive_during_slow_publish(client):
    """A slow publish of job A must not block queries/stops of job B."""
    slow_id = submit(
        client, *gz_pair(PAIR_R1, PAIR_R2), slow_publish_seconds="20"
    ).json()["job_id"]

    # Wait for the slow job to be publishing.
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        payload = client.get(f"/jobs/{slow_id}").json()
        if payload.get("phase") == "publishing":
            break
        time.sleep(0.01)
    assert payload["phase"] == "publishing"

    # Another, normal job runs to completion while the first is publishing.
    other = submit(
        client,
        *gz_pair(
            fastq_record("b/1", "ACGT", HIGH),
            fastq_record("b/2", "ACGT", HIGH),
        ),
    ).json()["job_id"]
    other_done = wait_for(client, other, {"succeeded", "failed"})
    assert other_done["status"] == "succeeded", other_done

    # Polling the slow job and listing jobs stays fast.
    latencies = []
    for _ in range(5):
        start = time.monotonic()
        body = client.get(f"/jobs/{slow_id}").json()
        latencies.append(time.monotonic() - start)
        assert body["status"] in {"running", "stopping"}
        time.sleep(0.02)
    assert max(latencies) < 1.0, latencies

    # Stopping the slow job is accepted immediately and settles promptly.
    start = time.monotonic()
    stop = client.post(f"/jobs/{slow_id}/stop")
    assert stop.status_code == 200
    assert stop.json()["status"] == "stopping"
    done = wait_for(client, slow_id, {"stopped"}, timeout=10)
    assert done["status"] == "stopped"
    assert time.monotonic() - start < 5.0

    # The unrelated successful job is still downloadable.
    assert client.get(f"/jobs/{other}/download/r1").status_code == 200
    assert client.get(f"/jobs/{other}/download/r2").status_code == 200


def test_slow_publish_without_stop_still_succeeds_and_downloads(client):
    job_id = submit(
        client, *gz_pair(PAIR_R1, PAIR_R2), slow_publish_seconds="1"
    ).json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed", "stopped"}, timeout=15)
    assert done["status"] == "succeeded", done
    assert done["phase"] is None
    r1 = client.get(f"/jobs/{job_id}/download/r1")
    r2 = client.get(f"/jobs/{job_id}/download/r2")
    assert r1.status_code == r2.status_code == 200
    assert read_gz(r1.content) == PAIR_R1
    assert read_gz(r2.content) == PAIR_R2
    # Only the two final gzip files remain; inputs and parts are gone.
    remaining = sorted(p.name for p in job_dir(client, job_id).iterdir())
    assert remaining == ["R1.kept.fastq.gz", "R2.kept.fastq.gz"]


def test_invalid_fault_fields_rejected(client):
    resp = submit(
        client, *gz_pair(PAIR_R1, PAIR_R2), fail_r2_publish="maybe"
    )
    assert resp.status_code == 422
    resp = submit(
        client, *gz_pair(PAIR_R1, PAIR_R2), slow_publish_seconds="not-a-number"
    )
    assert resp.status_code == 422
    resp = submit(
        client, *gz_pair(PAIR_R1, PAIR_R2), slow_publish_seconds="999999"
    )
    assert resp.status_code == 422


def test_publish_failure_via_environment(monkeypatch, tmp_path: Path):
    monkeypatch.setenv("QC_FAIL_R2_PUBLISH", "true")
    app = create_app(data_dir=tmp_path / "data")
    with TestClient(app) as env_client:
        job_id = submit(
            env_client, *gz_pair(PAIR_R1, PAIR_R2)
        ).json()["job_id"]
        done = wait_for(env_client, job_id, {"failed"})
        assert done["error"]["code"] == "publish_failed"
        assert not (env_client.app.state.manager.jobs_dir / job_id).exists()
