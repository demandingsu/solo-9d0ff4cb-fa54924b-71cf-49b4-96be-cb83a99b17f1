from __future__ import annotations

import gzip
import io
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app.main import create_app


def fastq_record(name: str, seq: str, qual: str) -> str:
    return f"@{name}\n{seq}\n+\n{qual}\n"


def gz_bytes(text: str) -> bytes:
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        handle.write(text.encode("ascii"))
    return buffer.getvalue()


def gz_pair(
    r1: str, r2: str, *, truncate: str | None = None
) -> tuple[bytes, bytes]:
    data1 = gz_bytes(r1)
    data2 = gz_bytes(r2)
    if truncate == "r1":
        data1 = data1[:-6]
    elif truncate == "r2":
        data2 = data2[:-6]
    return data1, data2


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
    response = client.post(
        "/jobs",
        files={
            "r1": ("R1.fastq.gz", r1, "application/gzip"),
            "r2": ("R2.fastq.gz", r2, "application/gzip"),
        },
        data={"min_quality": str(threshold), **extra},
    )
    return response


def wait_for(client: TestClient, job_id: str, statuses, timeout=15.0) -> dict:
    import time

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        payload = client.get(f"/jobs/{job_id}").json()
        if payload["status"] in statuses:
            return payload
        time.sleep(0.01)
    raise AssertionError(f"job {job_id} did not reach {statuses}: {payload}")


def read_gz(data: bytes) -> str:
    return gzip.decompress(data).decode("ascii")


# ---------------------------------------------------------------------------
# Happy path / filtering rules
# ---------------------------------------------------------------------------

HIGH = "I" * 4


def test_keeps_pairs_and_downloads(client):
    r1 = fastq_record("read1/1", "ACGT", HIGH) + fastq_record(
        "read2/1", "TTTT", HIGH
    )
    r2 = fastq_record("read1/2", "GGGG", HIGH) + fastq_record(
        "read2/2", "CCCC", HIGH
    )
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    assert response.status_code == 202, response.text
    job = response.json()
    assert job["status"] in {"pending", "running", "succeeded"}

    done = wait_for(client, job["job_id"], {"succeeded", "failed"})
    assert done["status"] == "succeeded", done
    assert done["total_pairs"] == 2
    assert done["kept_pairs"] == 2
    assert done["discarded_pairs"] == 0
    assert done["downloads"] == {
        "r1": f"/jobs/{job['job_id']}/download/r1",
        "r2": f"/jobs/{job['job_id']}/download/r2",
    }

    file1 = client.get(done["downloads"]["r1"])
    file2 = client.get(done["downloads"]["r2"])
    assert file1.status_code == 200 and file2.status_code == 200
    assert file1.headers["content-type"] == "application/gzip"
    assert read_gz(file1.content) == r1
    assert read_gz(file2.content) == r2


def test_discards_whole_pair_when_either_end_low(client):
    # Pair 1: R1 low quality; pair 2: R2 low quality; pair 3: both fine.
    r1 = (
        fastq_record("a/1", "ACGT", "III#")
        + fastq_record("b/1", "ACGT", HIGH)
        + fastq_record("c/1", "ACGT", HIGH)
    )
    r2 = (
        fastq_record("a/2", "ACGT", HIGH)
        + fastq_record("b/2", "ACGT", "II#I")
        + fastq_record("c/2", "ACGT", HIGH)
    )
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed"})
    assert done["status"] == "succeeded", done
    assert done["total_pairs"] == 3
    assert done["kept_pairs"] == 1
    assert done["discarded_pairs"] == 2

    kept = read_gz(client.get(f"/jobs/{job_id}/download/r1").content)
    assert kept == fastq_record("c/1", "ACGT", HIGH)
    kept2 = read_gz(client.get(f"/jobs/{job_id}/download/r2").content)
    assert kept2 == fastq_record("c/2", "ACGT", HIGH)


def test_preserves_original_order(client):
    names = [f"read{i}" for i in range(1, 8)]
    r1_parts, r2_parts = [], []
    for index, name in enumerate(names):
        # Odd-indexed pairs are discarded.
        qual = HIGH if index % 2 == 0 else "III#"
        r1_parts.append(fastq_record(f"{name}/1", "ACGT", qual))
        r2_parts.append(fastq_record(f"{name}/2", "ACGT", HIGH))
    response = submit(
        client, *gz_pair("".join(r1_parts), "".join(r2_parts)), threshold=20
    )
    job_id = response.json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed"})
    assert done["kept_pairs"] == 4
    assert done["discarded_pairs"] == 3
    kept = read_gz(client.get(f"/jobs/{job_id}/download/r1").content)
    expected = "".join(
        fastq_record(f"{name}/1", "ACGT", HIGH)
        for name in names[::2]
    )
    assert kept == expected


def test_threshold_zero_keeps_low_quality(client):
    r1 = fastq_record("a/1", "ACGT", "!!!!")
    r2 = fastq_record("a/2", "ACGT", "!!!!")
    response = submit(client, *gz_pair(r1, r2), threshold=0)
    job_id = response.json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed"})
    assert done["status"] == "succeeded"
    assert done["kept_pairs"] == 1


def test_threshold_boundary_is_inclusive(client):
    # '5' is Phred 20 in Phred+33; a threshold of 20 must keep it.
    r1 = fastq_record("a/1", "ACGT", "5555")
    r2 = fastq_record("a/2", "ACGT", "IIII")
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed"})
    assert done["kept_pairs"] == 1

    response = submit(client, *gz_pair(r1, r2), threshold=21)
    job_id = response.json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed"})
    assert done["kept_pairs"] == 0
    assert done["discarded_pairs"] == 1


def test_comments_after_names_are_ignored(client):
    r1 = fastq_record("read1/1 some comment", "ACGT", HIGH)
    r2 = fastq_record("read1/2 other comment", "ACGT", HIGH)
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed"})
    assert done["status"] == "succeeded", done
    assert done["kept_pairs"] == 1


# ---------------------------------------------------------------------------
# Mandatory failures
# ---------------------------------------------------------------------------


def assert_failed(client, job_id, code):
    done = wait_for(client, job_id, {"succeeded", "failed", "stopped"})
    assert done["status"] == "failed", done
    assert done["error"]["code"] == code, done
    assert done["kept_pairs"] is None
    for side in ("r1", "r2"):
        response = client.get(f"/jobs/{job_id}/download/{side}")
        assert response.status_code == 409


def test_name_mismatch_fails(client):
    r1 = fastq_record("read1/1", "ACGT", HIGH) + fastq_record(
        "read2/1", "ACGT", HIGH
    )
    r2 = fastq_record("read1/2", "ACGT", HIGH) + fastq_record(
        "read3/2", "ACGT", HIGH
    )
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "pair_name_mismatch")


def test_missing_suffix_fails(client):
    r1 = fastq_record("read1", "ACGT", HIGH)
    r2 = fastq_record("read1/2", "ACGT", HIGH)
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "invalid_name_suffix")


def test_wrong_suffix_side_fails(client):
    r1 = fastq_record("read1/2", "ACGT", HIGH)
    r2 = fastq_record("read1/1", "ACGT", HIGH)
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "invalid_name_suffix")


def test_record_count_mismatch_fails(client):
    r1 = fastq_record("read1/1", "ACGT", HIGH)
    r2 = fastq_record("read1/2", "ACGT", HIGH) + fastq_record(
        "read2/2", "ACGT", HIGH
    )
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "pair_count_mismatch")


def test_seq_qual_length_mismatch_fails(client):
    r1 = fastq_record("read1/1", "ACGT", "III")
    r2 = fastq_record("read1/2", "ACGT", HIGH)
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "seq_qual_length_mismatch")


def test_truncated_gzip_fails(client):
    r1 = fastq_record("read1/1", "ACGT", HIGH)
    r2 = fastq_record("read1/2", "ACGT", HIGH)
    data1, data2 = gz_pair(r1, r2, truncate="r2")
    response = submit(client, data1, data2, threshold=20)
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "gzip_corrupt")


def test_trailing_garbage_in_gzip_fails(client):
    r1 = fastq_record("read1/1", "ACGT", HIGH)
    r2 = fastq_record("read1/2", "ACGT", HIGH)
    data1, data2 = gz_pair(r1, r2)
    response = submit(client, data1 + b"junk", data2, threshold=20)
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "gzip_corrupt")


def test_not_gzip_fails(client):
    response = submit(
        client, b"@read1/1\nACGT\n+\nIIII\n", gz_bytes(fastq_record("read1/2", "ACGT", HIGH))
    )
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "gzip_corrupt")


def test_invalid_quality_char_fails(client):
    r1 = fastq_record("read1/1", "ACGT", "III\x7f")
    r2 = fastq_record("read1/2", "ACGT", HIGH)
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "invalid_quality_char")


def test_malformed_record_fails(client):
    # Only three lines for the record.
    r1 = "@read1/1\nACGT\n+\n"
    r2 = fastq_record("read1/2", "ACGT", HIGH)
    response = submit(client, *gz_pair(r1, r2), threshold=20)
    job_id = response.json()["job_id"]
    assert_failed(client, job_id, "malformed_record")


def test_empty_upload_rejected(client):
    response = client.post(
        "/jobs",
        files={
            "r1": ("R1.fastq.gz", b"", "application/gzip"),
            "r2": ("R2.fastq.gz", gz_bytes(fastq_record("a/2", "ACGT", HIGH)), "application/gzip"),
        },
        data={"min_quality": "20"},
    )
    assert response.status_code == 400


# ---------------------------------------------------------------------------
# Request validation and routing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("threshold", ["-1", "94", "abc", ""])
def test_invalid_threshold_rejected(client, threshold):
    response = submit(
        client, gz_bytes(fastq_record("a/1", "A", "I")), gz_bytes(fastq_record("a/2", "A", "I")), threshold=threshold
    )
    assert response.status_code == 422


def test_missing_threshold_rejected(client):
    response = client.post(
        "/jobs",
        files={
            "r1": ("R1.fastq.gz", gz_bytes("x"), "application/gzip"),
            "r2": ("R2.fastq.gz", gz_bytes("x"), "application/gzip"),
        },
    )
    assert response.status_code == 422


def test_missing_file_rejected(client):
    response = client.post(
        "/jobs",
        files={"r1": ("R1.fastq.gz", b"x", "application/gzip")},
        data={"min_quality": "20"},
    )
    assert response.status_code == 422


def test_unknown_job_404(client):
    assert client.get("/jobs/nope").status_code == 404
    assert client.post("/jobs/nope/stop").status_code == 404
    assert client.get("/jobs/nope/download/r1").status_code == 404


def test_download_before_completion_conflicts(client, monkeypatch):
    import threading

    import app.jobs as jobs_module

    release = threading.Event()
    entered = threading.Event()
    real_run_qc = jobs_module.run_qc

    def blocking_run_qc(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        return real_run_qc(*args, **kwargs)

    monkeypatch.setattr(jobs_module, "run_qc", blocking_run_qc)
    r1 = fastq_record("a/1", "A", "I")
    r2 = fastq_record("a/2", "A", "I")
    job_id = submit(client, *gz_pair(r1, r2), threshold=20).json()["job_id"]
    assert entered.wait(10)
    status = client.get(f"/jobs/{job_id}").json()
    assert status["status"] == "running"
    assert client.get(f"/jobs/{job_id}/download/r1").status_code == 409
    release.set()
    done = wait_for(client, job_id, {"succeeded", "failed"})
    assert done["status"] == "succeeded"


def test_docs_available(client):
    response = client.get("/docs")
    assert response.status_code == 200
    assert "swagger" in response.text.lower()
    assert client.get("/openapi.json").status_code == 200


# ---------------------------------------------------------------------------
# Stop semantics
# ---------------------------------------------------------------------------


def _blocked_job(client, monkeypatch, tmp_path):
    import threading

    import app.jobs as jobs_module

    release = threading.Event()
    entered = threading.Event()
    real_run_qc = jobs_module.run_qc

    def blocking_run_qc(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        return real_run_qc(*args, **kwargs)

    monkeypatch.setattr(jobs_module, "run_qc", blocking_run_qc)
    r1 = fastq_record("a/1", "ACGT", HIGH)
    r2 = fastq_record("a/2", "ACGT", HIGH)
    job_id = submit(client, *gz_pair(r1, r2), threshold=20).json()["job_id"]
    assert entered.wait(10)
    return job_id, release


def test_stop_waits_for_worker_and_cleans_before_stopped(
    client, monkeypatch, tmp_path
):
    import time

    job_id, release = _blocked_job(client, monkeypatch, tmp_path)
    response = client.post(f"/jobs/{job_id}/stop")
    assert response.status_code == 200
    assert response.json()["status"] == "stopping"

    # While the worker is still alive the job must not report "stopped".
    for _ in range(20):
        assert client.get(f"/jobs/{job_id}").json()["status"] == "stopping"
        time.sleep(0.005)

    release.set()
    done = wait_for(client, job_id, {"stopped", "succeeded", "failed"})
    assert done["status"] == "stopped", done
    assert done["kept_pairs"] is None

    data_dir = client.app.state.manager.jobs_dir / job_id
    assert not data_dir.exists(), "temporary files must be cleaned before 'stopped'"
    assert client.get(f"/jobs/{job_id}/download/r1").status_code == 409
    assert client.get(f"/jobs/{job_id}/download/r2").status_code == 409

    # Stop is idempotent once stopped.
    again = client.post(f"/jobs/{job_id}/stop")
    assert again.status_code == 200
    assert again.json()["status"] == "stopped"


def test_stop_after_success_rejected(client):
    r1 = fastq_record("a/1", "ACGT", HIGH)
    r2 = fastq_record("a/2", "ACGT", HIGH)
    job_id = submit(client, *gz_pair(r1, r2), threshold=20).json()["job_id"]
    done = wait_for(client, job_id, {"succeeded", "failed"})
    assert done["status"] == "succeeded"
    response = client.post(f"/jobs/{job_id}/stop")
    assert response.status_code == 409
    assert wait_for(client, job_id, {"succeeded"})["status"] == "succeeded"
    assert client.get(f"/jobs/{job_id}/download/r1").status_code == 200


def test_stop_winning_over_failure(client, monkeypatch, tmp_path):
    import threading

    import app.jobs as jobs_module
    from app.qc import QCError

    entered = threading.Event()
    release = threading.Event()

    def failing_after_stop(*args, **kwargs):
        entered.set()
        assert release.wait(10)
        stop_event = kwargs.get("stop_event")
        if stop_event is not None and stop_event.is_set():
            from app.qc import JobStopped

            raise JobStopped()
        raise QCError("gzip_corrupt", "simulated")

    monkeypatch.setattr(jobs_module, "run_qc", failing_after_stop)
    r1 = fastq_record("a/1", "ACGT", HIGH)
    r2 = fastq_record("a/2", "ACGT", HIGH)
    job_id = submit(client, *gz_pair(r1, r2), threshold=20).json()["job_id"]
    assert entered.wait(10)
    assert client.post(f"/jobs/{job_id}/stop").status_code == 200
    release.set()
    done = wait_for(client, job_id, {"stopped", "failed"})
    assert done["status"] == "stopped"
    assert done["error"] is None


def test_stop_during_processing_of_large_file(client):
    """A genuinely running large job can be stopped and never succeeds."""
    import time

    names = 20000
    chunk = "".join(
        fastq_record(f"read{i}/1", "A" * 100, "I" * 100)
        for i in range(names)
    )
    chunk2 = "".join(
        fastq_record(f"read{i}/2", "T" * 100, "I" * 100)
        for i in range(names)
    )
    job_id = submit(
        client, gz_bytes(chunk), gz_bytes(chunk2), threshold=20
    ).json()["job_id"]
    time.sleep(0.01)
    response = client.post(f"/jobs/{job_id}/stop")
    assert response.status_code in {200, 409}
    if response.status_code == 200:
        done = wait_for(client, job_id, {"stopped", "succeeded", "failed"}, timeout=30)
        # The stop was accepted at some point; if success was already
        # published, the stop must have been rejected with 409.
        if response.json()["status"] != "succeeded":
            assert done["status"] == "stopped", done


def test_stop_before_worker_starts(client, monkeypatch):
    import threading

    import app.jobs as jobs_module

    release = threading.Event()
    real_run_qc = jobs_module.run_qc
    started = threading.Event()

    def slow_start(*args, **kwargs):
        started.set()
        assert release.wait(10)
        return real_run_qc(*args, **kwargs)

    monkeypatch.setattr(jobs_module, "run_qc", slow_start)
    r1 = fastq_record("a/1", "ACGT", HIGH)
    r2 = fastq_record("a/2", "ACGT", HIGH)
    job_id = submit(client, *gz_pair(r1, r2), threshold=20).json()["job_id"]
    assert started.wait(10)
    client.post(f"/jobs/{job_id}/stop")
    release.set()
    done = wait_for(client, job_id, {"stopped", "failed"})
    assert done["status"] == "stopped"


def test_failed_job_cannot_be_stopped(client):
    r1 = fastq_record("a/1", "ACGT", "III")
    r2 = fastq_record("a/2", "ACGT", HIGH)
    job_id = submit(client, *gz_pair(r1, r2), threshold=20).json()["job_id"]
    wait_for(client, job_id, {"failed"})
    assert client.post(f"/jobs/{job_id}/stop").status_code == 409


# ---------------------------------------------------------------------------
# Memory behaviour (guards against reading a whole file into memory)
# ---------------------------------------------------------------------------


def test_engine_does_not_read_whole_file(monkeypatch, tmp_path):
    from app import qc

    calls = {"readline": 0}
    original = qc._GzipLineReader.readline

    def counting(self):
        calls["readline"] += 1
        return original(self)

    monkeypatch.setattr(qc._GzipLineReader, "readline", counting)

    total = 500
    r1 = "".join(fastq_record(f"r{i}/1", "A", "I") for i in range(total))
    r2 = "".join(fastq_record(f"r{i}/2", "T", "I") for i in range(total))
    input1 = tmp_path / "R1.fastq.gz"
    input2 = tmp_path / "R2.fastq.gz"
    input1.write_bytes(gz_bytes(r1))
    input2.write_bytes(gz_bytes(r2))
    result = qc.run_qc(
        input1,
        input2,
        tmp_path / "o1.part",
        tmp_path / "o2.part",
        20,
        __import__("threading").Event(),
    )
    assert result.kept_pairs == total
    # Exactly 4 lines per record per side, plus one EOF-probe read per side.
    assert calls["readline"] == total * 4 * 2 + 2


def test_status_responds_while_large_job_runs(client):
    """Status queries stay responsive during streaming processing."""
    import time

    # ~120k pairs; large enough to still be running when polled.
    n = 60000
    body1 = gz_bytes(
        "".join(fastq_record(f"r{i}/1", "A" * 200, "I" * 200) for i in range(n))
    )
    body2 = gz_bytes(
        "".join(fastq_record(f"r{i}/2", "T" * 200, "I" * 200) for i in range(n))
    )
    job_id = submit(client, body1, body2, threshold=20).json()["job_id"]

    saw_running = False
    latencies = []
    deadline = time.monotonic() + 30
    while time.monotonic() < deadline:
        start = time.monotonic()
        payload = client.get(f"/jobs/{job_id}").json()
        latencies.append(time.monotonic() - start)
        if payload["status"] in {"running", "pending"}:
            saw_running = True
        if payload["status"] in {"succeeded", "failed", "stopped"}:
            break
        time.sleep(0.005)
    assert max(latencies) < 1.0, latencies
    final = wait_for(client, job_id, {"succeeded", "failed", "stopped"}, 30)
    if final["status"] == "succeeded":
        assert final["kept_pairs"] == n
