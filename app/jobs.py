"""In-process job manager for the paired-end FASTQ QC service.

Only one API process is assumed.  Every job owns a private directory under
``DATA_DIR/jobs``.  Uploads and partial outputs live there, and are removed
before a stop or failure is reported.  Successful jobs keep only the two
filtered gzip files.

Concurrency contract
--------------------
* A stop requested while a job is pending/running atomically moves it to
  ``stopping`` and sets its event; it can never subsequently become
  ``succeeded``.
* A job that already reached ``succeeded``/``failed``/``stopped`` rejects
  stop requests (the API maps that to HTTP 409).
* ``stopping`` is only turned into ``stopped`` after the worker has exited and
  its temporary files have been removed.
"""

from __future__ import annotations

import logging
import os
import shutil
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .qc import (
    JobStopped,
    PublishError,
    PublishFaults,
    QCError,
    QCResult,
    publish_outputs,
    run_qc,
)

logger = logging.getLogger(__name__)

PENDING = "pending"
RUNNING = "running"
STOPPING = "stopping"
SUCCEEDED = "succeeded"
FAILED = "failed"
STOPPED = "stopped"

# Worker sub-phases (surfaced only while a job is active).
PHASE_PROCESSING = "processing"
PHASE_PUBLISHING = "publishing"

TERMINAL_STATUSES = frozenset({SUCCEEDED, FAILED, STOPPED})
STOPPABLE_STATUSES = frozenset({PENDING, RUNNING})

MAX_SLOW_PUBLISH_SECONDS = 300.0


def _env_flag(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


def _env_float(name: str) -> float:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return 0.0
    try:
        return max(0.0, float(raw))
    except ValueError:
        logger.warning("ignoring invalid %s=%r", name, raw)
        return 0.0


class JobNotFound(Exception):
    pass


class JobConflict(Exception):
    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.message = message


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class Job:
    job_id: str
    directory: Path
    r1_filename: str
    r2_filename: str
    min_quality: int
    publish_faults: PublishFaults = field(default_factory=PublishFaults)
    status: str = PENDING
    phase: str | None = None
    created_at: str = field(default_factory=_now)
    started_at: str | None = None
    finished_at: str | None = None
    pairs_processed: int = 0
    total_pairs: int | None = None
    kept_pairs: int | None = None
    discarded_pairs: int | None = None
    error_code: str | None = None
    error_message: str | None = None
    stop_event: threading.Event = field(default_factory=threading.Event)
    thread: threading.Thread | None = None

    # -- paths -----------------------------------------------------------
    @property
    def input_r1(self) -> Path:
        return self.directory / "R1.fastq.gz"

    @property
    def input_r2(self) -> Path:
        return self.directory / "R2.fastq.gz"

    @property
    def part_r1(self) -> Path:
        return self.directory / "R1.kept.fastq.gz.part"

    @property
    def part_r2(self) -> Path:
        return self.directory / "R2.kept.fastq.gz.part"

    @property
    def output_r1(self) -> Path:
        return self.directory / "R1.kept.fastq.gz"

    @property
    def output_r2(self) -> Path:
        return self.directory / "R2.kept.fastq.gz"

    @property
    def downloads_ready(self) -> bool:
        return self.status == SUCCEEDED and self.output_r1.is_file() and self.output_r2.is_file()

    def to_dict(self) -> dict[str, Any]:
        result: dict[str, Any] = {
            "job_id": self.job_id,
            "status": self.status,
            "phase": self.phase if self.status not in TERMINAL_STATUSES else None,
            "min_quality": self.min_quality,
            "r1_filename": self.r1_filename,
            "r2_filename": self.r2_filename,
            "created_at": self.created_at,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "pairs_processed": self.pairs_processed,
            "total_pairs": self.total_pairs,
            "kept_pairs": self.kept_pairs,
            "discarded_pairs": self.discarded_pairs,
            "error": (
                {"code": self.error_code, "message": self.error_message}
                if self.error_code
                else None
            ),
            "downloads": None,
        }
        if self.downloads_ready:
            result["downloads"] = {
                "r1": f"/jobs/{self.job_id}/download/r1",
                "r2": f"/jobs/{self.job_id}/download/r2",
            }
        return result


class JobManager:
    def __init__(self, data_dir: Path) -> None:
        self.data_dir = Path(data_dir)
        self.jobs_dir = self.data_dir / "jobs"
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        self._jobs: dict[str, Job] = {}
        self._lock = threading.Lock()
        # Process-wide publish fault defaults (overridable per job).
        self._default_faults = PublishFaults(
            fail_r2=_env_flag("QC_FAIL_R2_PUBLISH"),
            slow_seconds=min(_env_float("QC_SLOW_PUBLISH_SECONDS"), MAX_SLOW_PUBLISH_SECONDS),
        )
        self._remove_stale_directories()

    # -- lifecycle -------------------------------------------------------
    def _remove_stale_directories(self) -> None:
        """Drop directories left behind by a previous process.

        Job state is in memory only, so such directories can no longer be
        downloaded and would otherwise accumulate forever.
        """
        for child in self.jobs_dir.iterdir():
            try:
                if child.is_dir():
                    shutil.rmtree(child, ignore_errors=True)
                else:
                    child.unlink(missing_ok=True)
            except OSError:  # pragma: no cover - best effort
                logger.warning("could not remove stale path %s", child)

    def create_job(
        self,
        r1_filename: str,
        r2_filename: str,
        min_quality: int,
        publish_faults: PublishFaults | None = None,
    ) -> Job:
        job_id = uuid.uuid4().hex
        directory = self.jobs_dir / job_id
        directory.mkdir(parents=True)
        job = Job(
            job_id=job_id,
            directory=directory,
            r1_filename=r1_filename,
            r2_filename=r2_filename,
            min_quality=min_quality,
            publish_faults=publish_faults or self._default_faults,
        )
        with self._lock:
            self._jobs[job_id] = job
        return job

    def abort_upload(self, job_id: str) -> None:
        """Discard a job whose upload did not complete."""
        with self._lock:
            job = self._jobs.pop(job_id, None)
        if job is not None:
            shutil.rmtree(job.directory, ignore_errors=True)

    def start(self, job_id: str) -> None:
        job = self.get(job_id)
        thread = threading.Thread(
            target=self._run_job,
            args=(job,),
            name=f"qc-{job_id[:8]}",
            daemon=True,
        )
        with self._lock:
            if job.status != PENDING:
                raise JobConflict("job has already been started")
            job.thread = thread
        thread.start()

    def shutdown(self) -> None:
        with self._lock:
            jobs = list(self._jobs.values())
        for job in jobs:
            job.stop_event.set()
        for job in jobs:
            thread = job.thread
            if thread is not None and thread.is_alive():
                thread.join(timeout=5)

    # -- queries ---------------------------------------------------------
    def get(self, job_id: str) -> Job:
        with self._lock:
            job = self._jobs.get(job_id)
        if job is None:
            raise JobNotFound(job_id)
        return job

    def list_jobs(self) -> list[Job]:
        with self._lock:
            jobs = list(self._jobs.values())
        return sorted(jobs, key=lambda j: j.created_at, reverse=True)

    def snapshot(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                raise JobNotFound(job_id)
            return job.to_dict()

    def request_stop(self, job_id: str) -> dict[str, Any]:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                raise JobNotFound(job_id)
            if job.status in STOPPABLE_STATUSES:
                job.status = STOPPING
                job.stop_event.set()
            elif job.status == STOPPING:
                job.stop_event.set()
            elif job.status == STOPPED:
                pass  # idempotent
            else:  # succeeded / failed
                raise JobConflict(
                    f"job is already {job.status}; it can no longer be stopped"
                )
            return job.to_dict()

    def output_paths(self, job_id: str) -> tuple[Path, Path]:
        with self._lock:
            job = self._jobs.get(job_id)
            if job is None:
                raise JobNotFound(job_id)
            if job.status != SUCCEEDED:
                raise JobConflict(
                    f"job status is {job.status!r}; downloads are only "
                    "available after successful completion"
                )
            if not (job.output_r1.is_file() and job.output_r2.is_file()):
                raise JobConflict("result files are no longer available")
            return job.output_r1, job.output_r2

    # -- worker ----------------------------------------------------------
    def _run_job(self, job: Job) -> None:
        with self._lock:
            if job.status == STOPPING or job.stop_event.is_set():
                # A stop request landed while the worker was still starting.
                # The worker has now effectively exited, so it is responsible
                # for cleaning up and publishing ``stopped``.
                stop_early = True
            elif job.status == PENDING:
                job.status = RUNNING
                job.phase = PHASE_PROCESSING
                job.started_at = _now()
                stop_early = False
            else:  # pragma: no cover - upload aborted / defensive
                return
        if stop_early:
            self._finish_stopped(job)
            return

        try:
            result = run_qc(
                job.input_r1,
                job.input_r2,
                job.part_r1,
                job.part_r2,
                job.min_quality,
                job.stop_event,
                progress_cb=lambda pairs: self._set_progress(job, pairs),
            )
        except JobStopped:
            self._finish_stopped(job)
        except QCError as exc:
            self._finish_failed(job, exc.code, exc.message)
        except Exception as exc:  # pragma: no cover - unexpected bug
            logger.exception("unexpected error in job %s", job.job_id)
            self._finish_failed(job, "internal_error", str(exc))
        else:
            self._finish_completed(job, result)

    def _set_progress(self, job: Job, pairs: int) -> None:
        with self._lock:
            if pairs > job.pairs_processed:
                job.pairs_processed = pairs

    def _finish_completed(self, job: Job, result: QCResult) -> None:
        # Publish *outside* the job lock: a slow publish must not block
        # status queries or stop requests for this or any other job.  The
        # publish itself is stop-aware and rolls back on any failure.
        with self._lock:
            if job.status == STOPPING or job.stop_event.is_set():
                self._finish_stopped(job)
                return
            job.phase = PHASE_PUBLISHING

        try:
            publish_outputs(
                job.part_r1,
                job.part_r2,
                job.output_r1,
                job.output_r2,
                job.stop_event,
                faults=job.publish_faults,
            )
        except JobStopped:
            # Stop accepted while publishing.  The publish action has now
            # exited and cleaned up; only then may we report "stopped".
            self._finish_stopped(job)
            return
        except PublishError as exc:
            # One end could not be published.  publish_outputs already
            # removed every staged/temporary result, so settle the job into
            # an explainable failed terminal state (unless a stop raced and
            # won, in which case the stop contract takes precedence).
            with self._lock:
                stop_won = job.status == STOPPING or job.stop_event.is_set()
            if stop_won:
                self._finish_stopped(job)
            else:
                self._finish_failed(job, exc.code, exc.message)
            return
        except Exception as exc:  # pragma: no cover - defensive safety net
            logger.exception("unexpected publish error in job %s", job.job_id)
            with self._lock:
                stop_won = job.status == STOPPING or job.stop_event.is_set()
            if stop_won:
                self._finish_stopped(job)
            else:
                self._finish_failed(job, "publish_failed", str(exc))
            return

        # Both files are on disk.  Re-check the stop flag under the lock and
        # only then commit "succeeded": an accepted stop can never flip back
        # to success, and once "succeeded" is visible the pair is complete.
        with self._lock:
            if job.status == STOPPING or job.stop_event.is_set():
                stop_wins = True
            else:
                stop_wins = False
                job.total_pairs = result.total_pairs
                job.kept_pairs = result.kept_pairs
                job.discarded_pairs = result.discarded_pairs
                job.pairs_processed = result.total_pairs
                job.phase = None
                job.finished_at = _now()
                job.status = SUCCEEDED
        if stop_wins:
            self._finish_stopped(job)
            return
        self._cleanup_inputs(job)

    def _finish_failed(self, job: Job, code: str, message: str) -> None:
        with self._lock:
            if job.status == STOPPING or job.stop_event.is_set():
                stop_wins = True
            else:
                stop_wins = False
                job.status = FAILED
                job.phase = None
                job.error_code = code
                job.error_message = message
                job.finished_at = _now()
        if stop_wins:
            self._finish_stopped(job)
            return
        self._cleanup_all(job)

    def _finish_stopped(self, job: Job) -> None:
        # Temporary inputs, partial outputs and any staged results are
        # removed *before* the job is reported as stopped, as required by
        # the stop contract.  The publish action has already exited by now.
        self._cleanup_all(job)
        with self._lock:
            job.phase = None
            job.total_pairs = None
            job.kept_pairs = None
            job.discarded_pairs = None
            job.error_code = None
            job.error_message = None
            job.finished_at = _now()
            job.status = STOPPED

    def _cleanup_inputs(self, job: Job) -> None:
        for path in (job.input_r1, job.input_r2, job.part_r1, job.part_r2):
            try:
                path.unlink(missing_ok=True)
            except OSError:  # pragma: no cover - best effort
                logger.warning("could not remove %s", path)

    def _cleanup_all(self, job: Job) -> None:
        shutil.rmtree(job.directory, ignore_errors=True)
