"""Streaming paired-end FASTQ quality-control engine.

The engine never loads a whole FASTQ file into memory.  Records are read one
pair at a time through ``gzip`` + a buffered reader, and kept pairs are written
straight into gzip output streams.

Rules implemented here:

* a record is exactly four lines;
* R1 names must end in ``/1`` and R2 names in ``/2``; after trimming the
  suffixes the two identifiers must be equal;
* the sequence and quality lines of a record must have equal length;
* quality characters must be printable Phred+33 (ASCII 33..126);
* a pair is discarded when either end contains a base with quality below the
  threshold;
* pair order is preserved;
* gzip corruption (truncated stream, bad CRC, trailing garbage, bad header)
  is a hard failure;
* a pair-count or identifier mismatch is a hard failure.
"""

from __future__ import annotations

import gzip
import io
import logging
import os
import sys
import threading
import time
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

logger = logging.getLogger(__name__)

# Phred+33 printable range.
MIN_QUALITY_BYTE = 33
MAX_QUALITY_BYTE = 126
MIN_PHRED = 0
MAX_PHRED = MAX_QUALITY_BYTE - MIN_QUALITY_BYTE  # 93

MAX_LINE_BYTES = 1 << 20  # 1 MiB is far beyond any sane FASTQ line
READ_BUFFER_BYTES = 1 << 20
PROGRESS_EVERY_PAIRS = 256


class QCError(Exception):
    """A deterministic input error that must fail the job."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code
        self.message = message


class JobStopped(Exception):
    """Raised inside the worker as soon as a stop request is observed."""


@dataclass(frozen=True)
class QCResult:
    total_pairs: int
    kept_pairs: int
    discarded_pairs: int


ProgressCallback = Callable[[int], None]


def _looks_like_truncation(exc: BaseException) -> bool:
    return isinstance(exc, EOFError)


class _GzipLineReader:
    """Bounded, stop-aware line reader over one gzip file."""

    def __init__(
        self, path: Path, side: str, stop_event: threading.Event
    ) -> None:
        self.path = path
        self.side = side
        self.stop_event = stop_event
        try:
            if path.stat().st_size == 0:
                raise QCError(
                    "gzip_corrupt", f"{side}: uploaded file is empty"
                )
            self._raw = open(path, "rb")
        except OSError as exc:
            raise QCError(
                "unreadable_input", f"{side}: cannot open input file: {exc}"
            ) from exc
        try:
            self._gz = gzip.GzipFile(fileobj=self._raw, mode="rb")
            self._reader = io.BufferedReader(
                self._gz, buffer_size=READ_BUFFER_BYTES
            )
        except Exception as exc:  # pragma: no cover - defensive
            self._raw.close()
            raise QCError(
                "gzip_corrupt", f"{side}: cannot open gzip stream: {exc}"
            ) from exc

    def _decode_error(self, exc: BaseException) -> QCError:
        if _looks_like_truncation(exc):
            return QCError(
                "gzip_corrupt",
                f"{self.side}: gzip stream is truncated "
                "(missing end-of-stream marker)",
            )
        return QCError(
            "gzip_corrupt", f"{self.side}: gzip stream is corrupt: {exc}"
        )

    def readline(self) -> bytes:
        if self.stop_event.is_set():
            raise JobStopped()
        try:
            line = self._reader.readline(MAX_LINE_BYTES)
            if len(line) >= MAX_LINE_BYTES and not line.endswith(b"\n"):
                raise QCError(
                    "malformed_record",
                    f"{self.side}: FASTQ line exceeds {MAX_LINE_BYTES} bytes",
                )
            return line
        except QCError:
            raise
        except (EOFError, gzip.BadGzipFile, zlib.error, OSError) as exc:
            raise self._decode_error(exc) from exc
        except ValueError as exc:
            raise QCError(
                "gzip_corrupt", f"{self.side}: invalid gzip stream: {exc}"
            ) from exc

    def close(self) -> None:
        for closer in (self._reader, self._raw):
            try:
                closer.close()
            except Exception:
                # Closing a partially consumed gzip stream can itself raise a
                # CRC/EOF error; the real error has already been reported.
                pass


class _GzipWriter:
    """Small wrapper that writes bytes through gzip and closes cleanly."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self._raw = open(path, "wb")
        self._gz = gzip.GzipFile(
            fileobj=self._raw, mode="wb", compresslevel=6, mtime=0
        )

    def write(self, data: bytes) -> None:
        self._gz.write(data)

    def close(self) -> None:
        try:
            self._gz.close()
        finally:
            self._raw.close()


def _strip_eol(line: bytes) -> bytes:
    if line.endswith(b"\n"):
        line = line[:-1]
        if line.endswith(b"\r"):
            line = line[:-1]
    elif line.endswith(b"\r"):
        line = line[:-1]
    return line


def _with_eol(line: bytes) -> bytes:
    """Return *line* terminated by exactly one ``\\n``."""
    body = _strip_eol(line)
    return body + b"\n"


def _parse_record_name(
    header: bytes, expected_suffix: bytes, side: str, pair_number: int
) -> bytes:
    if not header.startswith(b"@"):
        raise QCError(
            "malformed_record",
            f"{side}: record {pair_number} header does not start with '@'",
        )
    body = _strip_eol(header)[1:]
    tokens = body.split()
    if not tokens:
        raise QCError(
            "malformed_record",
            f"{side}: record {pair_number} has an empty read name",
        )
    name = tokens[0]
    if not name.endswith(expected_suffix):
        suffix = expected_suffix.decode("ascii")
        raise QCError(
            "invalid_name_suffix",
            f"{side}: record {pair_number} read name "
            f"{name.decode('utf-8', 'replace')!r} does not end with {suffix!r}",
        )
    identifier = name[: -len(expected_suffix)]
    if not identifier:
        raise QCError(
            "malformed_record",
            f"{side}: record {pair_number} read name has no identifier "
            "before the /1 or /2 suffix",
        )
    return identifier


def _check_record(
    lines: list[bytes],
    side: str,
    pair_number: int,
    expected_suffix: bytes,
) -> tuple[bytes, int]:
    """Validate one record and return ``(identifier, min_phred)``."""
    header, sequence_line, plus_line, quality_line = lines

    identifier = _parse_record_name(
        header, expected_suffix, side, pair_number
    )

    if not plus_line.startswith(b"+"):
        raise QCError(
            "malformed_record",
            f"{side}: record {pair_number} third line does not start with '+'",
        )

    sequence = _strip_eol(sequence_line)
    quality = _strip_eol(quality_line)
    if len(sequence) != len(quality):
        raise QCError(
            "seq_qual_length_mismatch",
            f"{side}: record {pair_number} sequence length "
            f"({len(sequence)}) differs from quality length ({len(quality)})",
        )

    min_phred = MAX_PHRED
    for byte in quality:
        if byte < MIN_QUALITY_BYTE or byte > MAX_QUALITY_BYTE:
            raise QCError(
                "invalid_quality_char",
                f"{side}: record {pair_number} contains a non-Phred+33 "
                f"quality byte (0x{byte:02x}); only ASCII 33..126 is allowed",
            )
        phred = byte - MIN_QUALITY_BYTE
        if phred < min_phred:
            min_phred = phred
    return identifier, min_phred


def run_qc(
    r1_path: Path,
    r2_path: Path,
    output_r1_part: Path,
    output_r2_part: Path,
    threshold: int,
    stop_event: threading.Event,
    *,
    progress_cb: ProgressCallback | None = None,
) -> QCResult:
    """Filter paired FASTQ files.

    ``output_*_part`` files are complete, valid gzip files when this function
    returns.  The caller renames them into place only after a successful
    return, so a partially written result is never published.
    """
    if not (MIN_PHRED <= threshold <= MAX_PHRED):
        raise QCError(
            "invalid_threshold",
            f"quality threshold must be between {MIN_PHRED} and {MAX_PHRED}",
        )

    reader1 = reader2 = writer1 = writer2 = None
    try:
        reader1 = _GzipLineReader(r1_path, "R1", stop_event)
        reader2 = _GzipLineReader(r2_path, "R2", stop_event)
        writer1 = _GzipWriter(output_r1_part)
        writer2 = _GzipWriter(output_r2_part)

        kept = 0
        discarded = 0
        pair_number = 0
        while True:
            if stop_event.is_set():
                raise JobStopped()

            header1 = reader1.readline()
            header2 = reader2.readline()
            if not header1 and not header2:
                break
            if not header1 or not header2:
                missing = "R2" if header1 else "R1"
                other = "R1" if missing == "R2" else "R2"
                raise QCError(
                    "pair_count_mismatch",
                    f"{other} has more records than {missing}; "
                    "the two ends contain different numbers of reads",
                )

            pair_number += 1
            lines1 = [
                header1,
                reader1.readline(),
                reader1.readline(),
                reader1.readline(),
            ]
            lines2 = [
                header2,
                reader2.readline(),
                reader2.readline(),
                reader2.readline(),
            ]
            for side, lines in (("R1", lines1), ("R2", lines2)):
                if any(line == b"" for line in lines[1:]):
                    raise QCError(
                        "malformed_record",
                        f"{side}: record {pair_number} has fewer than four "
                        "lines (truncated FASTQ record)",
                    )

            id1, min_q1 = _check_record(lines1, "R1", pair_number, b"/1")
            id2, min_q2 = _check_record(lines2, "R2", pair_number, b"/2")
            if id1 != id2:
                raise QCError(
                    "pair_name_mismatch",
                    f"record {pair_number}: R1 read name "
                    f"{id1.decode('utf-8', 'replace')!r} does not match R2 "
                    f"read name {id2.decode('utf-8', 'replace')!r}",
                )

            if min_q1 < threshold or min_q2 < threshold:
                discarded += 1
            else:
                kept += 1
                if stop_event.is_set():
                    raise JobStopped()
                writer1.write(b"".join(_with_eol(line) for line in lines1))
                writer2.write(b"".join(_with_eol(line) for line in lines2))

            if progress_cb is not None and (
                pair_number % PROGRESS_EVERY_PAIRS == 0
            ):
                progress_cb(pair_number)

        if stop_event.is_set():
            raise JobStopped()
        return QCResult(
            total_pairs=pair_number,
            kept_pairs=kept,
            discarded_pairs=discarded,
        )
    finally:
        # If the body already raised, keep that error: a secondary close
        # failure (e.g. a truncated gzip raises on close) must not mask it.
        failed = sys.exc_info()[0] is not None
        for reader in (reader1, reader2):
            if reader is not None:
                reader.close()
        for writer in (writer1, writer2):
            if writer is not None:
                try:
                    writer.close()
                except Exception:
                    if not failed:
                        raise


@dataclass(frozen=True)
class PublishFaults:
    """Deterministic fault injection for the paired-result publish step.

    ``fail_r2`` makes the R2 publish raise after R1 has been moved into
    place, exercising the rollback path.  ``slow_seconds`` stretches the
    publish so responsiveness/stop behaviour can be observed.  The injected
    slumber is interruptible by the job's stop event.
    """

    fail_r2: bool = False
    slow_seconds: float = 0.0


class PublishError(Exception):
    """A result file could not be published; the job must fail cleanly."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.code = "publish_failed"
        self.message = message


def _interruptible_sleep(seconds: float, stop_event: threading.Event) -> None:
    """Sleep *seconds* but return as soon as the stop event is observed."""
    deadline = time.monotonic() + seconds
    while True:
        if stop_event.is_set():
            return
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return
        stop_event.wait(timeout=min(0.05, remaining))


def _rollback_publish(part_r1: Path, part_r2: Path, final_r1: Path, final_r2: Path) -> None:
    """Undo a partially completed publish and remove all result artifacts.

    Best effort: any file that has already reached its final name is moved
    back to its ``.part`` name when possible, then every part/final path is
    deleted so no half-published result can ever be downloaded.
    """
    for final, part in ((final_r1, part_r1), (final_r2, part_r2)):
        try:
            if final.is_file():
                if part.exists():
                    final.unlink()
                else:
                    final.replace(part)
        except OSError:
            logger.warning("could not roll back published file %s", final)
    for path in (part_r1, part_r2, final_r1, final_r2):
        try:
            path.unlink(missing_ok=True)
        except OSError:  # pragma: no cover - best effort
            logger.warning("could not remove leftover result file %s", path)


def publish_outputs(
    part_r1: Path,
    part_r2: Path,
    final_r1: Path,
    final_r2: Path,
    stop_event: threading.Event,
    *,
    faults: PublishFaults | None = None,
) -> None:
    """Publish the paired gzip results atomically as a pair.

    Both files are moved into their download names inside this call.  If the
    R2 publish fails for any reason -- an injected fault or a real
    ``OSError`` -- the already-published R1 file is rolled back and every
    result/temporary file is removed, then :class:`PublishError` is raised so
    the caller can settle the job into an explainable failed terminal state.

    A stop request observed before either move, between the two moves, or
    during the (interruptible) slow-publish window raises
    :class:`JobStopped` after cleaning up; a job that was stopped can never
    end up published.
    """
    faults = faults or PublishFaults()

    if stop_event.is_set():
        _rollback_publish(part_r1, part_r2, final_r1, final_r2)
        raise JobStopped()
    if faults.slow_seconds > 0:
        _interruptible_sleep(faults.slow_seconds, stop_event)
        if stop_event.is_set():
            _rollback_publish(part_r1, part_r2, final_r1, final_r2)
            raise JobStopped()

    try:
        os.replace(part_r1, final_r1)
    except OSError as exc:
        _rollback_publish(part_r1, part_r2, final_r1, final_r2)
        raise PublishError(f"R1 result could not be published: {exc}") from exc

    if stop_event.is_set():
        _rollback_publish(part_r1, part_r2, final_r1, final_r2)
        raise JobStopped()

    if faults.fail_r2:
        _rollback_publish(part_r1, part_r2, final_r1, final_r2)
        raise PublishError(
            "R2 result publish failed after R1 was staged; rolled R1 back"
        )

    try:
        os.replace(part_r2, final_r2)
    except OSError as exc:
        _rollback_publish(part_r1, part_r2, final_r1, final_r2)
        raise PublishError(f"R2 result could not be published: {exc}") from exc
