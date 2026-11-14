"""Generate deterministic local FASTQ.gz samples for manual testing.

Produces, under ``samples/``:

* ``good_R1.fastq.gz`` / ``good_R2.fastq.gz`` -- valid paired reads.
* ``bad_R1.fastq.gz`` / ``bad_R2.fastq.gz`` -- pairs that fail QC (low
  quality bases in one or both ends), used to show discarding behaviour.
* ``truncated_R2.fastq.gz`` -- a gzip stream that ends mid-record; pairing it
  with ``good_R1.fastq.gz`` makes the job fail.

Run from the repository root: ``python tools/make_samples.py``.
"""

from __future__ import annotations

import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SAMPLES = ROOT / "samples"


def record(name: str, sequence: str, qualities: str) -> str:
    assert len(sequence) == len(qualities)
    return f"@{name}\n{sequence}\n+\n{qualities}\n"


def write_gzip(path: Path, text: str, *, mtime: int = 0) -> None:
    with gzip.GzipFile(
        filename="", fileobj=open(path, "wb"), mode="wb", mtime=mtime
    ) as handle:
        handle.write(text.encode("ascii"))


def write_truncated(path: Path, text: str) -> None:
    """Write a valid gzip stream then cut it short mid-member."""
    raw = gzip.compress(text.encode("ascii"), mtime=0)
    # Keep the header and most of the body but drop the tail (CRC/size) and
    # start of the end-of-stream marker.
    truncated = raw[: max(10, len(raw) * 3 // 4)]
    path.write_bytes(truncated)


HIGH = "I"  # Phred 40
OK = "+"    # Phred 10
LOW = "#"   # Phred 2

GOOD_PAIRS = [
    ("read0001", "ACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIII"),
    ("read0002", "TTTTGGGGCCCCAAAATTTT", "IIIIIIIIIIIIIIIIIIII"),
    ("read0003", "GATTACAGATTACAGATTAC", "IIIIIIIIIIIIIIIIIIII"),
]

# Pair 4: R1 has one low-quality base (below threshold 20) -> whole pair drops.
# Pair 5: R2 has low quality -> whole pair drops.
# Pair 6: both ends low -> drop.
BAD_PAIRS = [
    (
        "read0004",
        "ACGTACGTACGTACGTACGT",
        "IIIIIIIIIIIIIIIIIII#",
        "ACGTACGTACGTACGTACGT",
        "IIIIIIIIIIIIIIIIIIII",
    ),
    (
        "read0005",
        "ACGTACGTACGTACGTACGT",
        "IIIIIIIIIIIIIIIIIIII",
        "ACGTACGTACGTACGTACGT",
        "IIIIIIIIIIIIIIIIII+#",
    ),
    (
        "read0006",
        "ACGTACGTACGTACGTACGT",
        "IIIIIIIIIIIIIIIIIII#",
        "ACGTACGTACGTACGTACGT",
        "IIIIIIIIIIIIIIIIIII#",
    ),
]

# A mixed file: first three pass, next three fail, last one passes. This
# demonstrates order preservation and per-pair filtering in one job.
MIXED_PASS_LAST = ("read0007", "ACGTACGTACGTACGTACGT", "IIIIIIIIIIIIIIIIIIII")


def build() -> None:
    SAMPLES.mkdir(parents=True, exist_ok=True)

    r1 = "".join(record(f"{name}/1", seq, qual) for name, seq, qual in GOOD_PAIRS)
    r2 = "".join(
        record(f"{name}/2", seq, qual) for name, seq, qual in GOOD_PAIRS
    )
    write_gzip(SAMPLES / "good_R1.fastq.gz", r1)
    write_gzip(SAMPLES / "good_R2.fastq.gz", r2)

    r1_parts, r2_parts = [], []
    for name, seq1, qual1, seq2, qual2 in BAD_PAIRS:
        r1_parts.append(record(f"{name}/1", seq1, qual1))
        r2_parts.append(record(f"{name}/2", seq2, qual2))
    write_gzip(SAMPLES / "bad_R1.fastq.gz", "".join(r1_parts))
    write_gzip(SAMPLES / "bad_R2.fastq.gz", "".join(r2_parts))

    mixed_r1, mixed_r2 = [], []
    for name, seq, qual in GOOD_PAIRS:
        mixed_r1.append(record(f"{name}/1", seq, qual))
        mixed_r2.append(record(f"{name}/2", seq, qual))
    for name, seq1, qual1, seq2, qual2 in BAD_PAIRS:
        mixed_r1.append(record(f"{name}/1", seq1, qual1))
        mixed_r2.append(record(f"{name}/2", seq2, qual2))
    name, seq, qual = MIXED_PASS_LAST
    mixed_r1.append(record(f"{name}/1", seq, qual))
    mixed_r2.append(record(f"{name}/2", seq, qual))
    write_gzip(SAMPLES / "mixed_R1.fastq.gz", "".join(mixed_r1))
    write_gzip(SAMPLES / "mixed_R2.fastq.gz", "".join(mixed_r2))

    # R1 valid; R2 truncated at the end of the gzip stream.
    write_gzip(SAMPLES / "truncated_R1.fastq.gz", r1)
    write_truncated(SAMPLES / "truncated_R2.fastq.gz", r2)

    print(f"samples written to {SAMPLES}")


if __name__ == "__main__":
    build()
