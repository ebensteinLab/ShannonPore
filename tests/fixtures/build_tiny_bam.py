"""Build a tiny FASTA + indexed BAM with MM/ML methylation tags for tests.

Produces a test fixture that modkit can extract per-CpG methylation calls
from. The reference contains a controlled number of CpG sites; each read
carries an MM tag (`C+m?` for 5mC, `C+h?` for 5hmC) and an ML tag with
known probabilities — letting downstream tests assert exact entropy /
MML values.

Usage from a pytest fixture:

    fa, bam = build_tiny_fixture(tmp_path)
"""

from __future__ import annotations

from pathlib import Path

import pysam

# Reference: a single 200-bp contig with 8 CpG sites at known positions.
# Constructed so every CpG is on the + strand and well-separated.
_CHROM = "chr_test"
_REF_SEQ = (
    "AAAA"
    "CG"
    "AAAA"  # CpG at pos 4
    "AAAA"
    "CG"
    "AAAA"  # CpG at pos 14
    "AAAA"
    "CG"
    "AAAA"  # CpG at pos 24
    "AAAA"
    "CG"
    "AAAA"  # CpG at pos 34
    "AAAA"
    "CG"
    "AAAA"  # CpG at pos 44
    "AAAA"
    "CG"
    "AAAA"  # CpG at pos 54
    "AAAA"
    "CG"
    "AAAA"  # CpG at pos 64
    "AAAA"
    "CG"
    "AAAA" + "A" * (200 - 80)  # CpG at pos 74
)
CPG_POSITIONS = [4, 14, 24, 34, 44, 54, 64, 74]


def write_reference_fasta(out_path: Path) -> Path:
    """Write a 200-bp single-contig FASTA and its .fai index."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write(f">{_CHROM}\n{_REF_SEQ}\n")
    pysam.faidx(str(out_path))
    return out_path


def _ml_score(prob: float) -> int:
    """Convert a [0,1] probability to a 0–255 ML-tag integer."""
    return max(0, min(255, int(round(prob * 255))))


def _build_one_read(
    bam: pysam.AlignmentFile,
    *,
    name: str,
    seq: str,
    methylation_pattern: list[float],
) -> None:
    """Write one read whose MM/ML tags encode the given 5mC pattern.

    `methylation_pattern[i]` is the 5mC probability for the i-th C in the
    read sequence. Length must equal the count of 'C' in `seq`.
    """
    a = pysam.AlignedSegment(bam.header)
    a.query_name = name
    a.query_sequence = seq
    a.flag = 0  # mapped, primary, forward
    a.reference_id = 0
    a.reference_start = 0
    a.cigar = [(0, len(seq))]  # 0 = M (match)
    a.mapping_quality = 60
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))

    # MM:Z: "C+m?,0,0,...;" with as many zero-skips as there are C bases.
    n_c = sum(1 for b in seq if b == "C")
    if n_c != len(methylation_pattern):
        raise ValueError(
            f"methylation_pattern length {len(methylation_pattern)} != "
            f"number of C bases in seq ({n_c})"
        )
    skips = ",".join(["0"] * n_c)
    mm_tag = f"C+m?,{skips};"
    ml_tag = [_ml_score(p) for p in methylation_pattern]

    a.set_tag("MM", mm_tag, value_type="Z")
    a.set_tag("ML", ml_tag)
    bam.write(a)


def write_methylation_bam(
    out_path: Path,
    fasta_path: Path,
    *,
    n_reads: int = 20,
    methylated_fraction: float = 0.5,
    seed: int = 0,
) -> Path:
    """Write `n_reads` reads each spanning the full reference. Each CpG
    is methylated in `methylated_fraction` of the reads (independently).
    """
    import random

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rng = random.Random(seed)
    fa = pysam.FastaFile(str(fasta_path))
    seq = fa.fetch(_CHROM)

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"LN": len(seq), "SN": _CHROM}],
    }

    with pysam.AlignmentFile(str(out_path), "wb", header=header) as bam:
        for i in range(n_reads):
            # Per-C methylation probabilities for THIS read.
            n_c = sum(1 for b in seq if b == "C")
            pattern: list[float] = []
            for _ in range(n_c):
                # Methylated → high prob, unmethylated → low prob
                pattern.append(0.9 if rng.random() < methylated_fraction else 0.05)
            _build_one_read(
                bam,
                name=f"read_{i:04d}",
                seq=seq,
                methylation_pattern=pattern,
            )

    pysam.sort("-o", str(out_path), str(out_path))
    pysam.index(str(out_path))
    return out_path


def build_tiny_fixture(
    tmp_path: Path,
    *,
    n_reads: int = 30,
    methylated_fraction: float = 0.5,
    seed: int = 0,
) -> tuple[Path, Path]:
    """Convenience wrapper: write FASTA + BAM into ``tmp_path``."""
    fa = write_reference_fasta(tmp_path / "tiny_ref.fa")
    bam = write_methylation_bam(
        tmp_path / "tiny.bam",
        fa,
        n_reads=n_reads,
        methylated_fraction=methylated_fraction,
        seed=seed,
    )
    return fa, bam


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        fa, bam = build_tiny_fixture(Path(td))
        print(f"Wrote: {fa}\n        {bam}")
