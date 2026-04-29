"""Golden-file test for the File Preparation pipeline.

Builds tiny BAMs and TSVs whose expected ME / MML / coverage bedgraphs
are known *a priori* (from first principles, not from a recorded run),
then runs the orchestrator end-to-end and asserts that the bedgraphs
match the predicted values exactly.

Why three fixtures?
-------------------
1. **Fully methylated**  — every CpG on every read is methylated. Each
   bin's MML must be 1.0 and ME must be 0.0 (single, identical pattern
   across reads → zero Shannon entropy).
2. **Fully unmethylated** — mirror image: MML = 0.0, ME = 0.0.
3. **Single-read TSV path** — the same fixture but fed in as a
   pre-computed TSV, exercising the input_kind="tsv" branch.

These three together cover:
  * BAM input → modkit extract → entropy (the full file-prep pipeline)
  * TSV input → entropy (the modkit-skipping branch)
  * Numerical correctness (golden MML / ME / coverage values)
  * Bedgraph schema (chrom · start · end · value, tab-separated)
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pysam = pytest.importorskip("pysam")

from tests.fixtures.build_tiny_bam import (  # noqa: E402
    CPG_POSITIONS,
    write_methylation_bam,
    write_reference_fasta,
)

V4_DIR = Path(__file__).resolve().parents[2]
MODKIT = shutil.which("modkit")

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(MODKIT is None, reason="modkit binary not on PATH"),
]

# Pipeline parameters used across the golden tests. k=2 over 8 CpGs = 4 bins.
K = 2
MIN_COV = 4
N_READS = 12
THREADS = 2


# ─────────────────────────── helpers ──────────────────────────────────────

def _build_uniform_bam(out_dir: Path, *, methylated: bool) -> tuple[Path, Path]:
    """Build a FASTA + BAM where every read is uniformly (un)methylated."""
    fa = write_reference_fasta(out_dir / "ref.fa")
    bam = write_methylation_bam(
        out_dir / "uniform.bam", fa,
        n_reads=N_READS,
        methylated_fraction=1.0 if methylated else 0.0,
        seed=0,
    )
    return fa, bam


def _expected_bin_edges() -> list[tuple[int, int]]:
    """Same logic the WG pipeline uses: bin = (cpg[i], cpg[i+k-1] + 2)."""
    edges: list[tuple[int, int]] = []
    n_full = len(CPG_POSITIONS) // K
    for i in range(n_full):
        start = CPG_POSITIONS[i * K]
        end = CPG_POSITIONS[i * K + (K - 1)] + 2
        edges.append((start, end))
    return edges


def _read_bedgraph(p: Path) -> pd.DataFrame:
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(columns=["chrom", "start", "end", "value"])
    return pd.read_csv(
        p, sep="\t", header=None, names=["chrom", "start", "end", "value"],
    )


def _orchestrate(fa: Path, *, bam: Path | None, tsv: Path | None,
                 out_dir: Path, label: str) -> object:
    """Drive the v4 orchestrator without going through the CLI."""
    from src.pipelines.orchestrator import run_pipeline
    from src.state import SampleSpec

    if bam is not None:
        spec = SampleSpec(label=label, input_kind="bam", bam_path=bam)
    else:
        assert tsv is not None
        spec = SampleSpec(label=label, input_kind="tsv", tsv_path=tsv)

    target, _ = run_pipeline(
        target_spec=spec, control_spec=None,
        out_dir=out_dir, fasta=fa,
        entropy_mode="true_mc",
        cpgs_per_bin=K, min_coverage=MIN_COV, methyl_threshold=0.5,
        threads=THREADS, chroms="chr_test",
    )
    return target


# ─────────────────────────── tests ────────────────────────────────────────

def test_bam_fully_methylated_yields_mml_one_and_zero_entropy(
    tmp_path: Path,
) -> None:
    """Fully methylated reads → every bin MML=1.0, ME=0.0, coverage=N_READS."""
    fa, bam = _build_uniform_bam(tmp_path / "fixture", methylated=True)
    out_dir = tmp_path / "out"
    result = _orchestrate(fa, bam=bam, tsv=None, out_dir=out_dir, label="meth")

    edges = _expected_bin_edges()

    cov = _read_bedgraph(result.coverage_bedgraph)
    mml = _read_bedgraph(result.mml_bedgraph)
    me = _read_bedgraph(result.me_bedgraph)

    assert len(cov) == len(edges), f"expected {len(edges)} cov rows"
    assert (cov["chrom"] == "chr_test").all()
    assert list(zip(cov["start"], cov["end"], strict=False)) == edges
    assert (cov["value"] == N_READS).all(), \
        f"expected coverage {N_READS} per bin, got {cov['value'].tolist()}"

    assert len(mml) == len(edges)
    assert (mml["value"] == 1.0).all(), \
        f"expected MML=1.0 per bin, got {mml['value'].tolist()}"

    assert len(me) == len(edges)
    assert (me["value"] == 0.0).all(), \
        f"expected ME=0.0 per bin, got {me['value'].tolist()}"


def test_bam_fully_unmethylated_yields_mml_zero_and_zero_entropy(
    tmp_path: Path,
) -> None:
    """Fully unmethylated reads → every bin MML=0.0, ME=0.0, coverage=N_READS."""
    fa, bam = _build_uniform_bam(tmp_path / "fixture", methylated=False)
    out_dir = tmp_path / "out"
    result = _orchestrate(fa, bam=bam, tsv=None, out_dir=out_dir, label="unmeth")

    cov = _read_bedgraph(result.coverage_bedgraph)
    mml = _read_bedgraph(result.mml_bedgraph)
    me = _read_bedgraph(result.me_bedgraph)

    assert (cov["value"] == N_READS).all()
    assert (mml["value"] == 0.0).all()
    assert (me["value"] == 0.0).all()


def test_tsv_input_skips_modkit_and_reproduces_golden_values(
    tmp_path: Path,
) -> None:
    """Pre-computed TSV path: modkit runs once, then re-running with the
    TSV must produce the SAME bedgraphs as the BAM-input run."""
    fa, bam = _build_uniform_bam(tmp_path / "fixture", methylated=True)

    # First run: BAM → orchestrator generates a TSV alongside the bedgraphs.
    bam_dir = tmp_path / "out_bam"
    bam_result = _orchestrate(fa, bam=bam, tsv=None, out_dir=bam_dir,
                              label="from_bam")
    assert bam_result.tsv_path is not None
    assert bam_result.tsv_path.exists() and bam_result.tsv_path.stat().st_size > 0

    # Second run: feed the same TSV; expect identical entropy bedgraphs.
    tsv_dir = tmp_path / "out_tsv"
    tsv_result = _orchestrate(
        fa, bam=None, tsv=bam_result.tsv_path, out_dir=tsv_dir, label="from_tsv",
    )

    bam_cov = _read_bedgraph(bam_result.coverage_bedgraph)
    tsv_cov = _read_bedgraph(tsv_result.coverage_bedgraph)
    pd.testing.assert_frame_equal(bam_cov, tsv_cov)

    bam_mml = _read_bedgraph(bam_result.mml_bedgraph)
    tsv_mml = _read_bedgraph(tsv_result.mml_bedgraph)
    pd.testing.assert_frame_equal(bam_mml, tsv_mml)

    bam_me = _read_bedgraph(bam_result.me_bedgraph)
    tsv_me = _read_bedgraph(tsv_result.me_bedgraph)
    pd.testing.assert_frame_equal(bam_me, tsv_me)


def test_bam_partial_methylation_produces_nonzero_entropy(tmp_path: Path) -> None:
    """50% methylated reads → ME must be strictly between 0 and 1, MML
    near 0.5, all bins covered."""
    fa = write_reference_fasta(tmp_path / "ref.fa")
    bam = write_methylation_bam(
        tmp_path / "partial.bam", fa,
        n_reads=N_READS, methylated_fraction=0.5, seed=42,
    )
    out_dir = tmp_path / "out"
    result = _orchestrate(fa, bam=bam, tsv=None, out_dir=out_dir, label="partial")

    cov = _read_bedgraph(result.coverage_bedgraph)
    mml = _read_bedgraph(result.mml_bedgraph)
    me = _read_bedgraph(result.me_bedgraph)

    edges = _expected_bin_edges()
    assert len(cov) == len(edges)
    assert (cov["value"] == N_READS).all()

    # MML should hover around 0.5 (within ±0.4 of expectation given small N)
    assert ((mml["value"] >= 0.0) & (mml["value"] <= 1.0)).all()
    assert 0.1 <= mml["value"].mean() <= 0.9

    # ME must be > 0 because reads disagree, and ≤ 1 by construction.
    assert (me["value"] > 0.0).any(), "expected at least one bin with non-zero ME"
    assert ((me["value"] >= 0.0) & (me["value"] <= 1.0)).all()


def test_cli_run_on_bam_produces_same_bedgraphs_as_orchestrator(
    tmp_path: Path,
) -> None:
    """Sanity: invoking `shannonpore run` from Bash must produce identical
    bedgraphs to calling the Python orchestrator directly."""
    fa, bam = _build_uniform_bam(tmp_path / "fixture", methylated=True)

    py_dir = tmp_path / "py_run"
    py_result = _orchestrate(fa, bam=bam, tsv=None, out_dir=py_dir, label="py")

    cli_dir = tmp_path / "cli_run"
    proc = subprocess.run(
        [
            sys.executable, "-m", "src.cli", "run",
            "--bam", str(bam), "--label", "cli",
            "--fasta", str(fa),
            "--out-dir", str(cli_dir),
            "--threads", str(THREADS),
            "--mode", "true_mc",
            "--cpgs-per-bin", str(K),
            "--min-coverage", str(MIN_COV),
            "--chroms", "chr_test",
        ],
        cwd=str(V4_DIR), capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"CLI failed: {proc.stderr[-400:]}"

    cli_prefix = cli_dir / "cli_true_mc"
    for suffix in ("coverage", "me", "mml"):
        py_df = _read_bedgraph(getattr(py_result, f"{suffix}_bedgraph"))
        cli_df = _read_bedgraph(Path(f"{cli_prefix}.{suffix}.bedgraph"))
        pd.testing.assert_frame_equal(py_df, cli_df)
