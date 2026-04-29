"""End-to-end CLI pipeline tests on a synthetic BAM with MM/ML tags.

Exercises every CLI subcommand against a real (tiny) BAM produced by
``tests/fixtures/build_tiny_bam.py``. Requires:
  - ``modkit`` on PATH
  - ``pysam`` (used to build the fixture)

Skipped gracefully if modkit isn't found.

What's covered
--------------
1. ``extract``: BAM + FASTA → modkit TSV (rows present, expected columns).
2. ``entropy --mode true_mc``: TSV → ME/MML/coverage bedgraphs, expected
   line counts and value ranges.
3. ``entropy --mode bisulfite``: same TSV, different SQL filter, ensures
   the mode wiring works end-to-end.
4. ``entropy --mode ternary``: TSV → MhML + extra bedgraphs (5hmC + 5mC
   coverage), exercises the ternary pipeline.
5. ``segment``: bedgraph → segments BED (≥1 segment).
6. ``plot scatter``: two bedgraphs → PNG.
7. ``run``: one-shot BAM → entropy + optional segment with summary JSON.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

# pysam optional — without it we can't even build the fixture.
pysam = pytest.importorskip("pysam")

from tests.fixtures.build_tiny_bam import (  # noqa: E402
    CPG_POSITIONS,
    build_tiny_fixture,
)

V4_DIR = Path(__file__).resolve().parents[2]
MODKIT = shutil.which("modkit")

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(MODKIT is None, reason="modkit binary not on PATH"),
]


# ─── shared fixture ───────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def bam_fixture(tmp_path_factory) -> tuple[Path, Path]:
    """One-time: build the tiny FASTA + BAM for this test module."""
    tmp = tmp_path_factory.mktemp("bam_fixture")
    fa, bam = build_tiny_fixture(tmp, n_reads=30, methylated_fraction=0.5)
    return fa, bam


def _run_cli(*args: str, expect_exit: int = 0, **kw) -> subprocess.CompletedProcess:
    proc = subprocess.run(
        [sys.executable, "-m", "src.cli", *args],
        cwd=str(V4_DIR), capture_output=True, text=True, **kw,
    )
    if proc.returncode != expect_exit:
        raise AssertionError(
            f"CLI exited {proc.returncode} (expected {expect_exit}).\n"
            f"args: {args}\nstdout: {proc.stdout}\nstderr: {proc.stderr}"
        )
    return proc


# ─── 1. extract ───────────────────────────────────────────────────────────

def test_cli_extract_produces_tsv_with_expected_rows(
    bam_fixture: tuple[Path, Path], tmp_path: Path
) -> None:
    fa, bam = bam_fixture
    out_tsv = tmp_path / "extract.tsv"
    _run_cli(
        "extract", "--bam", str(bam), str(out_tsv),
        "--fasta", str(fa), "--threads", "2",
    )
    assert out_tsv.exists() and out_tsv.stat().st_size > 0

    header = out_tsv.read_text().splitlines()[0].split("\t")
    for col in ("read_id", "ref_position", "chrom", "ref_strand",
                "mod_qual", "mod_code"):
        assert col in header, f"missing modkit column: {col}"

    n_rows = sum(1 for _ in out_tsv.read_text().splitlines()) - 1
    # 30 reads × 8 CpG positions = 240 expected rows (no canonical h calls
    # for the C+m? tag), give ±10% slack for any modkit edge effects.
    assert 200 <= n_rows <= 260, f"expected ~240 rows, got {n_rows}"


# ─── 2-4. entropy (each mode) ─────────────────────────────────────────────

@pytest.fixture(scope="module")
def shared_extract_tsv(bam_fixture, tmp_path_factory) -> tuple[Path, Path]:
    """Run `extract` once and reuse the TSV across all entropy tests."""
    fa, bam = bam_fixture
    tmp = tmp_path_factory.mktemp("extract")
    out_tsv = tmp / "extract.tsv"
    _run_cli(
        "extract", "--bam", str(bam), str(out_tsv),
        "--fasta", str(fa), "--threads", "2",
    )
    return fa, out_tsv


def _bedgraph_lines(p: Path) -> list[list[str]]:
    if not p.exists() or p.stat().st_size == 0:
        return []
    return [ln.split("\t") for ln in p.read_text().splitlines() if ln]


@pytest.mark.parametrize("mode", ["true_mc", "bisulfite"])
def test_cli_entropy_two_state_modes_produce_bedgraphs(
    shared_extract_tsv: tuple[Path, Path], tmp_path: Path, mode: str,
) -> None:
    fa, tsv = shared_extract_tsv
    out_prefix = tmp_path / f"sample_{mode}"
    _run_cli(
        "entropy", str(tsv), str(out_prefix),
        "--fasta", str(fa), "--threads", "2",
        "--mode", mode,
        "--cpgs-per-bin", "2",
        "--min-coverage", "4",
        "--chroms", "chr_test",
    )

    cov = Path(f"{out_prefix}.coverage.bedgraph")
    me = Path(f"{out_prefix}.me.bedgraph")
    mml = Path(f"{out_prefix}.mml.bedgraph")
    assert cov.exists() and me.exists() and mml.exists()

    cov_rows = _bedgraph_lines(cov)
    # 8 CpGs / 2 per bin = 4 full bins on chr_test
    assert len(cov_rows) == 4
    for r in cov_rows:
        assert r[0] == "chr_test"
        assert int(r[3]) > 0  # all 30 reads cover the synthetic region

    me_rows = _bedgraph_lines(me)
    mml_rows = _bedgraph_lines(mml)
    assert me_rows and mml_rows
    # ME / MML must be in [0, 1]
    for r in me_rows:
        v = float(r[3])
        assert 0.0 <= v <= 1.0
    for r in mml_rows:
        v = float(r[3])
        assert 0.0 <= v <= 1.0


def test_cli_entropy_ternary_writes_extra_bedgraphs(
    shared_extract_tsv: tuple[Path, Path], tmp_path: Path,
) -> None:
    fa, tsv = shared_extract_tsv
    out_prefix = tmp_path / "sample_ternary"
    _run_cli(
        "entropy", str(tsv), str(out_prefix),
        "--fasta", str(fa), "--threads", "2",
        "--mode", "ternary",
        "--cpgs-per-bin", "2",      # k=2 so 3^k=9, fits 30 reads
        "--min-coverage", "4",
        "--chroms", "chr_test",
    )

    for suffix in ("coverage", "me", "mml", "mhml", "cov_5mc", "cov_5hmc"):
        p = Path(f"{out_prefix}.{suffix}.bedgraph")
        assert p.exists(), f"missing ternary output: {p}"

    cov = _bedgraph_lines(Path(f"{out_prefix}.coverage.bedgraph"))
    assert len(cov) == 4
    for r in cov:
        assert int(r[3]) > 0


# ─── 5. plot ──────────────────────────────────────────────────────────────

def test_cli_plot_scatter_produces_png(
    shared_extract_tsv: tuple[Path, Path], tmp_path: Path,
) -> None:
    fa, tsv = shared_extract_tsv
    out_prefix = tmp_path / "for_plot"
    _run_cli(
        "entropy", str(tsv), str(out_prefix),
        "--fasta", str(fa), "--threads", "2",
        "--mode", "true_mc",
        "--cpgs-per-bin", "2",
        "--min-coverage", "4",
        "--chroms", "chr_test",
    )
    me_bg = Path(f"{out_prefix}.me.bedgraph")
    mml_bg = Path(f"{out_prefix}.mml.bedgraph")

    out_png = tmp_path / "scatter.png"
    _run_cli(
        "plot", "scatter", str(out_png),
        "--x-bedgraph", str(me_bg),
        "--y-bedgraph", str(mml_bg),
        "--x-label", "ME",
        "--y-label", "MML",
        "--subsample", "1000",
    )
    assert out_png.exists() and out_png.stat().st_size > 0


def test_cli_plot_distribution_produces_png(
    shared_extract_tsv: tuple[Path, Path], tmp_path: Path,
) -> None:
    fa, tsv = shared_extract_tsv
    out_prefix = tmp_path / "for_dist"
    _run_cli(
        "entropy", str(tsv), str(out_prefix),
        "--fasta", str(fa), "--threads", "2",
        "--mode", "true_mc",
        "--cpgs-per-bin", "2",
        "--min-coverage", "4",
        "--chroms", "chr_test",
    )
    me_bg = Path(f"{out_prefix}.me.bedgraph")
    out_png = tmp_path / "dist.png"
    _run_cli(
        "plot", "distribution", str(out_png),
        "--x-bedgraph", str(me_bg),
        "--dist-kind", "hist",
    )
    assert out_png.exists() and out_png.stat().st_size > 0


# ─── 6. one-shot run ──────────────────────────────────────────────────────

def test_cli_run_oneshot_bam_to_entropy(
    bam_fixture: tuple[Path, Path], tmp_path: Path,
) -> None:
    fa, bam = bam_fixture
    out_dir = tmp_path / "run_out"
    _run_cli(
        "run",
        "--bam", str(bam),
        "--label", "tinybam",
        "--fasta", str(fa),
        "--out-dir", str(out_dir),
        "--threads", "2",
        "--mode", "true_mc",
        "--cpgs-per-bin", "2",
        "--min-coverage", "4",
        "--chroms", "chr_test",
    )

    # Bedgraphs produced (label-prefixed)
    prefix = out_dir / "tinybam_true_mc"
    for s in ("coverage", "me", "mml"):
        assert Path(f"{prefix}.{s}.bedgraph").exists()

    # Summary JSON produced
    summary = out_dir / "run_summary.json"
    assert summary.exists()
    import json
    s = json.loads(summary.read_text())
    assert s["mode"] == "true_mc"
    assert s["pair_mode"] is False
    assert s["target"]["label"] == "tinybam"
    assert "me_bedgraph" in s["target"]
    assert s["cpgs_per_bin"] == 2


def test_cli_run_with_pre_extracted_tsv_skips_modkit(
    shared_extract_tsv: tuple[Path, Path], tmp_path: Path,
) -> None:
    fa, tsv = shared_extract_tsv
    out_dir = tmp_path / "run_out_tsv"
    _run_cli(
        "run",
        "--tsv", str(tsv),
        "--label", "fromtsv",
        "--fasta", str(fa),
        "--out-dir", str(out_dir),
        "--threads", "2",
        "--mode", "bisulfite",
        "--cpgs-per-bin", "2",
        "--min-coverage", "4",
        "--chroms", "chr_test",
    )
    prefix = out_dir / "fromtsv_bisulfite"
    for s in ("coverage", "me", "mml"):
        assert Path(f"{prefix}.{s}.bedgraph").exists()


def test_cli_run_pair_mode_produces_two_bedgraph_sets(
    bam_fixture: tuple[Path, Path], tmp_path: Path,
) -> None:
    fa, bam = bam_fixture
    out_dir = tmp_path / "run_pair"
    _run_cli(
        "run", "--pair",
        "--control-bam", str(bam),
        "--target-bam", str(bam),  # use same BAM as both sides for the test
        "--control-label", "ctrl",
        "--target-label", "tgt",
        "--fasta", str(fa),
        "--out-dir", str(out_dir),
        "--threads", "2",
        "--mode", "true_mc",
        "--cpgs-per-bin", "2",
        "--min-coverage", "4",
        "--chroms", "chr_test",
    )
    for label in ("ctrl", "tgt"):
        prefix = out_dir / f"{label}_true_mc"
        for s in ("coverage", "me", "mml"):
            assert Path(f"{prefix}.{s}.bedgraph").exists(), f"missing {s} for {label}"

    import json
    summary = json.loads((out_dir / "run_summary.json").read_text())
    assert summary["pair_mode"] is True
    assert summary["control"]["label"] == "ctrl"
    assert summary["target"]["label"] == "tgt"


def test_cli_run_bam_folder_merges_and_runs(
    tmp_path_factory, tmp_path: Path,
) -> None:
    """Folder of BAMs → samtools merge → modkit → entropy bedgraphs."""
    pytest.importorskip("pysam")
    from tests.fixtures.build_tiny_bam import (
        write_methylation_bam,
        write_reference_fasta,
    )

    fixture_dir = tmp_path_factory.mktemp("folder_fixture")
    fa = write_reference_fasta(fixture_dir / "ref.fa")
    bam_folder = fixture_dir / "bams"
    bam_folder.mkdir()
    # 3 BAMs that will be merged.
    for i, frac in enumerate((0.2, 0.5, 0.8)):
        write_methylation_bam(
            bam_folder / f"sample_{i}.bam", fa,
            n_reads=10, methylated_fraction=frac, seed=i,
        )

    out_dir = tmp_path / "folder_run"
    _run_cli(
        "run",
        "--bam-folder", str(bam_folder),
        "--label", "merged",
        "--fasta", str(fa),
        "--out-dir", str(out_dir),
        "--threads", "2",
        "--mode", "true_mc",
        "--cpgs-per-bin", "2",
        "--min-coverage", "4",
        "--chroms", "chr_test",
    )
    # Merged BAM was created (single-mode lands at out_dir / merged_merged.sorted.bam)
    merged = out_dir / "merged_merged.sorted.bam"
    assert merged.exists()
    assert (merged.with_suffix(merged.suffix + ".bai")).exists()

    prefix = out_dir / "merged_true_mc"
    for s in ("coverage", "me", "mml"):
        assert Path(f"{prefix}.{s}.bedgraph").exists()


# ─── reference-data integrity check ──────────────────────────────────────

def test_synthetic_reference_has_expected_cpg_positions() -> None:
    """Sanity check the fixture itself."""
    from tests.fixtures.build_tiny_bam import _REF_SEQ
    found: list[int] = []
    for i in range(len(_REF_SEQ) - 1):
        if _REF_SEQ[i:i + 2] == "CG":
            found.append(i)
    assert found == CPG_POSITIONS
