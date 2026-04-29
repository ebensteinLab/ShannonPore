"""pytest configuration: ensure v4 root is importable as `src.*`."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    d = ROOT / "tests" / "fixtures"
    d.mkdir(parents=True, exist_ok=True)
    return d


@pytest.fixture()
def tiny_bedgraph(tmp_path: Path) -> Path:
    """A 4-column bedgraph with deterministic values."""
    p = tmp_path / "tiny.bedgraph"
    p.write_text(
        "chr1\t10\t20\t0.10\n"
        "chr1\t20\t30\t0.20\n"
        "chr1\t30\t40\t0.30\n"
        "chr1\t40\t50\t0.40\n"
        "chr2\t10\t20\t0.50\n"
    )
    return p


@pytest.fixture()
def tiny_gtf(tmp_path: Path) -> Path:
    """A minimal GTF with 2 genes."""
    p = tmp_path / "tiny.gtf"
    p.write_text(
        '1\thavana\tgene\t100\t200\t.\t+\t.\tgene_id "ENSG001"; gene_name "GENE_A";\n'
        '1\thavana\tgene\t300\t500\t.\t-\t.\tgene_id "ENSG002"; gene_name "GENE_B";\n'
    )
    return p


@pytest.fixture()
def tiny_roi_bed(tmp_path: Path) -> Path:
    p = tmp_path / "tiny_roi.bed"
    p.write_text("chr1\t100\t200\tregion_A\nchr2\t50\t150\tregion_B\n")
    return p
