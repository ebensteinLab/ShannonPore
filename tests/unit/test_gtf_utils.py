"""Unit tests for GTF utilities."""

from __future__ import annotations

import pandas as pd
import pytest

from src.io.gtf_utils import get_genes_for_region, load_genes_from_gtf


@pytest.mark.unit
def test_load_two_genes(tiny_gtf) -> None:
    df = load_genes_from_gtf(str(tiny_gtf))
    assert len(df) == 2
    assert set(df["target"]) == {"GENE_A", "GENE_B"}
    # Column 1 is ensured to be prefixed with "chr"
    assert all(c.startswith("chr") for c in df["chr"])


@pytest.mark.unit
def test_gtf_start_is_zero_based(tiny_gtf) -> None:
    df = load_genes_from_gtf(str(tiny_gtf))
    # GTF is 1-based; loader converts to 0-based start
    row = df[df["target"] == "GENE_A"].iloc[0]
    assert int(row["start"]) == 99  # 100 (1-based) → 99 (0-based)
    assert int(row["end"]) == 200


@pytest.mark.unit
def test_get_genes_for_region_filters_by_overlap(tiny_gtf) -> None:
    df = load_genes_from_gtf(str(tiny_gtf))
    matched = get_genes_for_region(df, "chr1", (150, 350))
    # GENE_A (99-200) and GENE_B (299-500) both overlap [150, 350]
    assert set(matched["gene_name"]) == {"GENE_A", "GENE_B"}


@pytest.mark.unit
def test_get_genes_for_region_empty_input() -> None:
    df = pd.DataFrame(columns=["chr", "start", "end", "target", "strand"])
    out = get_genes_for_region(df, "chr1", (0, 1000))
    assert out.empty


@pytest.mark.unit
def test_load_empty_gtf_returns_empty_df(tmp_path) -> None:
    p = tmp_path / "empty.gtf"
    p.write_text("# comment only\n")
    df = load_genes_from_gtf(str(p))
    assert df.empty
