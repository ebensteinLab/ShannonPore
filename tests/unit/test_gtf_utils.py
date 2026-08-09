"""Unit tests for GTF utilities."""

from __future__ import annotations

import pandas as pd
import pytest

from shannonpore.io.gtf_utils import (
    GeneStructure,
    find_gene_by_name,
    get_gene_structures_for_region,
    get_genes_for_region,
    load_gene_structures,
    load_genes_from_gtf,
    search_gene_names,
)


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


# ─── New API: GeneStructure / load_gene_structures / search ───────────────


@pytest.mark.unit
def test_load_gene_structures_returns_structures(tiny_gtf_with_exons) -> None:
    structs = load_gene_structures(str(tiny_gtf_with_exons))
    assert len(structs) == 2
    assert all(isinstance(g, GeneStructure) for g in structs)
    by_name = {g.name: g for g in structs}
    assert set(by_name) == {"GENE_A", "GENE_B"}
    # GENE_A on + strand: GTF 5000-5400 → 0-based start 4999, end 5400.
    ga = by_name["GENE_A"]
    assert ga.chrom == "chr1"
    assert ga.strand == "+"
    assert ga.start == 4999
    assert ga.end == 5400
    # Two non-overlapping exons → preserved as two separate entries.
    assert len(ga.exons) == 2


@pytest.mark.unit
def test_load_gene_structures_caches_result(tiny_gtf_with_exons) -> None:
    """Same path + promoter_upstream → same tuple object (lru_cache)."""
    a = load_gene_structures(str(tiny_gtf_with_exons))
    b = load_gene_structures(str(tiny_gtf_with_exons))
    assert a is b


@pytest.mark.unit
def test_promoter_respects_strand(tiny_gtf_with_exons) -> None:
    structs = load_gene_structures(str(tiny_gtf_with_exons), promoter_upstream=500)
    ga = next(g for g in structs if g.name == "GENE_A")  # + strand
    gb = next(g for g in structs if g.name == "GENE_B")  # − strand
    # + strand: promoter ends at TSS (gene start); width matches upstream.
    assert ga.promoter is not None
    ps, pe = ga.promoter
    assert pe == ga.start
    assert pe - ps == 500
    # − strand: promoter starts at gene end (the TSS for −).
    assert gb.promoter is not None
    qs, qe = gb.promoter
    assert qs == gb.end
    assert qe - qs == 500


@pytest.mark.unit
def test_promoter_upstream_changes_cache_key(tiny_gtf_with_exons) -> None:
    a = load_gene_structures(str(tiny_gtf_with_exons), promoter_upstream=1000)
    b = load_gene_structures(str(tiny_gtf_with_exons), promoter_upstream=500)
    # Different cache key → different tuple.
    assert a is not b


@pytest.mark.unit
def test_find_gene_by_name_case_insensitive(tiny_gtf_with_exons) -> None:
    structs = load_gene_structures(str(tiny_gtf_with_exons))
    assert find_gene_by_name(structs, "gene_a") is not None
    assert find_gene_by_name(structs, "GENE_A") is not None
    assert find_gene_by_name(structs, "  Gene_A  ") is not None


@pytest.mark.unit
def test_find_gene_by_name_returns_none_for_missing(tiny_gtf_with_exons) -> None:
    structs = load_gene_structures(str(tiny_gtf_with_exons))
    assert find_gene_by_name(structs, "NOT_A_GENE") is None
    assert find_gene_by_name(structs, "") is None
    assert find_gene_by_name(structs, "   ") is None


@pytest.mark.unit
def test_search_gene_names_prefix_then_substring(tiny_gtf_with_exons) -> None:
    structs = load_gene_structures(str(tiny_gtf_with_exons))
    out = search_gene_names(structs, "GENE")
    names = [g.name for g in out]
    assert names == ["GENE_A", "GENE_B"]  # both prefix-match


@pytest.mark.unit
def test_search_gene_names_empty_query(tiny_gtf_with_exons) -> None:
    structs = load_gene_structures(str(tiny_gtf_with_exons))
    assert search_gene_names(structs, "") == []
    assert search_gene_names(structs, "   ") == []


@pytest.mark.unit
def test_get_gene_structures_for_region_overlap(tiny_gtf_with_exons) -> None:
    structs = load_gene_structures(str(tiny_gtf_with_exons))
    # GENE_A spans 4999-5400; ask for [5100, 5300] overlap (inside the gene).
    matched = get_gene_structures_for_region(structs, "chr1", (5100, 5300))
    assert {g.name for g in matched} == {"GENE_A"}
    # Both genes when the window covers them.
    matched_all = get_gene_structures_for_region(structs, "chr1", (0, 7000))
    assert {g.name for g in matched_all} == {"GENE_A", "GENE_B"}
    # Wrong chrom → empty.
    assert get_gene_structures_for_region(structs, "chrX", (0, 7000)) == []
