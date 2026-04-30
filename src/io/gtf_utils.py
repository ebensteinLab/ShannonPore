"""Stream GTF (optionally gzipped) into a DataFrame of gene/transcript
coordinates for plotting."""

from __future__ import annotations

import functools
import gzip
import logging
import re
from dataclasses import dataclass, field
from typing import IO

import pandas as pd

logger = logging.getLogger(__name__)


def _open_text_maybe_gz(path: str) -> IO[str]:
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace", encoding="utf-8")


_ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"\s*;')


def _parse_gtf_attributes(attr_str: str) -> dict[str, str]:
    return {m.group(1): m.group(2) for m in _ATTR_RE.finditer(attr_str)}


def load_genes_from_gtf(gtf_path: str) -> pd.DataFrame:
    """Stream a GTF (optionally .gz) and return columns
    [chr, start (0-based), end (1-based), target (gene name), strand].
    """
    if not gtf_path:
        return pd.DataFrame(columns=["chr", "start", "end", "target", "strand"])

    genes: list[dict[str, object]] = []
    seen: set[tuple] = set()

    with _open_text_maybe_gz(gtf_path) as f:
        for ln in f:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            chrom, _src, feature, start_s, end_s, _score, strand, _frame, attrs = parts[:9]

            if not str(chrom).startswith("chr"):
                chrom = "chr" + str(chrom)
            if feature not in ("gene", "transcript"):
                continue

            try:
                start_1based = int(start_s)
                end_1based = int(end_s)
            except ValueError:
                continue
            if end_1based <= start_1based:
                continue

            attr_map = _parse_gtf_attributes(attrs)
            target = attr_map.get("gene_name") or attr_map.get("gene_id")
            if not target:
                continue

            key = (chrom, target, start_1based, end_1based)
            if key in seen:
                continue
            seen.add(key)

            genes.append(
                {
                    "chr": chrom,
                    "start": int(start_1based - 1),
                    "end": int(end_1based),
                    "target": str(target),
                    "strand": strand,
                }
            )

    if not genes:
        return pd.DataFrame(columns=["chr", "start", "end", "target", "strand"])
    return pd.DataFrame(genes)[["chr", "start", "end", "target", "strand"]]


# ─── Gene structure (exons + promoter) for richer track plots ─────────────


@dataclass(frozen=True)
class GeneStructure:
    """One gene's plotted geometry: span, strand, exon list, promoter."""

    chrom: str
    name: str
    start: int  # 0-based inclusive
    end: int  # 1-based exclusive (so end - start = length)
    strand: str
    exons: tuple[tuple[int, int], ...] = field(default_factory=tuple)
    promoter: tuple[int, int] | None = None  # (start, end), 0-based


def _promoter_for(start: int, end: int, strand: str, upstream: int) -> tuple[int, int]:
    """1-kb (default) upstream window relative to TSS."""
    if strand == "-":
        return (end, end + upstream)
    return (max(0, start - upstream), start)


@functools.lru_cache(maxsize=4)
def _load_gene_structures_cached(
    gtf_path: str, promoter_upstream: int
) -> tuple[GeneStructure, ...]:
    """Stream GTF once, group exons by gene, attach promoter window."""
    if not gtf_path:
        return ()

    gene_meta: dict[str, dict] = {}  # name -> {chrom, start, end, strand}
    gene_exons: dict[str, list[tuple[int, int]]] = {}

    with _open_text_maybe_gz(gtf_path) as f:
        for ln in f:
            if not ln or ln.startswith("#"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _src, feature, start_s, end_s, _score, strand, _frame, attrs = parts[:9]
            if feature not in ("gene", "exon"):
                continue
            if not str(chrom).startswith("chr"):
                chrom = "chr" + str(chrom)
            try:
                s_1 = int(start_s)
                e_1 = int(end_s)
            except ValueError:
                continue
            if e_1 <= s_1:
                continue
            attr_map = _parse_gtf_attributes(attrs)
            name = attr_map.get("gene_name") or attr_map.get("gene_id")
            if not name:
                continue

            s_0 = s_1 - 1  # convert to 0-based
            if feature == "gene":
                meta = gene_meta.setdefault(
                    name,
                    {"chrom": chrom, "start": s_0, "end": e_1, "strand": strand},
                )
                # If duplicates exist (e.g. RefSeq paralogs), widen the span.
                meta["start"] = min(meta["start"], s_0)
                meta["end"] = max(meta["end"], e_1)
            else:  # exon
                gene_exons.setdefault(name, []).append((s_0, e_1))
                # Some GTFs lack 'gene' rows — synthesise from exon spans.
                meta = gene_meta.setdefault(
                    name,
                    {"chrom": chrom, "start": s_0, "end": e_1, "strand": strand},
                )
                meta["start"] = min(meta["start"], s_0)
                meta["end"] = max(meta["end"], e_1)

    out: list[GeneStructure] = []
    for name, meta in gene_meta.items():
        # Merge overlapping/touching exons (RefSeq lists per-transcript exons).
        ex = sorted(gene_exons.get(name, []))
        merged: list[tuple[int, int]] = []
        for s, e in ex:
            if merged and s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        prom = _promoter_for(
            meta["start"],
            meta["end"],
            meta["strand"],
            promoter_upstream,
        )
        out.append(
            GeneStructure(
                chrom=meta["chrom"],
                name=name,
                start=meta["start"],
                end=meta["end"],
                strand=meta["strand"],
                exons=tuple(merged),
                promoter=prom,
            )
        )
    return tuple(out)


def load_gene_structures(
    gtf_path: str,
    *,
    promoter_upstream: int = 1000,
) -> tuple[GeneStructure, ...]:
    """Public wrapper that caches the parsed GTF in memory.

    Subsequent calls with the same path are O(1). Pass an empty path to
    get an empty tuple back.
    """
    return _load_gene_structures_cached(gtf_path or "", int(promoter_upstream))


def find_gene_by_name(
    structures: tuple[GeneStructure, ...],
    name: str,
) -> GeneStructure | None:
    """Case-insensitive exact-name lookup. Returns the first hit, or None."""
    if not name:
        return None
    target = name.strip().upper()
    if not target:
        return None
    for g in structures:
        if g.name.upper() == target:
            return g
    return None


def search_gene_names(
    structures: tuple[GeneStructure, ...],
    query: str,
    *,
    limit: int = 20,
) -> list[GeneStructure]:
    """Case-insensitive prefix / substring search for autocomplete.

    Exact matches rank first, then prefix matches, then substring.
    """
    if not query:
        return []
    q = query.strip().upper()
    if not q:
        return []
    exact: list[GeneStructure] = []
    prefix: list[GeneStructure] = []
    substr: list[GeneStructure] = []
    for g in structures:
        n = g.name.upper()
        if n == q:
            exact.append(g)
        elif n.startswith(q):
            prefix.append(g)
        elif q in n:
            substr.append(g)
        if len(exact) + len(prefix) + len(substr) >= limit * 3:
            break
    return (exact + prefix + substr)[:limit]


def get_gene_structures_for_region(
    structures: tuple[GeneStructure, ...],
    chrom: str,
    xlim: tuple[int, int],
) -> list[GeneStructure]:
    """Filter pre-parsed structures to those overlapping `chrom:xlim`."""
    if not structures:
        return []
    lo, hi = int(xlim[0]), int(xlim[1])
    return [g for g in structures if g.chrom == chrom and g.start <= hi and g.end >= lo]


def get_genes_for_region(gene_df: pd.DataFrame, chrom: str, xlim: tuple[int, int]) -> pd.DataFrame:
    """Return genes whose intervals overlap `[xlim[0], xlim[1]]` on `chrom`."""
    cols = ["chr", "start", "end", "strand", "gene_name"]
    if gene_df is None or gene_df.empty:
        return pd.DataFrame(columns=cols)

    df = gene_df[
        (gene_df["chr"] == chrom) & (gene_df["start"] <= xlim[1]) & (gene_df["end"] >= xlim[0])
    ].copy()
    if df.empty:
        return pd.DataFrame(columns=cols)
    df = df.rename(columns={"target": "gene_name"})
    return df[cols]
