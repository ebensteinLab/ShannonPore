"""Stream GTF (optionally gzipped) into a DataFrame of gene/transcript
coordinates for plotting."""

from __future__ import annotations

import gzip
import logging
import re
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

            genes.append({
                "chr": chrom,
                "start": int(start_1based - 1),
                "end": int(end_1based),
                "target": str(target),
                "strand": strand,
            })

    if not genes:
        return pd.DataFrame(columns=["chr", "start", "end", "target", "strand"])
    return pd.DataFrame(genes)[["chr", "start", "end", "target", "strand"]]


def get_genes_for_region(
    gene_df: pd.DataFrame, chrom: str, xlim: tuple[int, int]
) -> pd.DataFrame:
    """Return genes whose intervals overlap `[xlim[0], xlim[1]]` on `chrom`."""
    cols = ["chr", "start", "end", "strand", "gene_name"]
    if gene_df is None or gene_df.empty:
        return pd.DataFrame(columns=cols)

    df = gene_df[
        (gene_df["chr"] == chrom)
        & (gene_df["start"] <= xlim[1])
        & (gene_df["end"] >= xlim[0])
    ].copy()
    if df.empty:
        return pd.DataFrame(columns=cols)
    df = df.rename(columns={"target": "gene_name"})
    return df[cols]
