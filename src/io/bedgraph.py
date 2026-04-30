"""Bedgraph reader/writer helpers."""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

from src.constants import BEDGRAPH_COLS
from src.io.utils_io import safe_mkdir


def read_bedgraph(path: str | Path, chrom_filter: str | None = None) -> pd.DataFrame:
    """Read a 4-column bedgraph (chrom, start, end, value).

    Returns an empty DataFrame on EOF/empty file rather than raising.
    """
    try:
        df = pd.read_csv(path, sep="\t", header=None, names=list(BEDGRAPH_COLS))
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=list(BEDGRAPH_COLS))

    if chrom_filter:
        df = df[df["chrom"] == chrom_filter]
    return df


def write_bedgraph(df: pd.DataFrame, value_col: str, out_path: str | Path) -> None:
    """Write a 4-column bedgraph: chrom, start, end, value_col.

    Empty `df` produces an empty file (consistent with v3 behaviour).
    """
    safe_mkdir(os.path.dirname(str(out_path)))
    if df.empty:
        Path(out_path).write_text("")
        return
    df2 = df[["chrom", "start", "end", value_col]].copy()
    df2.to_csv(out_path, sep="\t", header=False, index=False)


def smooth_bedgraph(
    path: str | Path,
    chrom_filter: str | None,
    window_size: int,
) -> pd.DataFrame:
    """Read a bedgraph and rolling-mean the value column, returning
    columns [mid, smoothed_score].
    """
    df = read_bedgraph(path, chrom_filter=chrom_filter)
    if df.empty:
        return pd.DataFrame(columns=["mid", "smoothed_score"])

    df = df.rename(columns={"value": "score"})
    df["mid"] = (df["start"] + df["end"]) / 2
    df = df.sort_values("mid")
    df["smoothed_score"] = (
        df["score"].rolling(window=window_size, center=True, min_periods=1).mean()
    )
    return df
