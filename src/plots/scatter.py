"""Scatter / hexbin plots for control-vs-target comparisons.

Migrated from v3 `_show_scatter` block (lines 6060-7186 of app.py, ~1,126
lines). v4 splits by plot kind (hexbin, density-coloured, paired-overlay).
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def hexbin_density(
    df_x: pd.DataFrame,
    df_y: pd.DataFrame,
    *,
    x_label: str,
    y_label: str,
    title: str = "",
    out_path: Path | str | None = None,
    gridsize: int = 60,
    bins: str = "log",
) -> plt.Figure:
    """Hexbin scatter for two value series. df_x/df_y must share the same
    bedgraph row order (chrom/start/end). Misaligned series are dropped.
    """
    common = df_x.merge(df_y, on=["chrom", "start", "end"], suffixes=("_x", "_y"))
    if common.empty:
        logger.warning("hexbin_density: no overlapping rows; returning empty fig")
        fig, _ = plt.subplots()
        return fig

    fig, ax = plt.subplots(figsize=(7, 6))
    hb = ax.hexbin(
        common["value_x"], common["value_y"],
        gridsize=gridsize, bins=bins, cmap="viridis", mincnt=1,
    )
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if title:
        ax.set_title(title)
    fig.colorbar(hb, ax=ax, label="log10(N)" if bins == "log" else "N")
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=200)
    return fig


def paired_scatter(
    df_x: pd.DataFrame,
    df_y: pd.DataFrame,
    *,
    x_label: str,
    y_label: str,
    point_color: str = "#1f77b4",
    alpha: float = 0.4,
    subsample: int | None = 50_000,
    out_path: Path | str | None = None,
) -> plt.Figure:
    """Plain scatter with optional subsampling for performance."""
    common = df_x.merge(df_y, on=["chrom", "start", "end"], suffixes=("_x", "_y"))
    if subsample is not None and len(common) > subsample:
        common = common.sample(subsample, random_state=42)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        common["value_x"], common["value_y"],
        s=6, alpha=alpha, color=point_color,
    )
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=1)
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=200)
    return fig


def region_subset_scatter(
    df_x: pd.DataFrame,
    df_y: pd.DataFrame,
    region_bed: pd.DataFrame,
    *,
    x_label: str,
    y_label: str,
    out_path: Path | str | None = None,
) -> plt.Figure:
    """Scatter restricted to rows overlapping any region in `region_bed`
    (columns: chr, start, end)."""
    keep = []
    for _, r in region_bed.iterrows():
        mask = (
            (df_x["chrom"] == r["chr"])
            & (df_x["end"] > r["start"])
            & (df_x["start"] < r["end"])
        )
        keep.append(df_x[mask])
    if not keep:
        return plt.subplots(figsize=(7, 6))[0]
    sub_x = pd.concat(keep)
    sub_y = df_y[df_y[["chrom", "start", "end"]].apply(tuple, axis=1).isin(
        sub_x[["chrom", "start", "end"]].apply(tuple, axis=1)
    )]
    return paired_scatter(
        sub_x, sub_y, x_label=x_label, y_label=y_label,
        point_color="#ff7f0e", alpha=0.6, subsample=None, out_path=out_path,
    )


def ternary_prevalence_panel(
    prevalence_df: pd.DataFrame,
    out_path: Path | str | None = None,
) -> plt.Figure:
    """Stacked-bar panel of C / 5mC / 5hmC prevalence per sample, faceted
    by Dataset. Mirrors the user's selection in app.py at v3.
    """
    datasets = prevalence_df["Dataset"].unique().tolist()
    fig, axes = plt.subplots(1, max(len(datasets), 1), figsize=(5 * len(datasets), 5))
    if len(datasets) == 1:
        axes = [axes]

    for ax, ds in zip(axes, datasets):
        sub = prevalence_df[prevalence_df["Dataset"] == ds].reset_index(drop=True)
        if sub.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(ds)
            continue
        x = np.arange(len(sub))
        ax.bar(x, sub["5hmC_pct"], color="#27ae60", label="5hmC", edgecolor="white")
        ax.bar(x, sub["5mC_pct"], bottom=sub["5hmC_pct"],
               color="#3498db", label="5mC", edgecolor="white")
        ax.bar(
            x, sub["C_pct"],
            bottom=sub["5hmC_pct"] + sub["5mC_pct"],
            color="#bdc3c7", label="Unmodified C", edgecolor="white",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(sub["Sample"], rotation=15, ha="right", fontsize=9)
        ax.set_ylim(0, 100)
        ax.set_ylabel("% of CpG calls")
        ax.set_title(ds, fontsize=12, fontweight="bold")
        ax.legend(loc="upper right", fontsize=8)

    fig.suptitle(
        "C / 5mC / 5hmC Prevalence per Sample (from ternary CpG calls)",
        fontsize=14, fontweight="bold",
    )
    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=200)
    return fig
