"""Track plots: MML / ME / coverage / gene panels.

Lifted from v3 `src/plotting_tracks.py`. v4 uses `src.io.bedgraph.smooth_bedgraph`
to drop the duplicated rolling-mean code.
"""

from __future__ import annotations

import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.io.bedgraph import smooth_bedgraph
from src.io.gtf_utils import get_genes_for_region
from src.io.utils_io import safe_mkdir

logger = logging.getLogger(__name__)


def plot_two_bedgraph_overlays(
    top_files: list[str],
    bottom_files: list[str],
    gene_df: pd.DataFrame | None,
    chrom: str,
    start: int,
    end: int,
    out_path: str,
    window_size: int,
) -> None:
    """3-panel figure: top=MML overlay, middle=ME overlay, bottom=gene track."""
    from scipy.interpolate import interp1d
    from scipy.ndimage import label

    safe_mkdir(os.path.dirname(out_path))
    xlim = (int(start), int(end))
    region_label = f"{chrom}:{xlim[0]:,}-{xlim[1]:,}"

    fig, axes = plt.subplots(
        3, 1, sharex=True, figsize=(16, 9),
        gridspec_kw={"height_ratios": [4, 4, 1.5]},
    )
    ax_top, ax_bottom, ax_genes = axes

    dfs_top = [smooth_bedgraph(f, chrom_filter=chrom, window_size=window_size) for f in top_files]
    for df_, lbl in zip(dfs_top, ["Control MML", "Target MML"]):
        if not df_.empty:
            ax_top.plot(df_["mid"], df_["smoothed_score"], label=lbl, linewidth=2)
    ax_top.set_ylabel("MML", fontsize=14)
    ax_top.legend()
    ax_top.set_ylim([-0.05, 1.05])

    dfs_bottom = [smooth_bedgraph(f, chrom_filter=chrom, window_size=window_size) for f in bottom_files]
    for df_, lbl in zip(dfs_bottom, ["Control ME", "Target ME"]):
        if not df_.empty:
            ax_bottom.plot(df_["mid"], df_["smoothed_score"], label=lbl, linewidth=2)
    ax_bottom.set_ylabel("ME", fontsize=14)
    ax_bottom.legend()
    ax_bottom.set_ylim([-0.05, 1.05])

    if all(not d.empty for d in dfs_top) and all(not d.empty for d in dfs_bottom):
        x_grid = np.linspace(xlim[0], xlim[1], num=2000)
        f_top1 = interp1d(dfs_top[0]["mid"], dfs_top[0]["smoothed_score"],
                          bounds_error=False, fill_value="extrapolate")
        f_top2 = interp1d(dfs_top[1]["mid"], dfs_top[1]["smoothed_score"],
                          bounds_error=False, fill_value="extrapolate")
        f_bot1 = interp1d(dfs_bottom[0]["mid"], dfs_bottom[0]["smoothed_score"],
                          bounds_error=False, fill_value="extrapolate")
        f_bot2 = interp1d(dfs_bottom[1]["mid"], dfs_bottom[1]["smoothed_score"],
                          bounds_error=False, fill_value="extrapolate")
        top_div = np.abs(f_top1(x_grid) - f_top2(x_grid))
        bot_div = np.abs(f_bot1(x_grid) - f_bot2(x_grid))
        cond = (top_div < 0.1) & (bot_div > 0.2)

        labeled, _ = label(cond)
        for i in range(1, labeled.max() + 1):
            idx = np.where(labeled == i)[0]
            if idx.size > 0:
                ax_top.axvspan(x_grid[idx[0]], x_grid[idx[-1]],
                               color="lightgrey", alpha=0.5, zorder=0)
                ax_bottom.axvspan(x_grid[idx[0]], x_grid[idx[-1]],
                                  color="lightgrey", alpha=0.5, zorder=0)

    if gene_df is not None and hasattr(gene_df, "empty") and (not gene_df.empty):
        df_genes = get_genes_for_region(gene_df, chrom, xlim)
        if not df_genes.empty:
            gene_rows: list[list[dict]] = []
            for _, g in df_genes.sort_values("start").iterrows():
                placed = False
                for row in gene_rows:
                    if int(g["start"]) > int(row[-1]["end"]):
                        row.append(g.to_dict())
                        placed = True
                        break
                if not placed:
                    gene_rows.append([g.to_dict()])

            ax_genes.set_ylim(-0.5, len(gene_rows))
            for i, row in enumerate(gene_rows):
                for g in row:
                    color = "royalblue" if g["strand"] == "+" else "darkorange"
                    ax_genes.broken_barh(
                        [(int(g["start"]), int(g["end"]) - int(g["start"]))],
                        (len(gene_rows) - 1 - i, 0.8),
                        facecolors=color,
                    )
                    if (int(g["end"]) - int(g["start"])) > (xlim[1] - xlim[0]) * 0.05:
                        ax_genes.text(
                            (int(g["start"]) + int(g["end"])) / 2,
                            len(gene_rows) - 0.6 - i,
                            g["gene_name"], ha="center", va="bottom", fontsize=9,
                        )

    ax_genes.set_yticks([])
    ax_genes.set_ylabel("Genes", fontsize=14)
    ax_genes.set_xlabel(f"Genomic position ({region_label})", fontsize=14)

    fig.align_ylabels(axes)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
