"""Genomic-region track plot: gene structure on top, ME and MML below.

The gene panel renders each overlapping gene with:
  * a thin intron line + small direction arrows that follow the strand,
  * thick rectangles for exons,
  * a red rectangle for the 1 kb upstream promoter,
  * the gene name annotated below.

The signal panels are smoothed with a uniform filter (window in bins)
and shaded under the line so the two samples are easy to compare.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

from src.io.bedgraph import read_bedgraph
from src.io.gtf_utils import (
    GeneStructure,
    get_gene_structures_for_region,
    load_gene_structures,
)
from src.io.utils_io import safe_mkdir

logger = logging.getLogger(__name__)


# ─── Gene structure panel ─────────────────────────────────────────────────

def _draw_gene(ax, gene: GeneStructure, *, y: float, plot_lo: int, plot_hi: int) -> None:
    """Draw one gene at vertical position ``y`` on ``ax``."""
    # Intron line spanning the whole gene
    ax.plot(
        [gene.start, gene.end], [y, y],
        color="#333333", lw=1.5, zorder=1,
    )
    # Direction arrows along the intron line
    span = gene.end - gene.start
    if span > 0:
        n_arrows = max(2, int(min(20, max(2, span // 800))))
        positions = np.linspace(gene.start, gene.end, n_arrows + 2)[1:-1]
        dx = max(50, span // 200) * (1 if gene.strand == "+" else -1)
        for ap in positions:
            ax.annotate(
                "",
                xy=(ap + dx, y), xytext=(ap, y),
                arrowprops={
                    "arrowstyle": "->",
                    "color": "#666666",
                    "lw": 0.8,
                },
            )
    # Exons
    for es, ee in gene.exons:
        ax.add_patch(plt.Rectangle(
            (es, y - 0.30), ee - es, 0.60,
            facecolor="#2c3e50", edgecolor="black", lw=0.8, zorder=3,
        ))
    # Promoter
    if gene.promoter is not None:
        ps, pe = gene.promoter
        ax.add_patch(plt.Rectangle(
            (ps, y - 0.30), pe - ps, 0.60,
            facecolor="#e74c3c", edgecolor="#c0392b", lw=0.8, zorder=3,
        ))
        # Promoter label only if it's actually visible in the window.
        if ps <= plot_hi and pe >= plot_lo:
            ax.text(
                (max(ps, plot_lo) + min(pe, plot_hi)) / 2, y + 0.55,
                "Promoter",
                ha="center", va="bottom",
                fontsize=8, color="#e74c3c", fontstyle="italic",
            )
    # Gene name
    ax.text(
        (gene.start + gene.end) / 2, y - 0.65,
        gene.name, ha="center", va="top",
        fontsize=11, fontweight="bold",
    )


def _stack_genes_by_row(genes: list[GeneStructure]) -> list[list[GeneStructure]]:
    """Pack overlapping genes onto separate rows (greedy left-to-right)."""
    rows: list[list[GeneStructure]] = []
    for g in sorted(genes, key=lambda x: x.start):
        placed = False
        for row in rows:
            if g.start > row[-1].end:
                row.append(g)
                placed = True
                break
        if not placed:
            rows.append([g])
    return rows


def _draw_gene_panel(
    ax, genes: list[GeneStructure], plot_lo: int, plot_hi: int,
) -> None:
    if not genes:
        ax.text(
            0.5, 0.5, "No genes in this window",
            transform=ax.transAxes,
            ha="center", va="center", fontsize=10, color="#888888",
        )
        ax.set_xlim(plot_lo, plot_hi)
        ax.set_ylim(-0.5, 1.5)
        ax.axis("off")
        return

    rows = _stack_genes_by_row(genes)
    n = len(rows)
    # Reserve 1.5 vertical units per row so the labels don't collide.
    for i, row in enumerate(rows):
        y = (n - i) * 1.5  # top row gets the highest y
        for g in row:
            _draw_gene(ax, g, y=y, plot_lo=plot_lo, plot_hi=plot_hi)
    ax.set_xlim(plot_lo, plot_hi)
    ax.set_ylim(0, n * 1.5 + 1.0)
    ax.axis("off")


# ─── Signal (ME / MML) panel ──────────────────────────────────────────────

def _read_region(
    path: str | Path | None, chrom: str, lo: int, hi: int,
) -> pd.DataFrame:
    if not path:
        return pd.DataFrame(columns=["chrom", "start", "end", "value"])
    df = read_bedgraph(str(path))
    if df.empty:
        return df
    # Half-open overlap so a bin straddling the window edge is kept rather
    # than silently dropped (BED bins can be wider than the padding).
    return df[
        (df["chrom"] == chrom)
        & (df["start"] < hi)
        & (df["end"] > lo)
    ].copy()


def _smooth(vals: np.ndarray, win: int) -> np.ndarray:
    if win > 1 and vals.size > win:
        return uniform_filter1d(vals, size=int(win))
    return vals


def _plot_signal(
    ax, *, ctrl_df: pd.DataFrame, case_df: pd.DataFrame,
    label_a: str, label_b: str, color_a: str, color_b: str,
    ylabel: str, plot_lo: int, plot_hi: int, smooth_win: int,
    ylim: tuple[float, float] | None = (0.0, 1.05),
) -> None:
    drew_anything = False
    y_max_seen = 0.0
    for df, label, color in (
        (ctrl_df, label_a, color_a),
        (case_df, label_b, color_b),
    ):
        if df.empty:
            continue
        df = df.sort_values("start")
        x = ((df["start"].to_numpy() + df["end"].to_numpy()) / 2.0)
        y = _smooth(df["value"].to_numpy(dtype=float), smooth_win)
        ax.plot(x, y, color=color, lw=1.8, alpha=0.9, label=label)
        ax.fill_between(x, y, alpha=0.15, color=color)
        drew_anything = True
        y_max_seen = max(y_max_seen, float(np.nanmax(y)) if y.size else 0.0)

    if not drew_anything:
        ax.text(
            0.5, 0.5, "No bins in region",
            transform=ax.transAxes,
            ha="center", va="center", fontsize=10, color="#888888",
        )
    ax.set_xlim(plot_lo, plot_hi)
    if ylim is not None:
        ax.set_ylim(*ylim)
    elif drew_anything:
        # 5 % headroom so the line doesn't touch the top of the panel.
        ax.set_ylim(0, max(1.0, y_max_seen) * 1.05)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(axis="y", alpha=0.2)
    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{v / 1e6:.3f} Mb")
    )


# ─── Main entry point ─────────────────────────────────────────────────────

def plot_region_tracks(
    *,
    chrom: str,
    start: int,
    end: int,
    control_mml: str | Path | None,
    control_me: str | Path | None,
    target_mml: str | Path | None,
    target_me: str | Path | None,
    control_coverage: str | Path | None = None,
    target_coverage: str | Path | None = None,
    label_a: str = "Control",
    label_b: str = "Target",
    color_a: str = "#2980b9",
    color_b: str = "#e67e22",
    gtf_path: str | Path | None = None,
    promoter_upstream: int = 1000,
    pad_bp: int = 2000,
    smooth_win: int = 5,
    out_path: str | Path,
) -> Path:
    """Render the 4-panel region track figure.

    Panels (top to bottom):
      1. gene structure — every gene overlapping the window, drawn with
         exons + promoter + strand-arrows;
      2. smoothed ME (control vs target);
      3. smoothed MML (control vs target);
      4. smoothed coverage (control vs target), if coverage bedgraphs
         are supplied — otherwise the panel is omitted.

    The window is widened by ``pad_bp`` on each side so the gene body
    fits comfortably without flush-cutting promoters / exons.
    """
    safe_mkdir(os.path.dirname(str(out_path)))
    plot_lo = max(0, int(start) - int(pad_bp))
    plot_hi = int(end) + int(pad_bp)

    structures = load_gene_structures(
        str(gtf_path) if gtf_path else "", promoter_upstream=promoter_upstream,
    )
    region_genes = get_gene_structures_for_region(
        structures, chrom, (plot_lo, plot_hi),
    )

    # Read all bedgraphs, restricted to the padded window.
    ctrl_mml_df = _read_region(control_mml, chrom, plot_lo, plot_hi)
    case_mml_df = _read_region(target_mml, chrom, plot_lo, plot_hi)
    ctrl_me_df = _read_region(control_me, chrom, plot_lo, plot_hi)
    case_me_df = _read_region(target_me, chrom, plot_lo, plot_hi)
    has_cov = bool(control_coverage) or bool(target_coverage)
    if has_cov:
        ctrl_cov_df = _read_region(control_coverage, chrom, plot_lo, plot_hi)
        case_cov_df = _read_region(target_coverage, chrom, plot_lo, plot_hi)

    # Tall enough for a gene panel that may stack a few rows + an
    # optional coverage panel at the bottom.
    n_gene_rows = max(1, len(_stack_genes_by_row(region_genes)))
    n_signal = 3 if has_cov else 2
    height_ratios = [max(1.0, 0.8 * n_gene_rows)] + [2] * n_signal
    fig, axes = plt.subplots(
        1 + n_signal, 1,
        figsize=(14, 2.0 + 0.9 * n_gene_rows + 2.5 * n_signal),
        height_ratios=height_ratios,
        sharex=False,
        gridspec_kw={"hspace": 0.10},
    )
    ax_genes = axes[0]
    ax_me = axes[1]
    ax_mml = axes[2]
    ax_cov = axes[3] if has_cov else None

    _draw_gene_panel(ax_genes, region_genes, plot_lo, plot_hi)
    title_chr = chrom if chrom.startswith("chr") else f"chr{chrom}"
    ax_genes.set_title(
        f"{title_chr}:{plot_lo:,}-{plot_hi:,}",
        fontsize=12, fontweight="bold",
    )

    _plot_signal(
        ax_me,
        ctrl_df=ctrl_me_df, case_df=case_me_df,
        label_a=label_a, label_b=label_b,
        color_a=color_a, color_b=color_b,
        ylabel="Methylation Entropy (ME)",
        plot_lo=plot_lo, plot_hi=plot_hi, smooth_win=smooth_win,
    )
    _plot_signal(
        ax_mml,
        ctrl_df=ctrl_mml_df, case_df=case_mml_df,
        label_a=label_a, label_b=label_b,
        color_a=color_a, color_b=color_b,
        ylabel="Mean Methylation (MML)",
        plot_lo=plot_lo, plot_hi=plot_hi, smooth_win=smooth_win,
    )

    if ax_cov is not None:
        _plot_signal(
            ax_cov,
            ctrl_df=ctrl_cov_df, case_df=case_cov_df,
            label_a=label_a, label_b=label_b,
            color_a=color_a, color_b=color_b,
            ylabel="Coverage (reads)",
            plot_lo=plot_lo, plot_hi=plot_hi, smooth_win=smooth_win,
            ylim=None,  # coverage isn't 0–1 — autoscale
        )

    # All signal panels share the x-axis with the gene panel.
    bottom_signal_ax = ax_cov if ax_cov is not None else ax_mml
    for ax in (ax_me, ax_mml) + ((ax_cov,) if ax_cov is not None else ()):
        ax.set_xlim(plot_lo, plot_hi)
    bottom_signal_ax.set_xlabel(f"Genomic position ({title_chr})", fontsize=10)

    # `tight_layout` warns about the gene panel's axis("off"); the
    # warning is cosmetic and `bbox_inches="tight"` on savefig handles
    # the actual cropping. Suppress to keep test logs clean.
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout()
    out = Path(out_path)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


# ─── Backwards-compatible wrapper kept for callers that still use the v3
#     signature (gene_df + window_size). Re-routes to plot_region_tracks. ──

def plot_two_bedgraph_overlays(
    top_files: list[str],
    bottom_files: list[str],
    gene_df: pd.DataFrame | None,  # ignored (kept for ABI compat)
    chrom: str,
    start: int,
    end: int,
    out_path: str,
    window_size: int,
) -> None:
    """Compat shim: forwards to ``plot_region_tracks`` with default colours.

    ``top_files`` is expected to be ``[control_mml, target_mml]`` and
    ``bottom_files`` ``[control_me, target_me]``. ``gene_df`` is ignored
    — the new path re-parses the GTF on demand via ``load_gene_structures``.
    """
    if len(top_files) != 2 or len(bottom_files) != 2:
        raise ValueError(
            "plot_two_bedgraph_overlays expects exactly two control/target paths."
        )
    plot_region_tracks(
        chrom=chrom, start=start, end=end,
        control_mml=top_files[0], target_mml=top_files[1],
        control_me=bottom_files[0], target_me=bottom_files[1],
        gtf_path=None,
        smooth_win=int(max(1, window_size)),
        out_path=out_path,
    )
