"""Comparison plots for control vs target methylation.

Three plot families, all sharing the same paired-bin merge step:

  1. ``me_mml_scatter``    — side-by-side ME and MML 2D-histograms
                             (control_value on x, target_value on y).
                             Toggle log/linear colour scale.

  2. ``triple_landscape``  — A | B | (B − A) panels of MML × ME density
                             with the theoretical binary-entropy arch
                             overlaid. Toggle log/linear.

  3. ``paired_landscape``  — paired-line scatter with direction arrows
                             showing how many bins shifted up vs down
                             after a configurable filter
                             (``MML`` / ``|dMML|`` / ``ME`` / ``|dME|``,
                             ``<`` or ``>`` a threshold).

Helper:
  ``load_paired_bedgraphs(...)`` joins four bedgraphs on
  (chrom, start, end) into a single ``pandas.DataFrame`` with columns
  ``mml_a``, ``me_a``, ``mml_b``, ``me_b``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.io.bedgraph import read_bedgraph

FilterDim = Literal["MML", "|dMML|", "ME", "|dME|"]
FilterOp = Literal["<", ">"]


# ─── Theoretical entropy arch ─────────────────────────────────────────────

def theoretical_entropy_binary(p: np.ndarray) -> np.ndarray:
    """Shannon entropy of a Bernoulli(p), normalised to log2 base.

    H(p) = -p log2(p) - (1-p) log2(1-p), with H(0) = H(1) = 0.
    """
    p = np.asarray(p, dtype=np.float64)
    out = np.zeros_like(p)
    m = (p > 0) & (p < 1)
    pm = p[m]
    out[m] = -(pm * np.log2(pm) + (1 - pm) * np.log2(1 - pm))
    return out


# ─── Bedgraph join ────────────────────────────────────────────────────────

def load_paired_bedgraphs(
    *,
    control_mml: str | Path,
    control_me: str | Path,
    target_mml: str | Path,
    target_me: str | Path,
) -> pd.DataFrame:
    """Inner-join four bedgraphs on (chrom, start, end). Returns columns
    [chrom, start, end, mml_a, me_a, mml_b, me_b] where _a is control
    and _b is target."""
    cols = ["chrom", "start", "end"]

    def _r(p, name):
        df = read_bedgraph(p)
        if df.empty:
            return pd.DataFrame(columns=cols + [name])
        return df.rename(columns={"value": name})[cols + [name]]

    df = (
        _r(control_mml, "mml_a")
        .merge(_r(control_me, "me_a"), on=cols, how="inner")
        .merge(_r(target_mml, "mml_b"), on=cols, how="inner")
        .merge(_r(target_me, "me_b"), on=cols, how="inner")
    )
    return df


# ─── Colour map helpers (lifted from the user's reference script) ─────────

def _seq_cmap(hex_color: str, lo: float = 0.2) -> mcolors.Colormap:
    rgb = mcolors.to_rgb(hex_color)
    r0 = 1.0 - lo * (1.0 - rgb[0])
    g0 = 1.0 - lo * (1.0 - rgb[1])
    b0 = 1.0 - lo * (1.0 - rgb[2])
    return mcolors.LinearSegmentedColormap.from_list(
        "seq", [(1, 1, 1), (r0, g0, b0), rgb], N=256,
    )


def _div_cmap(c_hex: str, t_hex: str) -> mcolors.Colormap:
    return mcolors.LinearSegmentedColormap.from_list(
        "div", [mcolors.to_rgb(c_hex), (1, 1, 1), mcolors.to_rgb(t_hex)], N=512,
    )


# ─── Plot 1: ME and MML control-vs-target scatters ───────────────────────

def me_mml_scatter(
    df: pd.DataFrame,
    *,
    label_a: str,
    label_b: str,
    log_scale: bool = True,
    gridsize: int = 80,
    axis_cap: float = 1.0,
    vmax_pct: float = 70.0,
    out_path: Path | str | None = None,
) -> plt.Figure:
    """Side-by-side 2D-histogram scatters: control vs target for MML
    and for ME. ``df`` must come from ``load_paired_bedgraphs``."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, metric, title in (
        (axes[0], "mml", "MML"),
        (axes[1], "me", "ME"),
    ):
        x = df[f"{metric}_a"].to_numpy()
        y = df[f"{metric}_b"].to_numpy()
        if x.size == 0:
            ax.text(0.5, 0.5, "no overlapping bins",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
            continue

        rng = [[0.0, axis_cap], [0.0, axis_cap]]
        H, _, _ = np.histogram2d(x, y, bins=int(gridsize), range=rng)
        H = H.T
        Hm = H.astype(float)
        Hm[Hm == 0] = np.nan

        if log_scale:
            Hd = np.log10(Hm)
            cbar_lbl = "log10(count)"
            finite = Hd[np.isfinite(Hd)]
            vmin = float(np.nanmin(finite)) if finite.size else 0.0
            vmax = float(np.percentile(finite, vmax_pct)) if finite.size else 1.0
        else:
            Hd = Hm
            cbar_lbl = "count"
            vmin = 1.0
            finite = Hd[np.isfinite(Hd)]
            vmax = float(np.percentile(finite, 50)) if finite.size else 1.0

        im = ax.imshow(
            Hd, origin="lower", extent=[0, axis_cap, 0, axis_cap],
            aspect="auto", cmap="viridis",
            vmin=vmin, vmax=max(vmax, vmin + 1e-6), interpolation="nearest",
        )
        # y = x reference line
        ax.plot([0, axis_cap], [0, axis_cap], "--", color="black", lw=1, alpha=0.6)

        ax.set_xlim(0, axis_cap)
        ax.set_ylim(0, axis_cap)
        ax.set_xlabel(f"{label_a} {title}")
        ax.set_ylabel(f"{label_b} {title}")
        scale = "log" if log_scale else "linear"
        ax.set_title(f"{title} — {label_a} vs {label_b} ({scale})", fontsize=12)
        ax.grid(alpha=0.15)
        fig.colorbar(im, ax=ax, label=cbar_lbl, fraction=0.046, pad=0.02)

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
    return fig


# ─── Plot 2: triple landscape (A | B | diff) with entropy arch ───────────

def triple_landscape(
    df: pd.DataFrame,
    *,
    label_a: str,
    label_b: str,
    color_a: str = "#1f77b4",
    color_b: str = "#ff7f0e",
    log_scale: bool = True,
    gridsize: int = 80,
    axis_cap: float = 1.0,
    vmax_pct: float = 70.0,
    diff_linthresh: float = 10.0,
    linear_vmax: float | None = None,
    out_path: Path | str | None = None,
) -> plt.Figure:
    """3-panel landscape: A | B | (B − A), each MML × ME with the
    theoretical binary-entropy arch overlaid.

    Adapted from the user's reference ``render_triple_scatter`` so the
    colour-mapping behaviour is preserved.
    """
    cap = float(axis_cap)
    nbins = int(gridsize)
    rng = [[0.0, cap], [0.0, cap]]

    mml_a = df["mml_a"].to_numpy()
    me_a = df["me_a"].to_numpy()
    mml_b = df["mml_b"].to_numpy()
    me_b = df["me_b"].to_numpy()

    H_a, _, _ = np.histogram2d(mml_a, me_a, bins=nbins, range=rng)
    H_b, _, _ = np.histogram2d(mml_b, me_b, bins=nbins, range=rng)
    H_a = H_a.T
    H_b = H_b.T

    H_a_m = H_a.astype(float)
    H_a_m[H_a_m == 0] = np.nan
    H_b_m = H_b.astype(float)
    H_b_m[H_b_m == 0] = np.nan
    H_diff_m = (H_b - H_a).astype(float)
    H_diff_m[H_diff_m == 0] = np.nan

    if log_scale:
        H_a_d = np.log10(H_a_m)
        H_b_d = np.log10(H_b_m)
        cbar_dens = "log10(count)"
        all_v = np.concatenate([
            H_a_d[np.isfinite(H_a_d)], H_b_d[np.isfinite(H_b_d)],
        ])
        s_vmin = float(np.nanmin(all_v)) if all_v.size else 0.0
        s_vmax = float(np.percentile(all_v, vmax_pct)) if all_v.size else 1.0
        H_diff_d = np.sign(H_diff_m) * np.log10(
            1 + np.abs(H_diff_m) / diff_linthresh
        )
        cbar_diff = f"symlog10(Δ, t={diff_linthresh:.0f})"
        d_vmax = np.log10(1 + 10 ** s_vmax / diff_linthresh)
    else:
        H_a_d = H_a_m
        H_b_d = H_b_m
        cbar_dens = "count"
        s_vmin = 1.0
        if linear_vmax is not None:
            s_vmax = float(linear_vmax)
        else:
            all_v = np.concatenate([
                H_a_d[np.isfinite(H_a_d)], H_b_d[np.isfinite(H_b_d)],
            ])
            s_vmax = float(np.percentile(all_v, 50)) if all_v.size else 1.0
        H_diff_d = H_diff_m
        cbar_diff = "Δ count"
        d_vmax = s_vmax

    d_vmax = max(d_vmax, 1e-6)

    cm_a = _seq_cmap(color_a)
    cm_b = _seq_cmap(color_b)
    diff_cm = _div_cmap(color_a, color_b)

    p = np.linspace(0.001, min(cap, 0.999), 400)
    arch = np.clip(theoretical_entropy_binary(p), 0, cap)

    extent = [0, cap, 0, cap]
    ticks = [t for t in (0, 0.25, 0.5, 0.75, 1.0) if t <= cap]
    scale_label = "log" if log_scale else "linear"

    fig, axes = plt.subplots(1, 3, figsize=(22, 7))

    im0 = axes[0].imshow(
        H_a_d, origin="lower", extent=extent, aspect="auto",
        cmap=cm_a, vmin=s_vmin, vmax=s_vmax, interpolation="nearest",
    )
    axes[0].plot(p, arch, "k", linewidth=2)
    axes[0].set_title(f"{label_a} — ME vs MML ({scale_label})", fontsize=13)
    fig.colorbar(im0, ax=axes[0], label=cbar_dens, fraction=0.046, pad=0.02)

    im1 = axes[1].imshow(
        H_b_d, origin="lower", extent=extent, aspect="auto",
        cmap=cm_b, vmin=s_vmin, vmax=s_vmax, interpolation="nearest",
    )
    axes[1].plot(p, arch, "k", linewidth=2)
    axes[1].set_title(f"{label_b} — ME vs MML ({scale_label})", fontsize=13)
    fig.colorbar(im1, ax=axes[1], label=cbar_dens, fraction=0.046, pad=0.02)

    im2 = axes[2].imshow(
        H_diff_d, origin="lower", extent=extent, aspect="auto",
        cmap=diff_cm, vmin=-d_vmax, vmax=d_vmax, interpolation="nearest",
    )
    axes[2].plot(p, arch, "k", linewidth=2)
    axes[2].set_title(
        f"Difference — {label_b} − {label_a} ({scale_label})", fontsize=13,
    )
    fig.colorbar(im2, ax=axes[2], label=cbar_diff, fraction=0.046, pad=0.02)

    for ax in axes:
        ax.set_xlim(0, cap)
        ax.set_ylim(0, cap)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel("Mean Methylation Level (MML)", fontsize=11)
        ax.set_ylabel("Methylation Entropy (ME)", fontsize=11)
        ax.grid(alpha=0.15)

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
    return fig


# ─── Plot 3: paired-bin landscape ────────────────────────────────────────

def _apply_filter(
    df: pd.DataFrame, dim: FilterDim, op: FilterOp, value: float,
) -> pd.DataFrame:
    """Apply a single ``dim {op} value`` filter to a paired bedgraph df.

    Allowed ``dim``:
      - ``"MML"``    — keep bins where mml_a is on the chosen side
                       (uses the control sample's value as anchor).
      - ``"ME"``     — same, on me_a.
      - ``"|dMML|"`` — keep bins where |mml_b - mml_a| {op} value.
      - ``"|dME|"``  — keep bins where |me_b - me_a| {op} value.
    """
    if dim == "MML":
        col = df["mml_a"]
    elif dim == "ME":
        col = df["me_a"]
    elif dim == "|dMML|":
        col = (df["mml_b"] - df["mml_a"]).abs()
    elif dim == "|dME|":
        col = (df["me_b"] - df["me_a"]).abs()
    else:
        raise ValueError(f"unknown filter dim {dim!r}")

    mask = (col < value) if op == "<" else (col > value)
    return df[mask]


def paired_landscape(
    df: pd.DataFrame,
    *,
    label_a: str,
    label_b: str,
    color_a: str = "#1f77b4",
    color_b: str = "#ff7f0e",
    filter_a_dim: FilterDim | None = "|dMML|",
    filter_a_op: FilterOp = "<",
    filter_a_value: float = 0.1,
    filter_b_dim: FilterDim | None = "|dME|",
    filter_b_op: FilterOp = ">",
    filter_b_value: float = 0.4,
    max_lines: int = 20_000,
    out_path: Path | str | None = None,
) -> plt.Figure:
    """Paired-line scatter (filtered).

    Each kept bin is drawn as a faint line from
    ``(mml_a, me_a) → (mml_b, me_b)``, plus two small dot scatters
    coloured by sample. The theoretical entropy arch is drawn on top.

    Bin counts (Δ ME ↑ / Δ ME ↓) appear in a caption below the axes
    rather than as in-plot arrows.

    Filters are AND-ed; pass ``filter_*_dim=None`` to disable.
    """
    sub = df.copy()
    if filter_a_dim is not None:
        sub = _apply_filter(sub, filter_a_dim, filter_a_op, filter_a_value)
    if filter_b_dim is not None:
        sub = _apply_filter(sub, filter_b_dim, filter_b_op, filter_b_value)

    n_bins = len(sub)
    # Slightly taller than 1:1 so the caption row below has breathing room.
    fig, ax = plt.subplots(figsize=(9, 7.6))

    if n_bins == 0:
        ax.text(
            0.5, 0.5, "no bins pass filter",
            ha="center", va="center", transform=ax.transAxes,
            fontsize=12, color="#a04040",
        )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Methylation level (MML)")
        ax.set_ylabel("Entropy (ME)")
        if out_path:
            fig.savefig(out_path, dpi=200, bbox_inches="tight")
        return fig

    plot_sub = (
        sub.sample(min(n_bins, max_lines), random_state=1)
        if n_bins > max_lines
        else sub
    )

    # Paired faint lines
    for _, row in plot_sub.iterrows():
        ax.plot(
            [row["mml_a"], row["mml_b"]],
            [row["me_a"], row["me_b"]],
            lw=0.5, alpha=0.03, color="black", zorder=1,
        )

    # Per-sample scatter
    ax.scatter(
        plot_sub["mml_a"], plot_sub["me_a"],
        s=6, alpha=0.3, color=color_a, label=label_a, zorder=3,
    )
    ax.scatter(
        plot_sub["mml_b"], plot_sub["me_b"],
        s=6, alpha=0.3, color=color_b, label=label_b, zorder=3,
    )

    # Theoretical entropy arch
    p = np.linspace(0, 1, 400)
    ax.plot(
        p, theoretical_entropy_binary(p),
        lw=2, color="black", alpha=0.9, zorder=5,
        label="Theoretical arch",
    )

    # Bin-count caption (outside the plot, so the points stay clean).
    dme = (sub["me_b"] - sub["me_a"]).to_numpy()
    n_up = int((dme > 0).sum())
    n_down = int((dme < 0).sum())
    n_zero = int((dme == 0).sum())

    # Title summarising filters
    title_parts = []
    if filter_a_dim is not None:
        title_parts.append(f"{filter_a_dim} {filter_a_op} {filter_a_value:g}")
    if filter_b_dim is not None:
        title_parts.append(f"{filter_b_dim} {filter_b_op} {filter_b_value:g}")
    filt_str = " AND ".join(title_parts) if title_parts else "all bins"
    ax.set_title(
        f"{label_a} → {label_b}\n{filt_str} | n={n_bins:,} bins",
        fontweight="bold",
    )
    ax.set_xlabel("Methylation level (MML)")
    ax.set_ylabel("Entropy (ME)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.05)
    ax.legend(loc="upper left")

    fig.tight_layout()
    # Make room at the bottom for the caption block.
    fig.subplots_adjust(bottom=0.18)
    caption = (
        f"Δ ME ↑ ({label_b} > {label_a}):  n = {n_up:,}        "
        f"Δ ME ↓ ({label_b} < {label_a}):  n = {n_down:,}"
    )
    if n_zero:
        caption += f"        Δ ME = 0:  n = {n_zero:,}"
    fig.text(
        0.5, 0.04, caption,
        ha="center", va="bottom",
        fontsize=11, fontweight="bold", color="#222222",
    )

    if out_path:
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
    return fig
