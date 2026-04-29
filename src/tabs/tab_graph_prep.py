"""Tab 2 — Graph Preparation.

Three plot families on top of the user's chosen samples:
  * tracks            — control & target MML/ME along a genomic window
  * ME / MML scatter  — 2D-histogram of control_value vs target_value
  * arch landscape    — A | B | (B − A) MML × ME density with the
                        theoretical entropy arch
  * paired landscape  — paired-line scatter with direction arrows on
                        a configurable bin filter
"""

from __future__ import annotations

import logging
from pathlib import Path

from src.config import RESULTS_DIR
from src.help_text import SEC_GP_REGION, SEC_GP_SAMPLES
from src.io.gtf_utils import load_genes_from_gtf
from src.io.utils_io import ensure_writable_dir
from src.plots.scatter import (
    load_paired_bedgraphs,
    me_mml_scatter,
    paired_landscape,
    triple_landscape,
)
from src.plots.theme import apply_default_style
from src.plots.tracks import plot_two_bedgraph_overlays
from src.state import get_state, update_section
from src.ui.error_handler import show_error

logger = logging.getLogger(__name__)


# ─── Sample input panel ───────────────────────────────────────────────────

def _track_inputs(state, side: str) -> None:
    import streamlit as st

    cur = getattr(state.graph_prep, side)
    st.markdown(
        f'<div style="font-family:\'IBM Plex Mono\',monospace;'
        f'font-size:0.78rem;text-transform:uppercase;letter-spacing:0.06em;'
        f'color:#7a7a78;margin-bottom:0.4rem;">{side}</div>',
        unsafe_allow_html=True,
    )
    name = st.text_input("label", value=cur.name, key=f"{side}_name")
    color = st.color_picker("colour", value=cur.color, key=f"{side}_color")
    mml = st.text_input(
        "MML bedgraph", value=str(cur.mml_path or ""), key=f"{side}_mml",
        placeholder="/path/to/sample.mml.bedgraph",
    )
    me = st.text_input(
        "ME bedgraph", value=str(cur.me_path or ""), key=f"{side}_me",
        placeholder="/path/to/sample.me.bedgraph",
    )
    cov = st.text_input(
        "coverage bedgraph", value=str(cur.coverage_path or ""),
        key=f"{side}_cov",
        placeholder="/path/to/sample.coverage.bedgraph",
    )

    from src.state import TrackPlot
    new_track = TrackPlot(
        name=name or cur.name, color=color or cur.color,
        mml_path=Path(mml) if mml else None,
        me_path=Path(me) if me else None,
        coverage_path=Path(cov) if cov else None,
    )
    update_section("graph_prep", **{side: new_track})


# ─── Loaders ──────────────────────────────────────────────────────────────

def _paired_paths_ready(gp) -> bool:
    return all([
        gp.control.mml_path, gp.control.me_path,
        gp.target.mml_path, gp.target.me_path,
    ])


@show_error(user_message="Could not load bedgraphs.")
def _load_paired(gp):
    if not _paired_paths_ready(gp):
        return None
    return load_paired_bedgraphs(
        control_mml=gp.control.mml_path,
        control_me=gp.control.me_path,
        target_mml=gp.target.mml_path,
        target_me=gp.target.me_path,
    )


# ─── Renderers (decorated with show_error so traceback shows on screen) ──

@show_error(user_message="Track plot failed.")
def _render_tracks(state, out_path: Path, window: int) -> Path | None:
    gp = state.graph_prep
    if not _paired_paths_ready(gp):
        return None
    gtf_df = (
        load_genes_from_gtf(str(gp.gtf_path))
        if gp.gtf_path and Path(gp.gtf_path).exists()
        else None
    )
    plot_two_bedgraph_overlays(
        top_files=[str(gp.control.mml_path), str(gp.target.mml_path)],
        bottom_files=[str(gp.control.me_path), str(gp.target.me_path)],
        gene_df=gtf_df,
        chrom=gp.region_chrom, start=gp.region_start, end=gp.region_end,
        out_path=str(out_path), window_size=window,
    )
    return out_path


@show_error(user_message="ME / MML scatter failed.")
def _render_me_mml_scatter(gp, out_path: Path, log_scale: bool):
    df = _load_paired(gp)
    if df is None or df.empty:
        return None
    return me_mml_scatter(
        df,
        label_a=gp.control.name or "Control",
        label_b=gp.target.name or "Target",
        log_scale=log_scale,
        out_path=out_path,
    )


@show_error(user_message="Arch landscape plot failed.")
def _render_triple_landscape(gp, out_path: Path, log_scale: bool):
    df = _load_paired(gp)
    if df is None or df.empty:
        return None
    return triple_landscape(
        df,
        label_a=gp.control.name or "Control",
        label_b=gp.target.name or "Target",
        color_a=gp.control.color,
        color_b=gp.target.color,
        log_scale=log_scale,
        out_path=out_path,
    )


@show_error(user_message="Paired landscape plot failed.")
def _render_paired_landscape(
    gp, out_path: Path, *,
    fa_dim, fa_op, fa_val, fb_dim, fb_op, fb_val, max_lines,
):
    df = _load_paired(gp)
    if df is None or df.empty:
        return None
    return paired_landscape(
        df,
        label_a=gp.control.name or "Control",
        label_b=gp.target.name or "Target",
        color_a=gp.control.color,
        color_b=gp.target.color,
        filter_a_dim=fa_dim, filter_a_op=fa_op, filter_a_value=fa_val,
        filter_b_dim=fb_dim, filter_b_op=fb_op, filter_b_value=fb_val,
        max_lines=max_lines,
        out_path=out_path,
    )


# ─── Tab help text ────────────────────────────────────────────────────────

_TAB_HELP = """\
Three plot families compare your control and target bedgraphs:

* **Track plot** — control & target MML/ME along a genomic region, with
  optional gene panel from a GTF.
* **ME / MML scatter** — 2D-histogram of control vs target for each
  metric. Toggle linear / log colour scale.
* **Arch landscape** — three panels (control, target, target − control)
  showing density on the MML × ME plane with the theoretical
  binary-entropy arch overlaid.
* **Paired landscape** — pairs of bins drawn as faint lines from
  control → target, with direction arrows showing how many bins moved
  up vs down in entropy. Filter by `MML`, `|dMML|`, `ME`, or `|dME|`
  with a configurable threshold.

Bedgraph paths are auto-populated after a successful run on the **File
Preparation** tab. Sample labels and colours flow through to every plot.
"""


# ─── Main render() ────────────────────────────────────────────────────────

def render() -> None:
    import streamlit as st

    apply_default_style()
    state = get_state()
    gp = state.graph_prep

    out_dir = Path(ensure_writable_dir(str(RESULTS_DIR / "graph_prep"), "Plots out"))
    update_section("graph_prep", plots_dir=out_dir)

    paths_loaded = _paired_paths_ready(gp)
    with st.expander(
        "📖  How to use this tab" if not paths_loaded else "How to use this tab",
        expanded=not paths_loaded,
    ):
        st.markdown(_TAB_HELP)
        if not paths_loaded:
            st.info(
                "💡 Run the **File Preparation** tab first — paths land "
                "here automatically."
            )

    def _section_caption(text: str) -> None:
        st.markdown(
            f'<div class="caption" style="margin:-0.4rem 0 0.6rem;">{text}</div>',
            unsafe_allow_html=True,
        )

    # ── 01 · samples ──────────────────────────────────────────────
    st.markdown("### 01 · samples")
    _section_caption(SEC_GP_SAMPLES)
    c1, c2 = st.columns(2)
    with c1:
        _track_inputs(state, "control")
    with c2:
        _track_inputs(state, "target")

    # ── 02 · region (track plot) ─────────────────────────────────
    st.markdown("### 02 · region (track plot)")
    _section_caption(SEC_GP_REGION)
    cc, cs, ce = st.columns([2, 2, 2])
    with cc:
        chrom = st.text_input("chromosome", value=gp.region_chrom, placeholder="chr1")
    with cs:
        start = st.number_input(
            "start", value=int(gp.region_start), min_value=0, step=1000,
        )
    with ce:
        end = st.number_input(
            "end", value=int(gp.region_end), min_value=0, step=1000,
        )
    update_section(
        "graph_prep",
        region_chrom=chrom or "",
        region_start=int(start), region_end=int(end),
    )
    gtf = st.text_input("GTF (optional, gene panel)", value=str(gp.gtf_path or ""))
    update_section("graph_prep", gtf_path=Path(gtf) if gtf else None)

    window = st.slider("smoothing window (bins)", 1, 200, 25)
    if st.button(
        "RENDER TRACK PLOT", type="primary", use_container_width=True,
        key="btn_tracks",
    ):
        if not chrom or end <= start:
            st.error("specify a valid chrom / start / end")
        else:
            out = out_dir / f"tracks_{chrom}_{int(start)}_{int(end)}.png"
            res = _render_tracks(state, out, window=window)
            if res and Path(res).exists():
                st.image(str(res))

    st.markdown("---")

    # ── 03 · ME / MML scatter ─────────────────────────────────────
    st.markdown("### 03 · ME / MML scatter")
    _section_caption(
        "Two side-by-side 2D-histograms: control vs target for MML and ME. "
        "The dashed line is y = x (perfect agreement)."
    )
    log_scale_03 = st.toggle(
        "log10 colour scale",
        value=True,
        key="me_mml_log",
        help="Off = linear count colour scale; on = log10 count.",
    )
    if st.button(
        "RENDER ME / MML SCATTER", type="primary", use_container_width=True,
        key="btn_me_mml_scatter",
    ):
        if not paths_loaded:
            st.error("control and target MML + ME bedgraphs all required")
        else:
            out = out_dir / f"me_mml_{'log' if log_scale_03 else 'linear'}.png"
            fig = _render_me_mml_scatter(gp, out, log_scale=log_scale_03)
            if fig is not None:
                st.pyplot(fig)

    st.markdown("---")

    # ── 04 · arch landscape ───────────────────────────────────────
    st.markdown("### 04 · arch landscape")
    _section_caption(
        "Three panels (control, target, target − control) showing density "
        "on the MML × ME plane with the theoretical binary-entropy arch "
        "overlaid."
    )
    log_scale_04 = st.toggle(
        "log10 colour scale",
        value=True,
        key="arch_log",
    )
    if st.button(
        "RENDER ARCH LANDSCAPE", type="primary", use_container_width=True,
        key="btn_arch",
    ):
        if not paths_loaded:
            st.error("control and target MML + ME bedgraphs all required")
        else:
            out = out_dir / f"arch_{'log' if log_scale_04 else 'linear'}.png"
            fig = _render_triple_landscape(gp, out, log_scale=log_scale_04)
            if fig is not None:
                st.pyplot(fig)

    st.markdown("---")

    # ── 05 · paired landscape ─────────────────────────────────────
    st.markdown("### 05 · paired landscape")
    _section_caption(
        "Paired-line scatter from control → target on the MML × ME plane. "
        "Configure two AND-ed filters to keep only the bins of interest. "
        "The two arrows show how many filtered bins gained vs lost entropy."
    )

    dim_options = ["off", "MML", "|dMML|", "ME", "|dME|"]
    op_options = ["<", ">"]

    # Filter A — defaults to "entropy-shifted" (|dMML| < 0.1)
    st.markdown("**Filter A**")
    fa1, fa2, fa3 = st.columns([2, 1, 2])
    with fa1:
        fa_dim = st.selectbox(
            "dimension", options=dim_options, index=2,  # |dMML|
            key="pl_fa_dim",
        )
    with fa2:
        fa_op = st.selectbox("op", options=op_options, index=0, key="pl_fa_op")
    with fa3:
        fa_val = st.number_input(
            "threshold", min_value=0.0, max_value=1.0,
            value=0.10, step=0.01, format="%.2f", key="pl_fa_val",
        )

    # Filter B — defaults to "entropy-shifted" (|dME| > 0.4)
    st.markdown("**Filter B**")
    fb1, fb2, fb3 = st.columns([2, 1, 2])
    with fb1:
        fb_dim = st.selectbox(
            "dimension", options=dim_options, index=4,  # |dME|
            key="pl_fb_dim",
        )
    with fb2:
        fb_op = st.selectbox("op", options=op_options, index=1, key="pl_fb_op")
    with fb3:
        fb_val = st.number_input(
            "threshold", min_value=0.0, max_value=1.0,
            value=0.40, step=0.01, format="%.2f", key="pl_fb_val",
        )

    max_lines = st.slider(
        "max paired lines drawn (subsample for speed)",
        1_000, 100_000, 20_000, step=1_000, key="pl_max_lines",
    )

    if st.button(
        "RENDER PAIRED LANDSCAPE", type="primary", use_container_width=True,
        key="btn_paired",
    ):
        if not paths_loaded:
            st.error("control and target MML + ME bedgraphs all required")
        else:
            tag = "_".join([
                f"{fa_dim}{fa_op}{fa_val:g}".replace("|", ""),
                f"{fb_dim}{fb_op}{fb_val:g}".replace("|", ""),
            ])
            out = out_dir / f"paired_{tag}.png"
            fig = _render_paired_landscape(
                gp, out,
                fa_dim=None if fa_dim == "off" else fa_dim,
                fa_op=fa_op, fa_val=float(fa_val),
                fb_dim=None if fb_dim == "off" else fb_dim,
                fb_op=fb_op, fb_val=float(fb_val),
                max_lines=int(max_lines),
            )
            if fig is not None:
                st.pyplot(fig)
