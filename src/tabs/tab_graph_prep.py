"""Tab 2 — Graph Preparation: track plots, scatter, distributions."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from src.config import RESULTS_DIR
from src.help_text import (
    SEC_GP_DISTRIBUTION,
    SEC_GP_REGION,
    SEC_GP_SAMPLES,
    SEC_GP_SCATTER,
)
from src.io.bedgraph import read_bedgraph
from src.io.gtf_utils import load_genes_from_gtf
from src.io.utils_io import ensure_writable_dir
from src.plots.distributions import metric_distribution
from src.plots.scatter import hexbin_density, paired_scatter
from src.plots.theme import apply_default_style
from src.plots.tracks import plot_two_bedgraph_overlays
from src.state import get_state, update_section
from src.ui.error_handler import show_error

logger = logging.getLogger(__name__)


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
    mml = st.text_input("MML bedgraph", value=str(cur.mml_path or ""), key=f"{side}_mml",
                        placeholder="/path/to/sample.mml.bedgraph")
    me = st.text_input("ME bedgraph", value=str(cur.me_path or ""), key=f"{side}_me",
                       placeholder="/path/to/sample.me.bedgraph")
    cov = st.text_input("coverage bedgraph", value=str(cur.coverage_path or ""),
                        key=f"{side}_cov",
                        placeholder="/path/to/sample.coverage.bedgraph")

    from src.state import TrackPlot
    new_track = TrackPlot(
        name=name or cur.name, color=color or cur.color,
        mml_path=Path(mml) if mml else None,
        me_path=Path(me) if me else None,
        coverage_path=Path(cov) if cov else None,
    )
    update_section("graph_prep", **{side: new_track})


@show_error(user_message="Track plot failed.")
def _render_tracks(state, out_path: Path, window: int) -> Path | None:
    gp = state.graph_prep
    if not (gp.control.mml_path and gp.target.mml_path
            and gp.control.me_path and gp.target.me_path):
        return None
    gtf_df = (load_genes_from_gtf(str(gp.gtf_path))
              if gp.gtf_path and Path(gp.gtf_path).exists() else None)
    plot_two_bedgraph_overlays(
        top_files=[str(gp.control.mml_path), str(gp.target.mml_path)],
        bottom_files=[str(gp.control.me_path), str(gp.target.me_path)],
        gene_df=gtf_df,
        chrom=gp.region_chrom, start=gp.region_start, end=gp.region_end,
        out_path=str(out_path), window_size=window,
    )
    return out_path


@show_error(user_message="Scatter plot failed.")
def _render_scatter(
    df_x_path: Path, df_y_path: Path, *, kind: str,
    x_label: str, y_label: str, out_path: Path, subsample: int,
):
    df_x = read_bedgraph(df_x_path)
    df_y = read_bedgraph(df_y_path)
    if df_x.empty or df_y.empty:
        return None
    if kind == "hexbin":
        return hexbin_density(df_x, df_y, x_label=x_label, y_label=y_label, out_path=out_path)
    return paired_scatter(
        df_x, df_y, x_label=x_label, y_label=y_label,
        subsample=subsample, out_path=out_path,
    )


_TAB_HELP = """\
Renders three plot families from your bedgraphs:

* **track plot** — control & target MML/ME along a genomic region, with
  optional gene panel from a GTF.
* **scatter / hexbin** — control vs target for one metric, density-coloured.
* **distribution** — violin / box / histogram of values for one bedgraph.

Bedgraph paths are auto-populated after a successful run on the **File
Preparation** tab. You can also paste paths in by hand.
"""


def render() -> None:
    import streamlit as st

    apply_default_style()
    state = get_state()
    gp = state.graph_prep

    out_dir = Path(ensure_writable_dir(str(RESULTS_DIR / "graph_prep"), "Plots out"))
    update_section("graph_prep", plots_dir=out_dir)

    # ── Smart help: only auto-expanded if no bedgraphs are loaded yet ──
    paths_loaded = bool(
        gp.target.me_path or gp.target.mml_path
        or gp.control.me_path or gp.control.mml_path
    )
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

    # ── Samples ──
    st.markdown("### 01 · samples")
    _section_caption(SEC_GP_SAMPLES)
    c1, c2 = st.columns(2)
    with c1:
        _track_inputs(state, "control")
    with c2:
        _track_inputs(state, "target")

    # ── Region ──
    st.markdown("### 02 · region")
    _section_caption(SEC_GP_REGION)
    cc, cs, ce = st.columns([2, 2, 2])
    with cc:
        chrom = st.text_input("chromosome", value=gp.region_chrom,
                              placeholder="chr1")
    with cs:
        start = st.number_input("start", value=int(gp.region_start),
                                min_value=0, step=1000)
    with ce:
        end = st.number_input("end", value=int(gp.region_end),
                              min_value=0, step=1000)
    update_section(
        "graph_prep",
        region_chrom=chrom or "",
        region_start=int(start), region_end=int(end),
    )
    gtf = st.text_input("GTF (optional, gene panel)", value=str(gp.gtf_path or ""))
    update_section("graph_prep", gtf_path=Path(gtf) if gtf else None)

    window = st.slider("smoothing window (bins)", 1, 200, 25)
    if st.button("RENDER TRACK PLOT", type="primary", use_container_width=True,
                 key="btn_tracks"):
        if not chrom or end <= start:
            st.error("specify a valid chrom / start / end")
        else:
            out = out_dir / f"tracks_{chrom}_{int(start)}_{int(end)}.png"
            res = _render_tracks(state, out, window=window)
            if res and Path(res).exists():
                st.image(str(res))

    st.markdown("---")
    # ── Scatter ──
    st.markdown("### 03 · scatter (control vs target)")
    _section_caption(SEC_GP_SCATTER)
    cs1, cs2, cs3 = st.columns([1, 1, 2])
    with cs1:
        metric = st.selectbox("metric", options=["MML", "ME"])
    with cs2:
        kind = st.selectbox("kind", options=["hexbin", "scatter"])
    with cs3:
        subsample = st.slider("subsample (scatter only)", 1000, 200_000,
                              gp.scatter_subsample, step=1000)
    update_section("graph_prep", scatter_subsample=int(subsample))
    if st.button("RENDER SCATTER", type="primary", use_container_width=True,
                 key="btn_scatter"):
        path_x = gp.control.mml_path if metric == "MML" else gp.control.me_path
        path_y = gp.target.mml_path if metric == "MML" else gp.target.me_path
        if not (path_x and path_y):
            st.error("both control and target paths required")
        else:
            out = out_dir / f"scatter_{metric}_{kind}.png"
            fig = _render_scatter(
                path_x, path_y, kind=kind,
                x_label=f"control {metric}", y_label=f"target {metric}",
                out_path=out, subsample=int(subsample),
            )
            if fig is not None:
                st.pyplot(fig)

    st.markdown("---")
    # ── Distributions ──
    st.markdown("### 04 · distribution")
    _section_caption(SEC_GP_DISTRIBUTION)
    cd1, cd2 = st.columns([3, 1])
    with cd1:
        dist_path = st.text_input(
            "bedgraph",
            value=str(gp.target.me_path or ""),
            placeholder="/path/to/sample.me.bedgraph",
        )
    with cd2:
        dist_kind = st.selectbox("plot", options=["violin", "box", "hist"])
    if st.button("RENDER DISTRIBUTION", type="primary", use_container_width=True,
                 key="btn_dist"):
        df = read_bedgraph(dist_path) if dist_path else pd.DataFrame()
        if df.empty:
            st.error("bedgraph empty or missing")
        else:
            df["sample"] = "data"
            fig = metric_distribution(df, "value", by="sample", kind=dist_kind)
            st.pyplot(fig)
