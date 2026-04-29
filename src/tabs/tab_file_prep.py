"""Tab 1 — File Preparation.

Streamlit UI that wires user inputs (single sample or control+target,
each one a BAM / TSV / folder-of-BAMs) to the pipeline orchestrator.

After a successful run, the produced bedgraph paths are written into
``state.graph_prep.{control,target}`` so the user can switch to Tab 2
and start plotting immediately.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path

from src.config import GENOMES, RESULTS_DIR, assets_for, ensure_dirs
from src.constants import (
    ENTROPY_MODE_HELP,
    ENTROPY_MODE_LABELS,
    ENTROPY_MODE_TERNARY,
    ENTROPY_MODES,
)
from src.help_text import (
    GUI_GETTING_STARTED,
    GUI_GLOSSARY,
    GUI_HELP_INPUT_KIND,
    GUI_HELP_K,
    GUI_HELP_METHYL_THRESHOLD,
    GUI_HELP_MIN_COV,
    GUI_HELP_PAIR_MODE,
    SEC_FP_BIN_PARAMS,
    SEC_FP_ENTROPY_MODE,
    SEC_FP_OUTPUT,
    SEC_FP_REFERENCE,
    SEC_FP_SAMPLES,
)
from src.io.utils_io import disk_free_gb, ensure_writable_dir
from src.pipelines.orchestrator import run_pipeline
from src.pipelines.ternary_entropy import required_coverage_for_k
from src.state import SampleSpec, TrackPlot, get_state, update_section
from src.ui.error_handler import show_error
from src.ui.progress import StreamlitProgress

logger = logging.getLogger(__name__)

INPUT_KIND_LABELS = {
    "bam": "single BAM",
    "tsv": "pre-computed TSV",
    "bam_folder": "folder of BAMs (auto-merge)",
}


def _sample_inputs(
    spec: SampleSpec, *, key_prefix: str, default_label: str,
) -> SampleSpec:
    """Render inputs for one sample. Returns an updated SampleSpec."""
    import streamlit as st

    label = st.text_input(
        "label", value=spec.label or default_label,
        key=f"{key_prefix}_label",
    )

    kind = st.radio(
        "input kind",
        options=["bam", "tsv", "bam_folder"],
        index=["bam", "tsv", "bam_folder"].index(spec.input_kind),
        format_func=lambda k: INPUT_KIND_LABELS[k],
        horizontal=True,
        key=f"{key_prefix}_kind",
        help=GUI_HELP_INPUT_KIND,
    )

    new = SampleSpec(label=label or default_label, input_kind=kind)
    if kind == "bam":
        bam = st.text_input(
            "BAM path",
            value=str(spec.bam_path or ""),
            placeholder="/path/to/sample.bam",
            key=f"{key_prefix}_bam",
        )
        new.bam_path = Path(bam) if bam else None
    elif kind == "tsv":
        tsv = st.text_input(
            "modkit TSV path",
            value=str(spec.tsv_path or ""),
            placeholder="/path/to/modkit_extract.tsv",
            key=f"{key_prefix}_tsv",
        )
        new.tsv_path = Path(tsv) if tsv else None
    else:  # bam_folder
        folder = st.text_input(
            "folder of BAMs",
            value=str(spec.bam_folder or ""),
            placeholder="/path/to/bams_dir",
            key=f"{key_prefix}_folder",
            help="All `*.bam` files found recursively will be merged, "
                 "sorted and indexed before modkit extract runs.",
        )
        new.bam_folder = Path(folder) if folder else None
        if folder and Path(folder).is_dir():
            try:
                from src.pipelines.bam_utils import find_bams
                found = find_bams(folder)
                st.markdown(
                    f'<div class="caption">found {len(found)} BAMs · '
                    f'will be merged → sorted → indexed</div>',
                    unsafe_allow_html=True,
                )
            except OSError:
                pass
    return new


@show_error(user_message="File preparation failed. See traceback below.")
def _run(out_dir: Path) -> dict:

    state = get_state()
    fp = state.file_prep
    fasta = (fp.custom_fasta if fp.custom_fasta
             else assets_for(fp.genome).fasta)
    out_dir = Path(ensure_writable_dir(str(out_dir), "Output dir"))
    progress_lines: list[str] = []

    def progress(msg: str) -> None:
        progress_lines.append(msg)
        logger.info("pipeline: %s", msg)

    target_result = None
    control_result = None
    with StreamlitProgress("shannonpore pipeline") as bar:
        def progress_and_bar(msg: str) -> None:
            progress(msg)
            bar.status(msg)

        target_result, control_result = run_pipeline(
            target_spec=fp.target,
            control_spec=fp.control if fp.pair_mode else None,
            out_dir=out_dir,
            fasta=fasta,
            entropy_mode=fp.entropy_mode,
            cpgs_per_bin=fp.cpgs_per_bin,
            min_coverage=fp.min_coverage,
            methyl_threshold=fp.methyl_threshold,
            threads=fp.modkit_threads,
            progress_cb=progress_and_bar,
            pct_cb=bar.update,
            status_cb=bar.status_line,  # modkit's in-place progress bar
        )

    summary = {
        "fasta": str(fasta),
        "genome": fp.genome,
        "entropy_mode": fp.entropy_mode,
        "pair_mode": fp.pair_mode,
        "target": _summarize(target_result),
        "control": _summarize(control_result) if control_result else None,
        "progress": progress_lines,
    }
    Path(out_dir, "run_summary.json").write_text(
        json.dumps(summary, default=str, indent=2)
    )

    update_section(
        "file_prep",
        last_run_id=datetime.now().strftime("%Y%m%d-%H%M%S"),
        last_run_summary=summary,
    )
    _autoload_graph_prep(target_result, control_result, fp)
    return summary


def _summarize(result) -> dict | None:
    if result is None:
        return None
    return {
        "label": result.label,
        "out_prefix": str(result.out_prefix),
        "bam_used": str(result.bam_used) if result.bam_used else None,
        "tsv_path": str(result.tsv_path) if result.tsv_path else None,
        "me_bedgraph": str(result.me_bedgraph),
        "mml_bedgraph": str(result.mml_bedgraph),
        "coverage_bedgraph": str(result.coverage_bedgraph),
        "mhml_bedgraph": (str(result.mhml_bedgraph)
                          if result.mhml_bedgraph else None),
    }


def _autoload_graph_prep(target_result, control_result, fp) -> None:
    """Wire produced bedgraph paths into state.graph_prep AND into the
    Streamlit-widget-keyed session_state so the Plotting tab is ready
    to go without copy-paste.

    Two writes are needed because Streamlit widgets with a ``key=``
    own their session_state slot: on first render they snapshot the
    ``value=`` argument, on subsequent renders they read from
    ``st.session_state[key]`` and *ignore* the new ``value=`` argument
    entirely. So updating the AppState dataclass alone won't move the
    UI — we have to force-overwrite the widget keys too. This works
    because Streamlit allows writing to a key before the widget
    instantiates (we run during the File Prep tab's button handler,
    which fires before the Graph Prep tab re-renders).
    """
    import streamlit as st

    target_track = TrackPlot(
        name=target_result.label,
        color="#ff7f0e",
        mml_path=target_result.mml_bedgraph,
        me_path=target_result.me_bedgraph,
        coverage_path=target_result.coverage_bedgraph,
    )
    update_section("graph_prep", target=target_track)
    _push_track_to_widgets(st, "target", target_track)

    if control_result is not None:
        control_track = TrackPlot(
            name=control_result.label,
            color="#1f77b4",
            mml_path=control_result.mml_bedgraph,
            me_path=control_result.me_bedgraph,
            coverage_path=control_result.coverage_bedgraph,
        )
        update_section("graph_prep", control=control_track)
        _push_track_to_widgets(st, "control", control_track)
    else:
        # Single-sample mode — make `control` a copy of target so the
        # Plotting tab still shows something sensible (user can swap).
        copied = TrackPlot(
            name=target_result.label,
            color="#1f77b4",
            mml_path=target_result.mml_bedgraph,
            me_path=target_result.me_bedgraph,
            coverage_path=target_result.coverage_bedgraph,
        )
        update_section("graph_prep", control=copied)
        _push_track_to_widgets(st, "control", copied)


def _push_track_to_widgets(st, side: str, track: TrackPlot) -> None:
    """Force-overwrite the Streamlit-widget-keyed slots that
    tab_graph_prep._track_inputs reads from.

    Keys must match those used in tab_graph_prep.py:
        f"{side}_name", f"{side}_color", f"{side}_mml",
        f"{side}_me",   f"{side}_cov"
    """
    st.session_state[f"{side}_name"] = track.name
    st.session_state[f"{side}_color"] = track.color
    st.session_state[f"{side}_mml"] = str(track.mml_path or "")
    st.session_state[f"{side}_me"] = str(track.me_path or "")
    st.session_state[f"{side}_cov"] = str(track.coverage_path or "")


def render() -> None:
    import streamlit as st

    ensure_dirs()
    state = get_state()
    fp = state.file_prep

    # ── Getting started — only auto-expanded for first-time users ──
    first_time = not fp.last_run_id
    with st.expander(
        "📖  How to use this tab" if first_time else "How to use this tab",
        expanded=first_time,
    ):
        st.markdown(GUI_GETTING_STARTED)
        st.markdown("---")
        st.markdown("**Glossary**")
        st.markdown(GUI_GLOSSARY)

    def _section_caption(text: str) -> None:
        st.markdown(
            f'<div class="caption" style="margin:-0.4rem 0 0.6rem;">{text}</div>',
            unsafe_allow_html=True,
        )

    # ── Sample design ──
    st.markdown("### 01 · sample design")
    _section_caption(SEC_FP_SAMPLES)
    pair_mode = st.toggle(
        "compare control vs target (pair mode)",
        value=fp.pair_mode,
        help=GUI_HELP_PAIR_MODE,
    )
    update_section("file_prep", pair_mode=pair_mode)

    if pair_mode:
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(
                '<div style="font-family:\'IBM Plex Mono\',monospace;'
                'font-size:0.78rem;text-transform:uppercase;letter-spacing:0.06em;'
                'color:#0f4c75;margin-bottom:0.3rem;">control</div>',
                unsafe_allow_html=True,
            )
            new_control = _sample_inputs(
                fp.control, key_prefix="ctrl", default_label="control",
            )
            update_section("file_prep", control=new_control)
        with c2:
            st.markdown(
                '<div style="font-family:\'IBM Plex Mono\',monospace;'
                'font-size:0.78rem;text-transform:uppercase;letter-spacing:0.06em;'
                'color:#b86b3a;margin-bottom:0.3rem;">target</div>',
                unsafe_allow_html=True,
            )
            new_target = _sample_inputs(
                fp.target, key_prefix="tgt", default_label="target",
            )
            update_section("file_prep", target=new_target)
    else:
        st.markdown(
            '<div style="font-family:\'IBM Plex Mono\',monospace;'
            'font-size:0.78rem;text-transform:uppercase;letter-spacing:0.06em;'
            'color:#7a7a78;margin-bottom:0.3rem;">single sample</div>',
            unsafe_allow_html=True,
        )
        new_target = _sample_inputs(
            fp.target, key_prefix="single", default_label="sample",
        )
        update_section("file_prep", target=new_target)

    col_disk, _ = st.columns([1, 3])
    with col_disk:
        st.metric("free disk", f"{disk_free_gb(str(RESULTS_DIR)):.0f} GiB")

    # ── Reference ──
    st.markdown("### 02 · reference")
    _section_caption(SEC_FP_REFERENCE)
    col_g, col_f = st.columns([1, 3])
    with col_g:
        genome = st.selectbox(
            "genome", options=list(GENOMES),
            index=list(GENOMES).index(fp.genome),
        )
    with col_f:
        custom = st.text_input(
            "custom FASTA (optional, overrides bundle)",
            value=str(fp.custom_fasta or ""),
            placeholder=str(assets_for(genome).fasta),
        )
    update_section(
        "file_prep",
        genome=genome,
        custom_fasta=Path(custom) if custom else None,
    )

    # ── Entropy mode ──
    st.markdown("### 03 · entropy mode")
    _section_caption(SEC_FP_ENTROPY_MODE)
    mode_idx = (
        ENTROPY_MODES.index(fp.entropy_mode)
        if fp.entropy_mode in ENTROPY_MODES else 0
    )
    mode = st.radio(
        "5hmC handling",
        options=ENTROPY_MODES,
        index=mode_idx,
        format_func=lambda m: ENTROPY_MODE_LABELS[m],
        label_visibility="collapsed",
    )
    st.markdown(
        f'<div class="caption">{ENTROPY_MODE_HELP[mode]}</div>',
        unsafe_allow_html=True,
    )

    # ── Bin parameters ──
    st.markdown("### 04 · bin parameters")
    _section_caption(SEC_FP_BIN_PARAMS)
    col_k, col_cov, col_t, col_th = st.columns(4)
    with col_k:
        cpgs = st.number_input(
            "CpGs per bin (k)",
            min_value=2, max_value=8, value=int(fp.cpgs_per_bin), step=1,
            help=GUI_HELP_K,
        )
    with col_cov:
        min_cov = st.number_input(
            "min coverage", min_value=1, max_value=500,
            value=int(fp.min_coverage), step=1,
            help=GUI_HELP_MIN_COV,
        )
    with col_t:
        threads = st.slider(
            "threads", 1, 64, fp.modkit_threads,
            help="Worker processes for the per-chromosome stage.",
        )
    with col_th:
        thresh = st.slider(
            "methyl threshold", 0.0, 1.0,
            float(fp.methyl_threshold), step=0.05,
            help=GUI_HELP_METHYL_THRESHOLD,
        )
    update_section(
        "file_prep",
        entropy_mode=mode, cpgs_per_bin=int(cpgs),
        min_coverage=int(min_cov), modkit_threads=threads,
        methyl_threshold=float(thresh),
    )

    if mode == ENTROPY_MODE_TERNARY:
        rec = required_coverage_for_k(int(cpgs))
        if int(min_cov) < rec:
            st.warning(
                f"ternary at k={int(cpgs)} ⇒ 3^k = {rec} possible patterns. "
                f"min coverage {int(min_cov)} is below 3^k — results will be "
                f"noisy. recommended ≥ {rec}× coverage."
            )
        else:
            st.info(f"ternary 3^k = {rec}× coverage requirement satisfied.")

    # ── Output ──
    st.markdown("### 05 · output")
    _section_caption(SEC_FP_OUTPUT)
    out_dir = st.text_input(
        "output directory",
        value=str(fp.output_dir or RESULTS_DIR),
    )
    update_section("file_prep", output_dir=Path(out_dir) if out_dir else None)

    st.markdown("---")
    col_run, col_status = st.columns([1, 3])
    with col_run:
        run_clicked = st.button("RUN PIPELINE", type="primary",
                                use_container_width=True)
    with col_status:
        if fp.last_run_id:
            st.markdown(
                f'<div class="caption">last run · {fp.last_run_id}</div>',
                unsafe_allow_html=True,
            )

    if run_clicked:
        if not out_dir:
            st.error("set an output directory first")
            return
        summary = _run(Path(out_dir))
        if summary:
            st.success(
                "pipeline complete · output paths auto-loaded into the "
                "Graph Preparation tab"
            )
            with st.expander("run summary (JSON)"):
                st.json(summary)
