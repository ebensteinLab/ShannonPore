"""Centralized session state for the v4 Streamlit app.

v3 used 107+ direct `st.session_state[...]` accesses with no schema. v4
defines every state key in `AppState` (a frozen dataclass for shape + types)
and mediates reads/writes through `SessionStore`. This keeps the session
schema discoverable, typed, and migrate-able.

Streamlit recreates the session_state dict on every rerun, so we keep ONE
key (`_app_state`) that holds the dataclass. Mutations replace the
dataclass with a new instance using `dataclasses.replace`.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from src.constants import PALETTE_CONTROL, PALETTE_TARGET


@dataclass
class SampleSpec:
    """One sample's input source. Either a single BAM, a single TSV, or
    a folder of BAMs to merge+sort+index into a single BAM."""

    label: str = "sample"
    input_kind: str = "bam"  # "bam" | "tsv" | "bam_folder"
    bam_path: Path | None = None
    tsv_path: Path | None = None
    bam_folder: Path | None = None


@dataclass
class FilePrepState:
    """Tab 1 — File Preparation.

    `pair_mode=False`: process `target` as a single sample.
    `pair_mode=True`:  process both `control` and `target`, producing
                       two parallel sets of bedgraphs.
    """

    target: SampleSpec = field(default_factory=lambda: SampleSpec(label="target"))
    control: SampleSpec = field(default_factory=lambda: SampleSpec(label="control"))
    pair_mode: bool = False
    genome: str = "hg38"
    custom_fasta: Path | None = None
    roi_text: str = ""
    roi_bed: Path | None = None
    use_whole_genome: bool = True
    modkit_threads: int = 8
    modkit_queue_size: int = 1000
    modkit_extra_args: str = ""
    include_bed: Path | None = None
    output_dir: Path | None = None
    cpgs_per_bin: int = 4
    min_coverage: int = 16
    methyl_threshold: float = 0.5
    entropy_mode: str = "true_mc"  # "true_mc" | "bisulfite" | "ternary"
    last_run_id: str | None = None
    last_run_summary: dict[str, Any] = field(default_factory=dict)


@dataclass
class TrackPlot:
    """One track-plot configuration.

    The default ``color`` matches the control side of the new palette;
    GraphPrepState's factory overrides target with PALETTE_TARGET.
    """

    name: str = ""
    color: str = PALETTE_CONTROL
    mml_path: Path | None = None
    me_path: Path | None = None
    coverage_path: Path | None = None


@dataclass
class GraphPrepState:
    """Tab 2 — Graph Preparation."""

    control: TrackPlot = field(
        default_factory=lambda: TrackPlot(name="Control", color=PALETTE_CONTROL),
    )
    target: TrackPlot = field(
        default_factory=lambda: TrackPlot(name="Target", color=PALETTE_TARGET),
    )
    genome: str = "hg38"
    gtf_path: Path | None = None  # custom override; None ⇒ use bundled gtf for `genome`
    region_chrom: str = ""
    region_start: int = 0
    region_end: int = 0
    plots_dir: Path | None = None
    scatter_subsample: int = 50_000
    scatter_color_by_density: bool = True


@dataclass
class AppState:
    """Top-level session state. Mutate via `SessionStore.update()`."""

    file_prep: FilePrepState = field(default_factory=FilePrepState)
    graph_prep: GraphPrepState = field(default_factory=GraphPrepState)
    sidebar_collapsed: bool = False
    last_error: str | None = None


# ─── Streamlit accessor ──────────────────────────────────────────────────

_STATE_KEY = "_app_state"


def get_state() -> AppState:
    """Return the single AppState instance for this session.

    Lazily initialises on first call. Imports streamlit here so non-UI tests
    can import this module without Streamlit installed.
    """
    import streamlit as st

    if _STATE_KEY not in st.session_state:
        st.session_state[_STATE_KEY] = AppState()
    return st.session_state[_STATE_KEY]


def update_state(**changes: Any) -> AppState:
    """Replace top-level fields on AppState. For nested updates use
    `update_section('file_prep', bam_path=...)`."""
    import streamlit as st

    new = replace(get_state(), **changes)
    st.session_state[_STATE_KEY] = new
    return new


def update_section(section: str, **changes: Any) -> AppState:
    """Replace fields on a sub-dataclass (file_prep, graph_prep)."""
    import streamlit as st

    state = get_state()
    if not hasattr(state, section):
        raise AttributeError(f"AppState has no section '{section}'")
    sub = replace(getattr(state, section), **changes)
    new = replace(state, **{section: sub})
    st.session_state[_STATE_KEY] = new
    return new


def reset_state() -> None:
    """Reset the session state. Useful for tests + 'Clear all' button."""
    import streamlit as st

    st.session_state[_STATE_KEY] = AppState()
