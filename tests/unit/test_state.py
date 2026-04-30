"""Unit tests for AppState dataclass."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from src.state import (
    AppState,
    FilePrepState,
    GraphPrepState,
    SampleSpec,
    TrackPlot,
)


@pytest.mark.unit
def test_appstate_defaults_are_isolated_per_instance() -> None:
    a = AppState()
    b = AppState()
    assert a.file_prep is not b.file_prep
    assert a.graph_prep.control is not b.graph_prep.control


@pytest.mark.unit
def test_replace_keeps_other_fields_unchanged() -> None:
    state = AppState()
    new_target = replace(
        state.file_prep.target,
        bam_path=Path("/tmp/foo.bam"),
    )
    new = replace(
        state,
        file_prep=replace(state.file_prep, target=new_target),
    )
    assert new.file_prep.target.bam_path == Path("/tmp/foo.bam")
    assert new.graph_prep == state.graph_prep
    assert state.file_prep.target.bam_path is None


@pytest.mark.unit
def test_track_plot_defaults() -> None:
    from src.constants import PALETTE_CONTROL

    t = TrackPlot()
    # Default colour is the central palette entry, not a hard-coded hex,
    # so a future palette change in src.constants flows through here.
    assert t.color == PALETTE_CONTROL
    assert t.mml_path is None and t.me_path is None and t.coverage_path is None


@pytest.mark.unit
def test_sample_spec_defaults() -> None:
    s = SampleSpec()
    assert s.input_kind == "bam"
    assert s.bam_path is None and s.tsv_path is None and s.bam_folder is None


@pytest.mark.unit
def test_file_prep_defaults() -> None:
    fp = FilePrepState()
    assert fp.target.input_kind == "bam"
    assert fp.control.input_kind == "bam"
    assert fp.target.label == "target"
    assert fp.control.label == "control"
    assert fp.pair_mode is False
    assert fp.modkit_threads == 8


@pytest.mark.unit
def test_graph_prep_distinct_control_target_colors() -> None:
    gp = GraphPrepState()
    assert gp.control.color != gp.target.color
