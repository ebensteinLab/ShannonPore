"""End-to-end smoke tests using `streamlit.testing.v1.AppTest`.

These tests exercise the actual Streamlit script without spinning up a
browser. They confirm that the page renders, all 3 tabs initialise without
errors, and headline widgets are present.
"""

from __future__ import annotations

import pytest

streamlit_testing = pytest.importorskip("streamlit.testing.v1")
AppTest = streamlit_testing.AppTest


pytestmark = pytest.mark.e2e


def test_app_renders_without_error() -> None:
    at = AppTest.from_file("app.py").run(timeout=30)
    assert not at.exception, f"App raised: {at.exception}"


def test_two_tabs_present() -> None:
    at = AppTest.from_file("app.py").run(timeout=30)
    # Two tabs after data-analysis was removed: File Preparation, Graph Preparation
    assert len(at.tabs) == 2


def test_sidebar_contains_reset_button() -> None:
    at = AppTest.from_file("app.py").run(timeout=30)
    button_labels = [b.label.lower() for b in at.sidebar.button]
    assert any("reset" in lbl for lbl in button_labels)


def test_appstate_initialised_in_session() -> None:
    at = AppTest.from_file("app.py").run(timeout=30)
    assert "_app_state" in at.session_state
    state = at.session_state["_app_state"]
    # Default: single-sample mode, target spec defaults to BAM input.
    assert state.file_prep.pair_mode is False
    assert state.file_prep.target.input_kind == "bam"
    assert state.file_prep.control.input_kind == "bam"
