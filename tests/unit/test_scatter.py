"""Unit tests for the comparison plots in src.plots.scatter.

Renders are exercised on synthetic DataFrames; we don't compare pixels,
just sanity-check the figure shape and that "fallback" code paths don't
crash on edge inputs.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from src.plots.scatter import (
    me_mml_scatter,
    paired_landscape,
    theoretical_entropy_binary,
    triple_landscape,
)


@pytest.fixture()
def paired_df() -> pd.DataFrame:
    rng = np.random.default_rng(42)
    n = 200
    return pd.DataFrame({
        "mml_a": rng.uniform(0, 1, n),
        "me_a": rng.uniform(0, 1, n),
        "mml_b": rng.uniform(0, 1, n),
        "me_b": rng.uniform(0, 1, n),
    })


@pytest.mark.unit
def test_theoretical_entropy_binary_endpoints_zero() -> None:
    p = np.array([0.0, 1.0])
    h = theoretical_entropy_binary(p)
    assert np.allclose(h, 0.0)


@pytest.mark.unit
def test_theoretical_entropy_binary_max_at_half() -> None:
    h = theoretical_entropy_binary(np.array([0.5]))
    assert pytest.approx(float(h[0]), abs=1e-9) == 1.0


@pytest.mark.unit
def test_me_mml_scatter_renders_two_panels(paired_df) -> None:
    fig = me_mml_scatter(paired_df, label_a="A", label_b="B", log_scale=True)
    # 2 image axes + 2 colorbar axes
    assert len(fig.axes) == 4
    plt.close(fig)


@pytest.mark.unit
def test_me_mml_scatter_handles_empty_df() -> None:
    df = pd.DataFrame(
        {"mml_a": [], "me_a": [], "mml_b": [], "me_b": []},
        dtype=float,
    )
    fig = me_mml_scatter(df, label_a="A", label_b="B")
    plt.close(fig)


@pytest.mark.unit
def test_triple_landscape_renders_three_panels(paired_df) -> None:
    fig = triple_landscape(
        paired_df, label_a="A", label_b="B", log_scale=False,
    )
    assert len(fig.axes) >= 3
    plt.close(fig)


@pytest.mark.unit
def test_paired_landscape_no_filter(paired_df) -> None:
    fig = paired_landscape(
        paired_df, label_a="A", label_b="B",
        filter_a_dim=None, filter_b_dim=None,
    )
    plt.close(fig)


@pytest.mark.unit
def test_paired_landscape_no_bins_pass_filter() -> None:
    """Filter forces zero-survivor branch — must not crash."""
    df = pd.DataFrame({
        "mml_a": [0.5], "me_a": [0.5],
        "mml_b": [0.5], "me_b": [0.5],
    })
    fig = paired_landscape(
        df, label_a="A", label_b="B",
        filter_a_dim="|dMML|", filter_a_op=">", filter_a_value=0.9,
        filter_b_dim=None,
    )
    # The fallback text should be drawn on the only axes.
    text_blobs = [t.get_text() for t in fig.axes[0].texts]
    assert any("no bins" in t.lower() for t in text_blobs)
    plt.close(fig)


@pytest.mark.unit
def test_paired_landscape_subsamples_when_over_max(paired_df) -> None:
    """max_lines smaller than n_bins must keep counts on full set but
    only draw `max_lines` segments — no crash, caption updates."""
    fig = paired_landscape(
        paired_df, label_a="A", label_b="B",
        filter_a_dim=None, filter_b_dim=None,
        max_lines=10,  # df has 200
    )
    plt.close(fig)
