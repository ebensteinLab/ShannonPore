"""Shared matplotlib/plotly theme."""

from __future__ import annotations

import matplotlib.pyplot as plt


def apply_default_style() -> None:
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titlesize": 13,
        "axes.labelsize": 12,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })
