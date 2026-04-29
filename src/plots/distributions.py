"""Distribution plots: violin, box, histogram of per-bin metrics."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def metric_distribution(
    df: pd.DataFrame,
    value_col: str,
    *,
    by: str | None = None,
    kind: str = "violin",
    out_path: Path | str | None = None,
) -> plt.Figure:
    """Distribution plot of `value_col`. If `by` given, group / facet by it."""
    fig, ax = plt.subplots(figsize=(7, 5))
    if kind == "violin":
        sns.violinplot(data=df, x=by, y=value_col, ax=ax, inner="quartile")
    elif kind == "box":
        sns.boxplot(data=df, x=by, y=value_col, ax=ax)
    elif kind == "hist":
        if by:
            for g, sub in df.groupby(by):
                ax.hist(sub[value_col].dropna(), bins=50, alpha=0.5, label=str(g))
            ax.legend()
        else:
            ax.hist(df[value_col].dropna(), bins=50)
        ax.set_xlabel(value_col)
        ax.set_ylabel("count")
    else:
        raise ValueError(f"Unknown distribution kind: {kind!r}")

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=200)
    return fig
