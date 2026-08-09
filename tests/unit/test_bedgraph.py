"""Unit tests for bedgraph reader/writer/smoother."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from shannonpore.io.bedgraph import read_bedgraph, smooth_bedgraph, write_bedgraph


@pytest.mark.unit
def test_read_basic(tiny_bedgraph: Path) -> None:
    df = read_bedgraph(tiny_bedgraph)
    assert list(df.columns) == ["chrom", "start", "end", "value"]
    assert len(df) == 5


@pytest.mark.unit
def test_read_with_chrom_filter(tiny_bedgraph: Path) -> None:
    df = read_bedgraph(tiny_bedgraph, chrom_filter="chr2")
    assert len(df) == 1
    assert df.iloc[0]["chrom"] == "chr2"


@pytest.mark.unit
def test_read_empty_file(tmp_path: Path) -> None:
    p = tmp_path / "empty.bedgraph"
    p.write_text("")
    df = read_bedgraph(p)
    assert df.empty


@pytest.mark.unit
def test_write_round_trip(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 10],
            "end": [10, 20],
            "metric": [0.5, 0.7],
        }
    )
    out = tmp_path / "rt.bedgraph"
    write_bedgraph(df, "metric", out)
    rt = read_bedgraph(out)
    assert list(rt["value"]) == [0.5, 0.7]


@pytest.mark.unit
def test_write_empty_produces_empty_file(tmp_path: Path) -> None:
    out = tmp_path / "empty.bedgraph"
    write_bedgraph(pd.DataFrame(), "value", out)
    assert out.exists()
    assert out.read_text() == ""


@pytest.mark.unit
def test_smooth_bedgraph_returns_smoothed_score(tiny_bedgraph: Path) -> None:
    df = smooth_bedgraph(tiny_bedgraph, chrom_filter="chr1", window_size=2)
    assert "smoothed_score" in df.columns
    assert "mid" in df.columns
    assert len(df) == 4
