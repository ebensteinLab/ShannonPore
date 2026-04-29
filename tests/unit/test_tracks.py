"""Unit tests for the 4-panel region track plot.

Exercises both 3-panel (no coverage) and 4-panel (with coverage) paths
on synthetic bedgraphs without needing modkit or any real BAMs.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from pathlib import Path

import pytest

from src.plots.tracks import plot_region_tracks


def _write_bedgraph(p: Path, *, chrom: str, value_floor: float) -> None:
    rows = []
    for i in range(50):
        s = 100 + 50 * i
        e = s + 50
        v = value_floor + 0.005 * i
        rows.append(f"{chrom}\t{s}\t{e}\t{v:.3f}")
    p.write_text("\n".join(rows) + "\n")


@pytest.mark.unit
def test_plot_region_tracks_three_panels_no_coverage(tmp_path: Path) -> None:
    chrom = "chr1"
    ctrl_mml = tmp_path / "ctrl.mml.bedgraph"
    ctrl_me = tmp_path / "ctrl.me.bedgraph"
    tgt_mml = tmp_path / "tgt.mml.bedgraph"
    tgt_me = tmp_path / "tgt.me.bedgraph"
    for p, floor in (
        (ctrl_mml, 0.10), (ctrl_me, 0.40),
        (tgt_mml, 0.30), (tgt_me, 0.60),
    ):
        _write_bedgraph(p, chrom=chrom, value_floor=floor)

    out = tmp_path / "tracks.png"
    plot_region_tracks(
        chrom=chrom, start=200, end=2500,
        control_mml=ctrl_mml, target_mml=tgt_mml,
        control_me=ctrl_me, target_me=tgt_me,
        gtf_path=None,
        out_path=out,
    )
    assert out.exists() and out.stat().st_size > 0
    # Confirm it's a real PNG by header magic.
    assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


@pytest.mark.unit
def test_plot_region_tracks_four_panels_with_coverage(tmp_path: Path) -> None:
    chrom = "chr1"
    paths = {
        name: tmp_path / f"{name}.bedgraph"
        for name in ("ctrl_mml", "ctrl_me", "tgt_mml", "tgt_me",
                     "ctrl_cov", "tgt_cov")
    }
    for name, p in paths.items():
        floor = 5.0 if "cov" in name else 0.2
        _write_bedgraph(p, chrom=chrom, value_floor=floor)

    out = tmp_path / "tracks_cov.png"
    plot_region_tracks(
        chrom=chrom, start=200, end=2500,
        control_mml=paths["ctrl_mml"], target_mml=paths["tgt_mml"],
        control_me=paths["ctrl_me"], target_me=paths["tgt_me"],
        control_coverage=paths["ctrl_cov"],
        target_coverage=paths["tgt_cov"],
        gtf_path=None,
        out_path=out,
    )
    assert out.exists() and out.stat().st_size > 0


@pytest.mark.unit
def test_plot_region_tracks_handles_no_overlapping_bins(tmp_path: Path) -> None:
    """If the requested region misses every bedgraph bin, the function
    must still emit a valid PNG with a fallback message in each panel."""
    chrom = "chr1"
    ctrl_mml = tmp_path / "ctrl.mml.bedgraph"
    ctrl_me = tmp_path / "ctrl.me.bedgraph"
    tgt_mml = tmp_path / "tgt.mml.bedgraph"
    tgt_me = tmp_path / "tgt.me.bedgraph"
    for p in (ctrl_mml, ctrl_me, tgt_mml, tgt_me):
        _write_bedgraph(p, chrom=chrom, value_floor=0.2)

    out = tmp_path / "tracks_empty.png"
    plot_region_tracks(
        chrom=chrom,
        start=10_000_000, end=10_001_000,  # far outside the synthetic bins
        control_mml=ctrl_mml, target_mml=tgt_mml,
        control_me=ctrl_me, target_me=tgt_me,
        gtf_path=None, pad_bp=0,
        out_path=out,
    )
    assert out.exists() and out.stat().st_size > 0
