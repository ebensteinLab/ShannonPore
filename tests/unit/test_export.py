"""Unit tests for the multi-format figure export helper.

Checks real file output per format (magic bytes, not pixels), DPI
scaling for raster formats, and input validation.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless

from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from src.plots.export import SUPPORTED_FORMATS, save_figure

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
JPG_MAGIC = b"\xff\xd8\xff"
PDF_MAGIC = b"%PDF"


@pytest.fixture()
def fig() -> plt.Figure:
    f, ax = plt.subplots(figsize=(4, 3))
    ax.plot([0, 1], [0, 1])
    yield f
    plt.close(f)


@pytest.mark.unit
def test_supported_formats_exact() -> None:
    assert set(SUPPORTED_FORMATS) == {"png", "jpg", "svg", "pdf"}


@pytest.mark.unit
def test_save_single_path_no_formats(fig, tmp_path: Path) -> None:
    """formats=None → save exactly out_path, format from its suffix."""
    out = tmp_path / "plot.png"
    written = save_figure(fig, out)
    assert written == [out]
    assert out.read_bytes()[:8] == PNG_MAGIC


@pytest.mark.unit
def test_save_all_formats(fig, tmp_path: Path) -> None:
    out = tmp_path / "plot.png"
    written = save_figure(fig, out, formats=("png", "jpg", "svg", "pdf"))
    assert [p.suffix for p in written] == [".png", ".jpg", ".svg", ".pdf"]
    by_suffix = {p.suffix: p for p in written}
    assert by_suffix[".png"].read_bytes()[:8] == PNG_MAGIC
    assert by_suffix[".jpg"].read_bytes()[:3] == JPG_MAGIC
    assert by_suffix[".pdf"].read_bytes()[:4] == PDF_MAGIC
    assert b"<svg" in by_suffix[".svg"].read_bytes()[:2048]


@pytest.mark.unit
def test_formats_share_stem_with_out_path(fig, tmp_path: Path) -> None:
    out = tmp_path / "myplot.png"
    written = save_figure(fig, out, formats=("svg", "pdf"))
    assert all(p.stem == "myplot" for p in written)
    assert all(p.parent == tmp_path for p in written)


@pytest.mark.unit
def test_dpi_scales_raster_output(fig, tmp_path: Path) -> None:
    from PIL import Image

    lo = save_figure(fig, tmp_path / "lo.png", dpi=100)[0]
    hi = save_figure(fig, tmp_path / "hi.png", dpi=300)[0]
    with Image.open(lo) as im_lo, Image.open(hi) as im_hi:
        assert im_hi.width > im_lo.width * 2
        assert im_hi.height > im_lo.height * 2


@pytest.mark.unit
def test_unknown_format_raises(fig, tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="tiff"):
        save_figure(fig, tmp_path / "plot.png", formats=("png", "tiff"))


@pytest.mark.unit
def test_unlisted_suffix_passes_through_to_matplotlib(fig, tmp_path: Path) -> None:
    """Back-compat: without ``formats``, any matplotlib-supported suffix
    (e.g. .eps) is honoured exactly as the old direct-savefig code did."""
    out = tmp_path / "plot.eps"
    written = save_figure(fig, out)
    assert written == [out]
    assert out.stat().st_size > 0


@pytest.mark.unit
def test_bare_string_format_accepted(fig, tmp_path: Path) -> None:
    """A single format passed as a bare str must not be iterated
    character-by-character."""
    written = save_figure(fig, tmp_path / "plot.png", formats="svg")
    assert [p.suffix for p in written] == [".svg"]


@pytest.mark.unit
def test_creates_missing_parent_dir(fig, tmp_path: Path) -> None:
    out = tmp_path / "deep" / "nested" / "plot.svg"
    written = save_figure(fig, out, formats=("svg",))
    assert written[0].exists()


@pytest.mark.unit
def test_empty_formats_raises(fig, tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        save_figure(fig, tmp_path / "plot.png", formats=())


@pytest.mark.unit
def test_jpeg_extension_normalised(fig, tmp_path: Path) -> None:
    """'jpeg' is accepted as an alias and written as .jpg."""
    written = save_figure(fig, tmp_path / "plot.png", formats=("jpeg",))
    assert written[0].suffix == ".jpg"
    assert written[0].read_bytes()[:3] == JPG_MAGIC


@pytest.mark.unit
@pytest.mark.parametrize("bad_dpi", [0, -50, 100_000])
def test_out_of_range_dpi_raises(fig, tmp_path: Path, bad_dpi: int) -> None:
    with pytest.raises(ValueError, match="dpi"):
        save_figure(fig, tmp_path / "plot.png", dpi=bad_dpi)
