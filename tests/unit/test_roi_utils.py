"""Unit tests for ROI parsers."""

from __future__ import annotations

from pathlib import Path

import pytest

from shannonpore.io.roi_utils import parse_roi_file_flexible, parse_roi_text_area


@pytest.mark.unit
def test_text_area_basic_three_columns() -> None:
    df = parse_roi_text_area("chr1\t100\t200\nchr2\t300\t400\n")
    assert df is not None
    assert list(df.columns) == ["chr", "start", "end", "target"]
    assert len(df) == 2
    assert df.iloc[0]["target"] == "chr1:100-200"


@pytest.mark.unit
def test_text_area_with_target_name() -> None:
    df = parse_roi_text_area("chr1\t100\t200\tFOO\n")
    assert df is not None
    assert df.iloc[0]["target"] == "FOO"


@pytest.mark.unit
def test_text_area_skips_blank_and_comments() -> None:
    df = parse_roi_text_area("# comment\n\nchr1\t1\t2\n")
    assert df is not None
    assert len(df) == 1


@pytest.mark.unit
def test_text_area_empty_returns_none() -> None:
    assert parse_roi_text_area("") is None
    assert parse_roi_text_area("   \n\n") is None


@pytest.mark.unit
def test_text_area_invalid_coords_returns_none() -> None:
    assert parse_roi_text_area("chr1\t200\t100\n") is None  # end <= start
    assert parse_roi_text_area("chr1\tabc\t200\n") is None  # non-int


@pytest.mark.unit
def test_text_area_too_few_columns_returns_none() -> None:
    assert parse_roi_text_area("chr1\t100\n") is None


@pytest.mark.unit
def test_file_flexible_with_header(tmp_path: Path) -> None:
    p = tmp_path / "roi.bed"
    p.write_text("chr\tstart\tend\tname\nchr1\t100\t200\tA\n")
    df = parse_roi_file_flexible(str(p))
    assert df is not None
    assert list(df.columns) == ["chr", "start", "end", "target"]
    assert df.iloc[0]["target"] == "A"


@pytest.mark.unit
def test_file_flexible_without_header(tmp_path: Path) -> None:
    p = tmp_path / "roi.bed"
    p.write_text("chr1\t100\t200\nchr2\t300\t400\n")
    df = parse_roi_file_flexible(str(p))
    assert df is not None
    assert len(df) == 2


@pytest.mark.unit
def test_file_flexible_missing_file_returns_none(tmp_path: Path) -> None:
    assert parse_roi_file_flexible(str(tmp_path / "no_such.bed")) is None
