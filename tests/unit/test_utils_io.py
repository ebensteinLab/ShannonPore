"""Unit tests for filesystem helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from shannonpore.io.utils_io import (
    disk_free_gb,
    ensure_writable_dir,
    safe_mkdir,
    validate_existing_path,
)


@pytest.mark.unit
def test_validate_existing_path_true_when_exists(tmp_path: Path) -> None:
    p = tmp_path / "x.txt"
    p.write_text("ok")
    assert validate_existing_path(str(p), "X") is True


@pytest.mark.unit
def test_validate_existing_path_returns_message_when_missing(tmp_path: Path) -> None:
    res = validate_existing_path(str(tmp_path / "no"), "X")
    assert isinstance(res, str)
    assert "does not exist" in res


@pytest.mark.unit
def test_validate_existing_path_empty_input() -> None:
    res = validate_existing_path("", "X")
    assert "empty" in str(res).lower()


@pytest.mark.unit
def test_safe_mkdir_creates(tmp_path: Path) -> None:
    p = tmp_path / "a" / "b"
    out = safe_mkdir(str(p))
    assert Path(out).exists()


@pytest.mark.unit
def test_ensure_writable_dir_writable(tmp_path: Path) -> None:
    out = ensure_writable_dir(str(tmp_path / "out"), "Out")
    assert Path(out).exists()


@pytest.mark.unit
def test_disk_free_gb_returns_positive_for_existing_path(tmp_path: Path) -> None:
    free = disk_free_gb(str(tmp_path))
    assert free > 0


@pytest.mark.unit
def test_disk_free_gb_returns_neg1_for_bogus_path() -> None:
    assert disk_free_gb("/no/such/path/should/exist/abc123xyz") == -1.0
