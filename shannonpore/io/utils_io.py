"""Filesystem utility helpers (mkdir, validation, disk-free)."""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)

PathLike = str | Path


def validate_existing_path(path_value: PathLike, description: str) -> bool | str:
    """Return True if path exists, otherwise an error message string."""
    if not path_value or str(path_value).strip() == "":
        return f"{description} path is empty."
    if not os.path.exists(path_value):
        return f"{description} path does not exist: {path_value}"
    return True


def ensure_writable_dir(path_value: PathLike, description: str) -> str:
    """Create the directory if missing and verify it is writable."""
    if not path_value or str(path_value).strip() == "":
        raise ValueError(f"{description} is empty.")
    abs_path = os.path.abspath(str(path_value))
    os.makedirs(abs_path, exist_ok=True)

    test_file = os.path.join(abs_path, ".write_test.tmp")
    try:
        with open(test_file, "w", encoding="utf-8") as f:
            f.write("ok")
        os.remove(test_file)
    except OSError as exc:
        raise ValueError(f"{description} is not writable: {abs_path}. Error: {exc}") from exc

    return abs_path


def safe_mkdir(path_value: PathLike) -> str:
    """Create a directory (and parents) and return its absolute path."""
    abs_path = os.path.abspath(str(path_value))
    os.makedirs(abs_path, exist_ok=True)
    return abs_path


def disk_free_gb(path_value: PathLike) -> float:
    """Return free disk space in GiB at `path_value`. -1.0 on error."""
    try:
        usage = shutil.disk_usage(str(path_value))
    except OSError as exc:
        logger.warning("disk_free_gb(%s) failed: %s", path_value, exc)
        return -1.0
    return float(usage.free) / (1024.0**3)
