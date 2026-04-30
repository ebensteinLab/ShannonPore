"""Integration tests for `shannonpore doctor` and `shannonpore selftest`.

Both subcommands are entrypoints to the post-install verification flow.
We exercise them as the user would after running `bash install.sh`.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

V4_DIR = Path(__file__).resolve().parents[2]


def _run_cli(*args: str, **kw) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "src.cli", *args],
        cwd=str(V4_DIR),
        capture_output=True,
        text=True,
        **kw,
    )


@pytest.mark.integration
def test_doctor_runs_and_reports() -> None:
    proc = _run_cli("doctor")
    # On hosts where pinned versions don't match exactly, doctor exits 1 —
    # that's expected. The point is that it RAN and produced a coherent
    # report.
    assert proc.returncode in (0, 1)
    out = proc.stdout
    for required_row in ("python", "streamlit", "pandas", "numpy", "modkit", "RESULTS_DIR"):
        assert required_row in out, f"doctor row missing: {required_row}"
    # Final summary line — one of: all passed, advisories-only (still OK),
    # or hard failures.
    assert "All checks passed" in out or "Install is functional" in out or "check(s) failed" in out


@pytest.mark.integration
def test_doctor_help() -> None:
    proc = _run_cli("doctor", "--help")
    assert proc.returncode == 0
    out = proc.stdout.lower()
    # `--help` shows the long description; check several phrases that
    # appear in either the short help= or description strings.
    assert "checks python" in out or "after install" in out


@pytest.mark.integration
@pytest.mark.skipif(shutil.which("modkit") is None, reason="modkit not on PATH")
def test_selftest_runs_synthetic_pipeline_end_to_end() -> None:
    proc = _run_cli("selftest")
    assert proc.returncode == 0, (
        f"selftest failed.\nstdout:\n{proc.stdout[-1000:]}\n" f"stderr:\n{proc.stderr[-1000:]}"
    )
    assert "selftest passed" in proc.stdout.lower()
    assert "true_mc" in proc.stdout
    assert "bisulfite" in proc.stdout


@pytest.mark.integration
def test_selftest_help_lists_include_ternary() -> None:
    proc = _run_cli("selftest", "--help")
    assert proc.returncode == 0
    assert "--include-ternary" in proc.stdout
