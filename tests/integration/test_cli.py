"""Integration tests for the nanoentropy CLI.

Invokes ``python -m src.cli`` (rather than the installed `nanoentropy`
script) so the tests pass without `pip install -e .` having been run.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from src import __version__

V4_DIR = Path(__file__).resolve().parents[2]

SUBCOMMANDS = (
    "extract", "entropy", "plot", "run", "doctor", "guide", "examples", "selftest",
)


def _run(*args: str, expect_exit: int = 0, **kw) -> subprocess.CompletedProcess:
    proc = subprocess.run(
        [sys.executable, "-m", "src.cli", *args],
        cwd=str(V4_DIR), capture_output=True, text=True, **kw,
    )
    if proc.returncode != expect_exit:
        raise AssertionError(
            f"CLI exited {proc.returncode} (expected {expect_exit}).\n"
            f"args: {args}\nstdout: {proc.stdout}\nstderr: {proc.stderr}"
        )
    return proc


@pytest.mark.integration
def test_version_flag_matches_package_version() -> None:
    proc = _run("--version")
    assert __version__ in proc.stdout


@pytest.mark.integration
def test_top_level_help_lists_all_subcommands() -> None:
    proc = _run("--help")
    for sub in SUBCOMMANDS:
        assert sub in proc.stdout, f"subcommand {sub!r} missing from --help"


@pytest.mark.integration
@pytest.mark.parametrize("sub", SUBCOMMANDS)
def test_each_subcommand_has_help(sub: str) -> None:
    proc = _run(sub, "--help")
    assert "usage:" in proc.stdout.lower()


@pytest.mark.integration
def test_entropy_modes_appear_in_help() -> None:
    proc = _run("entropy", "--help")
    for mode in ("true_mc", "bisulfite", "ternary"):
        assert mode in proc.stdout


@pytest.mark.integration
def test_run_requires_input(tmp_path: Path) -> None:
    proc = _run(
        "run", "--out-dir", str(tmp_path), "--genome", "hg38",
        expect_exit=2,
    )
    assert any(x in proc.stderr for x in ("--bam", "--tsv", "--bam-folder"))


@pytest.mark.integration
def test_run_pair_requires_both_sides(tmp_path: Path) -> None:
    proc = _run(
        "run", "--pair", "--out-dir", str(tmp_path),
        "--control-bam", "/no/ctrl.bam",  # missing target
        "--genome", "hg38",
        expect_exit=2,
    )
    assert "control" in proc.stderr.lower() or "target" in proc.stderr.lower()


@pytest.mark.integration
def test_guide_prints_walkthrough() -> None:
    proc = _run("guide")
    assert "quick start" in proc.stdout.lower()
    assert "true_mc" in proc.stdout
    assert "ternary" in proc.stdout
    assert "--bam-folder" in proc.stdout


@pytest.mark.integration
def test_examples_prints_recipes() -> None:
    proc = _run("examples")
    assert "nanoentropy run" in proc.stdout
    assert "--pair" in proc.stdout
    assert "--bam-folder" in proc.stdout


@pytest.mark.integration
def test_entropy_rejects_unknown_mode(tmp_path: Path) -> None:
    proc = _run(
        "entropy", "/no/such.tsv", str(tmp_path / "out"),
        "--mode", "not-a-mode",
        expect_exit=2,
    )
    assert "invalid choice" in proc.stderr.lower()
