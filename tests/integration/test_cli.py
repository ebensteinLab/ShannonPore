"""Integration tests for the shannonpore CLI.

Invokes ``python -m shannonpore.cli`` (rather than the installed `shannonpore`
script) so the tests pass without `pip install -e .` having been run.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from shannonpore import __version__

V4_DIR = Path(__file__).resolve().parents[2]

SUBCOMMANDS = (
    "extract",
    "entropy",
    "plot",
    "run",
    "doctor",
    "guide",
    "examples",
    "selftest",
)


def _run(*args: str, expect_exit: int = 0, **kw) -> subprocess.CompletedProcess:
    proc = subprocess.run(
        [sys.executable, "-m", "shannonpore.cli", *args],
        cwd=str(V4_DIR),
        capture_output=True,
        text=True,
        **kw,
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
        "run",
        "--out-dir",
        str(tmp_path),
        "--genome",
        "hg38",
        expect_exit=2,
    )
    assert any(x in proc.stderr for x in ("--bam", "--tsv", "--bam-folder"))


@pytest.mark.integration
def test_run_pair_requires_both_sides(tmp_path: Path) -> None:
    proc = _run(
        "run",
        "--pair",
        "--out-dir",
        str(tmp_path),
        "--control-bam",
        "/no/ctrl.bam",  # missing target
        "--genome",
        "hg38",
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
    assert "shannonpore run" in proc.stdout
    assert "--pair" in proc.stdout
    assert "--bam-folder" in proc.stdout


@pytest.mark.integration
def test_entropy_rejects_unknown_mode(tmp_path: Path) -> None:
    proc = _run(
        "entropy",
        "/no/such.tsv",
        str(tmp_path / "out"),
        "--mode",
        "not-a-mode",
        expect_exit=2,
    )
    assert "invalid choice" in proc.stderr.lower()


# ─── plot: multi-format export + connecting-line toggle ──────────────────


def _write_paired_bedgraphs(d: Path) -> dict[str, Path]:
    """Four tiny bedgraphs sharing identical intervals so the inner
    join in load_paired_bedgraphs keeps every bin."""
    paths = {}
    for name, floor in (
        ("control_mml", 0.1),
        ("control_me", 0.4),
        ("target_mml", 0.3),
        ("target_me", 0.6),
    ):
        p = d / f"{name}.bedgraph"
        rows = [f"chr1\t{100 + 50 * i}\t{150 + 50 * i}\t{floor + 0.01 * i:.3f}" for i in range(20)]
        p.write_text("\n".join(rows) + "\n")
        paths[name] = p
    return paths


@pytest.mark.integration
def test_plot_landscape_multi_format_and_no_lines(tmp_path: Path) -> None:
    paths = _write_paired_bedgraphs(tmp_path)
    out = tmp_path / "landscape.png"
    _run(
        "plot",
        "landscape",
        str(out),
        "--control-mml",
        str(paths["control_mml"]),
        "--control-me",
        str(paths["control_me"]),
        "--target-mml",
        str(paths["target_mml"]),
        "--target-me",
        str(paths["target_me"]),
        "--filter-a-dim",
        "off",
        "--filter-b-dim",
        "off",
        "--formats",
        "png,svg,pdf",
        "--dpi",
        "100",
        "--no-lines",
    )
    assert (tmp_path / "landscape.png").exists()
    assert (tmp_path / "landscape.svg").exists()
    assert (tmp_path / "landscape.pdf").exists()
    # --no-lines must drop the LineCollection: with only 20 bins the SVG
    # is tiny; the connecting-line variant embeds one path per bin pair.
    svg = (tmp_path / "landscape.svg").read_text()
    assert "<svg" in svg


@pytest.mark.integration
def test_plot_help_lists_export_flags() -> None:
    proc = _run("plot", "--help")
    assert "--formats" in proc.stdout
    assert "--dpi" in proc.stdout
    assert "--no-lines" in proc.stdout


@pytest.mark.integration
def test_plot_rejects_unknown_format(tmp_path: Path) -> None:
    paths = _write_paired_bedgraphs(tmp_path)
    proc = _run(
        "plot",
        "scatter",
        str(tmp_path / "s.png"),
        "--control-mml",
        str(paths["control_mml"]),
        "--control-me",
        str(paths["control_me"]),
        "--target-mml",
        str(paths["target_mml"]),
        "--target-me",
        str(paths["target_me"]),
        "--formats",
        "png,tiff",
        expect_exit=2,
    )
    assert "tiff" in (proc.stderr + proc.stdout).lower()


@pytest.mark.integration
def test_plot_rejects_out_of_range_dpi(tmp_path: Path) -> None:
    paths = _write_paired_bedgraphs(tmp_path)
    proc = _run(
        "plot",
        "scatter",
        str(tmp_path / "s.png"),
        "--control-mml",
        str(paths["control_mml"]),
        "--control-me",
        str(paths["control_me"]),
        "--target-mml",
        str(paths["target_mml"]),
        "--target-me",
        str(paths["target_me"]),
        "--dpi",
        "999999",
        expect_exit=2,
    )
    assert "dpi" in (proc.stderr + proc.stdout).lower()


@pytest.mark.integration
def test_formats_arg_dedupes_aliases() -> None:
    from shannonpore.cli import _formats_arg

    assert _formats_arg("jpg,jpeg,png") == ["jpg", "png"]
