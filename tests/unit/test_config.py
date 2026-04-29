"""Unit tests for env-var-resolved config."""

from __future__ import annotations

import importlib
import os
from pathlib import Path

import pytest


@pytest.mark.unit
def test_env_var_overrides_default_reference_dir(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    fake_ref = tmp_path / "ref"
    fake_ref.mkdir()
    monkeypatch.setenv("NANOENTROPY_REF_DIR", str(fake_ref))
    monkeypatch.setenv(
        "NANOENTROPY_RESULTS_DIR", str(tmp_path / "results")
    )

    import src.config as cfg
    importlib.reload(cfg)

    assert cfg.REFERENCE_DIR == fake_ref.resolve()
    assert "results" in str(cfg.RESULTS_DIR)


@pytest.mark.unit
def test_assets_for_unknown_genome_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    import src.config as cfg
    with pytest.raises(ValueError, match="Unknown genome"):
        cfg.assets_for("dm6")


@pytest.mark.unit
def test_genomes_registry_has_hg38_and_mm10() -> None:
    import src.config as cfg
    assert "hg38" in cfg.GENOMES
    assert "mm10" in cfg.GENOMES
    assert cfg.GENOMES["hg38"].name == "hg38"


@pytest.mark.unit
def test_no_hardcoded_home_ebensteinlab_paths_in_source() -> None:
    """v4 must never embed `/home/EbensteinLab/...` literals."""
    src_dir = Path(__file__).resolve().parents[2] / "src"
    offenders: list[str] = []
    for py in src_dir.rglob("*.py"):
        text = py.read_text()
        if "/home/EbensteinLab/" in text:
            offenders.append(str(py))
    assert not offenders, (
        "Hardcoded /home/EbensteinLab/ paths found in: " + ", ".join(offenders)
    )
