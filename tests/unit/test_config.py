"""Unit tests for env-var-resolved config."""

from __future__ import annotations

import importlib
from pathlib import Path

import pytest


@pytest.mark.unit
def test_env_var_overrides_default_reference_dir(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    fake_ref = tmp_path / "ref"
    fake_ref.mkdir()
    monkeypatch.setenv("SHANNONPORE_REF_DIR", str(fake_ref))
    monkeypatch.setenv(
        "SHANNONPORE_RESULTS_DIR", str(tmp_path / "results")
    )

    import src.config as cfg
    importlib.reload(cfg)

    assert fake_ref.resolve() == cfg.REFERENCE_DIR
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


@pytest.mark.unit
def test_reference_status_returns_per_genome_dict() -> None:
    import src.config as cfg
    importlib.reload(cfg)
    status = cfg.reference_status()
    assert set(status) >= {"hg38", "mm10"}
    for v in status.values():
        assert {"fasta", "fasta_index", "gtf_gz"} == set(v)
        assert all(isinstance(b, bool) for b in v.values())


@pytest.mark.unit
def test_ensure_genome_gtf_short_circuits_when_file_present(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """If the GTF is already on disk, no download attempt should happen."""
    monkeypatch.setenv("SHANNONPORE_REF_DIR", str(tmp_path))
    monkeypatch.setenv("SHANNONPORE_RESULTS_DIR", str(tmp_path / "out"))
    import src.config as cfg
    importlib.reload(cfg)

    fake_gtf = cfg.GENOMES["hg38"].gtf_gz
    fake_gtf.parent.mkdir(parents=True, exist_ok=True)
    fake_gtf.write_bytes(b"\x1f\x8bplaceholder")  # non-empty

    called = {"n": 0}

    def fail_download(url, dest, *, progress_cb=None):  # noqa: ANN001
        called["n"] += 1
        raise RuntimeError("should not be called")

    monkeypatch.setattr(cfg, "_download_to", fail_download)
    out = cfg.ensure_genome_gtf("hg38")
    assert out == fake_gtf
    assert called["n"] == 0


@pytest.mark.unit
def test_ensure_genome_gtf_downloads_when_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    monkeypatch.setenv("SHANNONPORE_REF_DIR", str(tmp_path))
    monkeypatch.setenv("SHANNONPORE_RESULTS_DIR", str(tmp_path / "out"))
    import src.config as cfg
    importlib.reload(cfg)

    expected = cfg.GENOMES["mm10"].gtf_gz
    assert not expected.exists()

    def fake_download(url, dest, *, progress_cb=None):  # noqa: ANN001
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        Path(dest).write_bytes(b"DOWNLOADED")

    monkeypatch.setattr(cfg, "_download_to", fake_download)
    out = cfg.ensure_genome_gtf("mm10")
    assert out == expected
    assert expected.exists()
    assert expected.read_bytes() == b"DOWNLOADED"


@pytest.mark.unit
def test_ensure_genome_gtf_failure_raises_runtime_error_and_cleans_up(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    monkeypatch.setenv("SHANNONPORE_REF_DIR", str(tmp_path))
    monkeypatch.setenv("SHANNONPORE_RESULTS_DIR", str(tmp_path / "out"))
    import src.config as cfg
    importlib.reload(cfg)

    expected = cfg.GENOMES["mm10"].gtf_gz

    # Simulate a partial download that died mid-stream — file exists with
    # bytes but the download function still raised.
    def half_download(url, dest, *, progress_cb=None):  # noqa: ANN001
        Path(dest).parent.mkdir(parents=True, exist_ok=True)
        Path(dest).write_bytes(b"HALFWRITTEN")
        raise OSError("network died")

    monkeypatch.setattr(cfg, "_download_to", half_download)
    with pytest.raises(RuntimeError, match="Failed to download"):
        cfg.ensure_genome_gtf("mm10")
    # The half-written file must NOT survive — otherwise next ensure_*
    # call would happily return a corrupt GTF.
    assert not expected.exists()
