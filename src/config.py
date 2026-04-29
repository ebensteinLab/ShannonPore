"""Path and runtime configuration resolved from environment variables.

All filesystem paths flow through this module. No source file should embed
absolute deployment-specific paths (e.g. lab-specific roots) — set the env
vars instead, or override the defaults below.
"""

from __future__ import annotations

import logging
import os
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


def _env_path(var: str, default: Path) -> Path:
    raw = os.environ.get(var)
    return Path(raw).expanduser().resolve() if raw else default


# Project root = directory containing this file's parent's parent.
PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent

REFERENCE_DIR: Path = _env_path(
    "SHANNONPORE_REF_DIR",
    PROJECT_ROOT / "reference_files",
)

RESULTS_DIR: Path = _env_path(
    "SHANNONPORE_RESULTS_DIR",
    PROJECT_ROOT / "results",
)


@dataclass(frozen=True)
class GenomeAssets:
    """Reference assets for one genome assembly."""

    name: str
    fasta: Path
    gtf_gz: Path
    gtf_url: str
    ccre_bb: Path | None = None


HG38 = GenomeAssets(
    name="hg38",
    fasta=REFERENCE_DIR / "hg38.fa",
    gtf_gz=REFERENCE_DIR / "hg38.ncbiRefSeq.gtf.gz",
    gtf_url=(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/"
        "hg38.ncbiRefSeq.gtf.gz"
    ),
    ccre_bb=REFERENCE_DIR / "hg38_encodeCcreCcreCombined.bb",
)

MM10 = GenomeAssets(
    name="mm10",
    fasta=REFERENCE_DIR / "mm10.fa",
    gtf_gz=REFERENCE_DIR / "mm10.ncbiRefSeq.gtf.gz",
    gtf_url=(
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/"
        "mm10.ncbiRefSeq.gtf.gz"
    ),
    ccre_bb=REFERENCE_DIR / "mm10_encodeCcreCombined.bb",
)

GENOMES: dict[str, GenomeAssets] = {"hg38": HG38, "mm10": MM10}


# DuckDB / chunked I/O defaults.
INGEST_CHUNK_ROWS: int = 1_000_000


def assets_for(genome: str) -> GenomeAssets:
    if genome not in GENOMES:
        raise ValueError(
            f"Unknown genome '{genome}'. Available: {sorted(GENOMES)}"
        )
    return GENOMES[genome]


def ensure_dirs() -> None:
    """Create RESULTS_DIR and REFERENCE_DIR if missing.

    REFERENCE_DIR will be created as an empty directory; population is
    handled by ``scripts/setup_references.sh`` (or the user's existing
    reference path via ``SHANNONPORE_REF_DIR``). ``reference_status``
    reports whether each expected genome's FASTA is actually present.
    """
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)


def _download_to(url: str, dest: Path, *, progress_cb: Callable[[float, str], None] | None = None) -> None:
    """Stream-download ``url`` to ``dest`` atomically (writes to ``.partial``
    then renames). Reports fractional progress via ``progress_cb`` if given.

    Network failures bubble up as ``urllib.error.URLError`` /
    ``ConnectionError`` so callers can decide whether to fall back.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".partial")
    if partial.exists():
        partial.unlink()

    logger.info("Downloading %s → %s", url, dest)
    with urllib.request.urlopen(url, timeout=30) as resp:  # noqa: S310 (trusted UCSC URL)
        total_bytes = int(resp.headers.get("Content-Length", "0") or 0)
        read = 0
        chunk = 1 << 16
        with open(partial, "wb") as f:
            while True:
                buf = resp.read(chunk)
                if not buf:
                    break
                f.write(buf)
                read += len(buf)
                if progress_cb and total_bytes > 0:
                    progress_cb(read / total_bytes, f"{read / 1e6:.1f} / {total_bytes / 1e6:.1f} MB")
                elif progress_cb:
                    progress_cb(0.0, f"{read / 1e6:.1f} MB")
    partial.rename(dest)


def ensure_genome_gtf(
    genome: str, *, progress_cb: Callable[[float, str], None] | None = None,
) -> Path:
    """Return the path to the GTF for ``genome``, downloading it on first use.

    GTFs are small (~30–40 MB) so we lazy-download them transparently — the
    user never has to remember to run ``setup_references.sh`` for the gene
    panel to work. FASTAs are still opt-in via the install script because
    they're multi-GB.

    Raises ``RuntimeError`` if the download fails (network down, mirror
    moved, etc.); callers can catch this and fall back to "no gene panel".
    """
    assets = assets_for(genome)
    if assets.gtf_gz.exists() and assets.gtf_gz.stat().st_size > 0:
        return assets.gtf_gz

    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    try:
        _download_to(assets.gtf_url, assets.gtf_gz, progress_cb=progress_cb)
    except Exception as exc:  # noqa: BLE001 — surface as a single retryable error
        # Clean up any half-written file so the next attempt re-tries.
        partial = assets.gtf_gz.with_suffix(assets.gtf_gz.suffix + ".partial")
        for p in (partial, assets.gtf_gz):
            if p.exists() and p.stat().st_size == 0:
                p.unlink()
        raise RuntimeError(
            f"Failed to download {assets.gtf_url}: {exc}. "
            "Check internet access or run "
            f"`bash scripts/setup_references.sh {genome}` manually."
        ) from exc
    return assets.gtf_gz


def reference_status() -> dict[str, dict[str, bool]]:
    """Return per-genome presence map for the bundled reference assets.

    Used by ``shannonpore doctor`` and the GUI sidebar to surface the
    "you still need to download FASTAs" state actionably rather than as
    an opaque "missing dir" warning.
    """
    out: dict[str, dict[str, bool]] = {}
    for name, assets in GENOMES.items():
        out[name] = {
            "fasta": assets.fasta.exists(),
            "fasta_index": Path(str(assets.fasta) + ".fai").exists(),
            "gtf_gz": assets.gtf_gz.exists(),
        }
    return out
