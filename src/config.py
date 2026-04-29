"""Path and runtime configuration resolved from environment variables.

All filesystem paths flow through this module. No source file should embed
absolute deployment-specific paths (e.g. lab-specific roots) — set the env
vars instead, or override the defaults below.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


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
    ccre_bb: Path | None = None


HG38 = GenomeAssets(
    name="hg38",
    fasta=REFERENCE_DIR / "hg38.fa",
    gtf_gz=REFERENCE_DIR / "hg38.ncbiRefSeq.gtf.gz",
    ccre_bb=REFERENCE_DIR / "hg38_encodeCcreCcreCombined.bb",
)

MM10 = GenomeAssets(
    name="mm10",
    fasta=REFERENCE_DIR / "mm10.fa",
    gtf_gz=REFERENCE_DIR / "mm10.ncbiRefSeq.gtf.gz",
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
