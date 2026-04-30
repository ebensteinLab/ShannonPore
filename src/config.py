"""Path and runtime configuration resolved from environment variables.

All filesystem paths flow through this module. No source file should embed
absolute deployment-specific paths (e.g. lab-specific roots) — set the env
vars instead, or override the defaults below.
"""

from __future__ import annotations

import gzip
import logging
import os
import shutil
import subprocess
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
    fasta_url: str
    gtf_gz: Path
    gtf_url: str
    ccre_bb: Path | None = None


HG38 = GenomeAssets(
    name="hg38",
    fasta=REFERENCE_DIR / "hg38.fa",
    fasta_url=("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"),
    gtf_gz=REFERENCE_DIR / "hg38.ncbiRefSeq.gtf.gz",
    gtf_url=(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/" "hg38.ncbiRefSeq.gtf.gz"
    ),
    ccre_bb=REFERENCE_DIR / "hg38_encodeCcreCcreCombined.bb",
)

MM10 = GenomeAssets(
    name="mm10",
    fasta=REFERENCE_DIR / "mm10.fa",
    fasta_url=("https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"),
    gtf_gz=REFERENCE_DIR / "mm10.ncbiRefSeq.gtf.gz",
    gtf_url=(
        "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/" "mm10.ncbiRefSeq.gtf.gz"
    ),
    ccre_bb=REFERENCE_DIR / "mm10_encodeCcreCombined.bb",
)

GENOMES: dict[str, GenomeAssets] = {"hg38": HG38, "mm10": MM10}


# DuckDB / chunked I/O defaults.
INGEST_CHUNK_ROWS: int = 1_000_000


def assets_for(genome: str) -> GenomeAssets:
    if genome not in GENOMES:
        raise ValueError(f"Unknown genome '{genome}'. Available: {sorted(GENOMES)}")
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


def _download_to(
    url: str, dest: Path, *, progress_cb: Callable[[float, str], None] | None = None
) -> None:
    """Stream-download ``url`` to ``dest`` atomically (writes to ``.partial``
    then renames). Reports fractional progress via ``progress_cb`` if given.

    Hardened against:
      * stalled mid-transfer reads — sets a socket-level default timeout so
        ``resp.read()`` cannot hang indefinitely;
      * leftover ``.partial`` files from a prior crashed run.

    Network failures bubble up as ``urllib.error.URLError`` /
    ``ConnectionError`` / ``TimeoutError`` so callers can decide whether to
    fall back.
    """
    import socket

    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".partial")
    if partial.exists():
        partial.unlink()

    logger.info("Downloading %s → %s", url, dest)
    # urllib's `timeout=` covers connect + initial response, not stalled
    # mid-stream reads. Set the global socket timeout for the duration
    # of the read loop so a hung server can't block us forever.
    prior_default_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(60)
    try:
        with urllib.request.urlopen(
            url, timeout=30
        ) as resp:  # noqa: S310 (trusted UCSC URL, see GenomeAssets.gtf_url)
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
                        progress_cb(
                            read / total_bytes, f"{read / 1e6:.1f} / {total_bytes / 1e6:.1f} MB"
                        )
                    elif progress_cb:
                        progress_cb(0.0, f"{read / 1e6:.1f} MB")
        partial.rename(dest)
    finally:
        socket.setdefaulttimeout(prior_default_timeout)


def ensure_genome_gtf(
    genome: str,
    *,
    progress_cb: Callable[[float, str], None] | None = None,
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
        # Clean up partials AND any non-empty incomplete final file so the
        # next attempt starts fresh. Previously this only removed zero-size
        # files, which left a half-written GTF on disk after a stall.
        partial = assets.gtf_gz.with_suffix(assets.gtf_gz.suffix + ".partial")
        for p in (partial, assets.gtf_gz):
            try:
                if p.exists():
                    p.unlink()
            except OSError:
                logger.warning("Could not remove %s during cleanup", p)
        raise RuntimeError(
            f"Failed to download {assets.gtf_url}: {exc}. "
            "Check internet access or run "
            f"`bash scripts/setup_references.sh {genome}` manually."
        ) from exc
    return assets.gtf_gz


def _gunzip_to(
    src_gz: Path, dest: Path, *, progress_cb: Callable[[float, str], None] | None = None
) -> None:
    """Stream-decompress ``src_gz`` to ``dest`` atomically (writes to
    ``.partial`` then renames). Reports fractional progress via
    ``progress_cb`` if given."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_suffix(dest.suffix + ".partial")
    if partial.exists():
        partial.unlink()

    total_bytes = int(src_gz.stat().st_size) if src_gz.exists() else 0
    read = 0
    chunk = 1 << 20  # 1 MiB
    with gzip.open(src_gz, "rb") as fin, open(partial, "wb") as fout:
        while True:
            buf = fin.read(chunk)
            if not buf:
                break
            fout.write(buf)
            read += len(buf)
            if progress_cb and total_bytes > 0:
                # Report compressed-bytes progress as a coarse indicator;
                # uncompressed total is unknown without a second pass.
                progress_cb(
                    min(1.0, read / max(1, total_bytes * 3)), f"unzipped {read / 1e9:.2f} GB"
                )
            elif progress_cb:
                progress_cb(0.0, f"unzipped {read / 1e9:.2f} GB")
    partial.rename(dest)


def ensure_genome_fasta(
    genome: str,
    *,
    progress_cb: Callable[[float, str], None] | None = None,
) -> Path:
    """Return the path to the FASTA for ``genome``, downloading + indexing
    on first use.

    FASTAs are big (~1 GB compressed, ~3 GB unzipped) so this can take
    several minutes on a slow link. The download is staged through a
    ``.fa.gz.partial`` file, gunzipped to ``.fa.partial``, then renamed
    atomically. ``samtools faidx`` runs at the end to produce the
    ``.fa.fai`` index modkit needs.

    Already-on-disk FASTA (with non-empty size) short-circuits to
    just ensuring the ``.fai`` is present.

    Raises ``RuntimeError`` on failure (network, gzip integrity,
    samtools missing). Cleans up partial files on any failure so the
    next attempt starts fresh.
    """
    assets = assets_for(genome)
    fa = assets.fasta
    fa_fai = Path(str(fa) + ".fai")

    # Fast path: FASTA + index already present.
    if fa.exists() and fa.stat().st_size > 0 and fa_fai.exists():
        return fa

    REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    fa_gz = fa.with_suffix(".fa.gz")
    partials = [
        fa_gz.with_suffix(fa_gz.suffix + ".partial"),
        fa.with_suffix(fa.suffix + ".partial"),
    ]

    def _cleanup_on_failure() -> None:
        for p in partials:
            try:
                if p.exists():
                    p.unlink()
            except OSError:
                logger.warning("Could not remove %s during cleanup", p)
        # Don't delete a healthy unzipped FASTA — only the gz we own.
        try:
            if fa_gz.exists():
                fa_gz.unlink()
        except OSError:
            logger.warning("Could not remove %s during cleanup", fa_gz)

    try:
        # 1. Download .fa.gz (skip if a healthy unzipped FASTA already exists).
        if not (fa.exists() and fa.stat().st_size > 0):
            if progress_cb:
                progress_cb(0.0, f"downloading {genome}.fa.gz (~1 GB)…")
            _download_to(assets.fasta_url, fa_gz, progress_cb=progress_cb)

            # 2. gunzip to .fa.
            if progress_cb:
                progress_cb(0.0, f"unzipping {genome}.fa (~3 GB)…")
            _gunzip_to(fa_gz, fa, progress_cb=progress_cb)

            # Drop the .gz once unzipped — saves ~1 GB.
            try:
                fa_gz.unlink()
            except OSError:
                logger.warning("Could not remove %s after unzip", fa_gz)

        # 3. samtools faidx for modkit.
        if not fa_fai.exists():
            if shutil.which("samtools") is None:
                raise RuntimeError(
                    "samtools not on PATH — cannot index FASTA. "
                    "Install via `bash install.sh` or `conda install -c bioconda samtools`."
                )
            if progress_cb:
                progress_cb(0.95, f"samtools faidx {genome}.fa…")
            subprocess.run(  # noqa: S603 (samtools is a vetted bioconda binary)
                ["samtools", "faidx", str(fa)],
                check=True,
                capture_output=True,
                text=True,
            )
            if progress_cb:
                progress_cb(1.0, "done")
    except Exception as exc:  # noqa: BLE001 — surface as a single retryable error
        _cleanup_on_failure()
        raise RuntimeError(
            f"Failed to provision {genome} FASTA from {assets.fasta_url}: {exc}. "
            "Check internet access and disk space, or run "
            f"`bash scripts/setup_references.sh {genome}` manually."
        ) from exc
    return fa


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
