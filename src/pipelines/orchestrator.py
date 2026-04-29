"""High-level orchestration: a sample (or pair of samples) through
the full pipeline. Used by the Streamlit tab AND the CLI so the two
share one code path.

A "sample" is described by ``SampleSpec`` (see ``src/state.py``):

- ``input_kind = "bam"``         single BAM
- ``input_kind = "tsv"``         pre-computed modkit TSV (skip extract)
- ``input_kind = "bam_folder"``  folder of BAMs → merge → sort → index

The orchestrator:

1. Resolves the sample's BAM (merging if it's a folder).
2. Runs ``modkit extract`` if needed.
3. Runs the entropy pipeline (true_mc / bisulfite / ternary).
4. Returns a ``SampleResult`` with all output paths and the run summary.
"""

from __future__ import annotations

import logging
import re
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from src.constants import ENTROPY_MODE_TERNARY
from src.pipelines.bam_utils import find_bams, merge_sort_index_bams
from src.pipelines.modkit_runner import run_modkit_extract_minimal
from src.pipelines.ternary_entropy import (
    required_coverage_for_k,
    run_whole_genome_ternary,
)
from src.pipelines.whole_genome_duckdb_pipeline import run_whole_genome_duckdb_only

logger = logging.getLogger(__name__)


@dataclass
class SampleResult:
    """Output paths produced for one sample."""

    label: str
    out_prefix: Path
    me_bedgraph: Path
    mml_bedgraph: Path
    coverage_bedgraph: Path
    mhml_bedgraph: Path | None = None  # ternary only
    bam_used: Path | None = None  # final BAM that fed modkit
    tsv_path: Path | None = None
    extras: dict = field(default_factory=dict)


def _resolve_bam(
    spec, out_dir: Path, threads: int,
    progress_cb: Callable[[str], None] | None,
    pct_cb: Callable[[float, str], None] | None,
) -> Path | None:
    """Return the BAM to feed modkit. None if the spec points at a TSV."""
    if spec.input_kind == "tsv":
        return None
    if spec.input_kind == "bam":
        if not spec.bam_path:
            raise ValueError(f"sample '{spec.label}': bam_path not set")
        return Path(spec.bam_path).expanduser().resolve()
    if spec.input_kind == "bam_folder":
        if not spec.bam_folder:
            raise ValueError(f"sample '{spec.label}': bam_folder not set")
        bams = find_bams(spec.bam_folder)
        if not bams:
            raise ValueError(
                f"sample '{spec.label}': no .bam files in {spec.bam_folder}"
            )
        if progress_cb:
            progress_cb(f"[bam] found {len(bams)} BAMs in {spec.bam_folder}")
        merged = out_dir / f"{spec.label}_merged.sorted.bam"
        return merge_sort_index_bams(
            bams, merged, threads=threads,
            progress_cb=progress_cb, pct_cb=pct_cb,
        )
    raise ValueError(f"unknown input_kind: {spec.input_kind!r}")


def _ensure_tsv(
    spec,
    out_dir: Path,
    fasta: Path,
    threads: int,
    progress_cb: Callable[[str], None] | None,
    pct_cb: Callable[[float, str], None] | None,
    status_cb: Callable[[str], None] | None = None,
) -> tuple[Path, Path | None]:
    """Return (tsv_path, bam_used_or_None). Runs modkit if needed.

    ``status_cb`` is called with ephemeral lines from modkit's progress
    bar (terminated by ``\\r``). When provided, this lets the UI render
    the live updating bar in place. When None, modkit's progress bar
    falls through to ``progress_cb`` so callers that don't distinguish
    still see it.
    """
    if spec.input_kind == "tsv":
        if not spec.tsv_path:
            raise ValueError(f"sample '{spec.label}': tsv_path not set")
        return Path(spec.tsv_path).expanduser().resolve(), None

    bam = _resolve_bam(spec, out_dir, threads, progress_cb, pct_cb)
    assert bam is not None
    tsv = out_dir / f"{spec.label}_modkit_extract.tsv"
    log_filepath = str(out_dir / f"{spec.label}_modkit.log")

    if pct_cb:
        pct_cb(0.0, "modkit extract")
    if progress_cb:
        progress_cb(f"[modkit] {spec.label}: BAM → TSV")

    run_modkit_extract_minimal(
        bam_path=str(bam),
        out_tsv_path=str(tsv),
        reference_fasta=str(fasta),
        threads=threads,
        log_filepath=log_filepath,
        stream_cb=progress_cb,
        status_cb=status_cb,
    )
    if pct_cb:
        pct_cb(1.0, "modkit extract done")
    return tsv, bam


_LABEL_SAFE = re.compile(r"[^A-Za-z0-9._-]")


def _safe_label(label: str) -> str:
    """Strip path-separator and shell-special characters from a sample
    label before using it as a path component. Caps at 64 chars."""
    cleaned = _LABEL_SAFE.sub("_", label).strip("_")
    return (cleaned or "sample")[:64]


def run_one_sample(
    spec,
    *,
    out_dir: Path,
    fasta: Path,
    entropy_mode: str,
    cpgs_per_bin: int,
    min_coverage: int,
    methyl_threshold: float,
    threads: int,
    chroms: str = "",
    force_ingest: bool = False,
    progress_cb: Callable[[str], None] | None = None,
    pct_cb: Callable[[float, str], None] | None = None,
    status_cb: Callable[[str], None] | None = None,
) -> SampleResult:
    """Run the full BAM(/folder)/TSV → entropy bedgraphs pipeline for
    one ``SampleSpec``."""

    # Defensive: never let a label like "../bad" escape out_dir.
    spec.label = _safe_label(spec.label)

    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = Path(fasta).expanduser().resolve()

    if entropy_mode == ENTROPY_MODE_TERNARY:
        rec = required_coverage_for_k(int(cpgs_per_bin))
        if int(min_coverage) < rec and progress_cb:
            progress_cb(
                f"[warn] ternary at k={cpgs_per_bin} ⇒ 3^k={rec}; "
                f"min_coverage={min_coverage} is below 3^k — results will "
                "be noisy. Recommended ≥ "
                f"{rec}× coverage."
            )

    tsv_path, bam_used = _ensure_tsv(
        spec, out_dir, fasta, threads, progress_cb, pct_cb,
        status_cb=status_cb,
    )

    out_prefix = out_dir / f"{spec.label}_{entropy_mode}"
    if entropy_mode == ENTROPY_MODE_TERNARY:
        run_whole_genome_ternary(
            tsv_path=str(tsv_path),
            fasta_path=str(fasta),
            out_prefix=str(out_prefix),
            work_dir=str(out_dir / f"{spec.label}_wg_work_ternary"),
            threads=int(threads),
            cpgs_per_bin=int(cpgs_per_bin),
            methyl_thresh=float(methyl_threshold),
            min_coverage=int(min_coverage),
            chroms=chroms,
            force_ingest=bool(force_ingest),
            progress_cb=progress_cb,
            pct_cb=pct_cb,
        )
    else:
        run_whole_genome_duckdb_only(
            tsv_path=str(tsv_path),
            fasta_path=str(fasta),
            out_prefix=str(out_prefix),
            work_dir=str(out_dir / f"{spec.label}_wg_work_{entropy_mode}"),
            threads=int(threads),
            cpgs_per_bin=int(cpgs_per_bin),
            methyl_thresh=float(methyl_threshold),
            min_coverage=int(min_coverage),
            chroms=chroms,
            force_ingest=bool(force_ingest),
            entropy_mode=entropy_mode,
            progress_cb=progress_cb,
            pct_cb=pct_cb,
        )

    res = SampleResult(
        label=spec.label,
        out_prefix=out_prefix,
        me_bedgraph=Path(f"{out_prefix}.me.bedgraph"),
        mml_bedgraph=Path(f"{out_prefix}.mml.bedgraph"),
        coverage_bedgraph=Path(f"{out_prefix}.coverage.bedgraph"),
        bam_used=bam_used,
        tsv_path=tsv_path,
    )
    if entropy_mode == ENTROPY_MODE_TERNARY:
        res.mhml_bedgraph = Path(f"{out_prefix}.mhml.bedgraph")
    return res


def run_pipeline(
    target_spec,
    control_spec=None,
    *,
    out_dir: Path,
    fasta: Path,
    entropy_mode: str = "true_mc",
    cpgs_per_bin: int = 4,
    min_coverage: int = 16,
    methyl_threshold: float = 0.5,
    threads: int = 8,
    chroms: str = "",
    force_ingest: bool = False,
    progress_cb: Callable[[str], None] | None = None,
    pct_cb: Callable[[float, str], None] | None = None,
    status_cb: Callable[[str], None] | None = None,
) -> tuple[SampleResult, SampleResult | None]:
    """Run one or two samples through the pipeline.

    If ``control_spec`` is None, only ``target_spec`` runs and the
    second tuple element is None. Returns ``(target_result, control_result)``.
    """
    if control_spec is None:
        target = run_one_sample(
            target_spec,
            out_dir=out_dir, fasta=fasta,
            entropy_mode=entropy_mode,
            cpgs_per_bin=cpgs_per_bin,
            min_coverage=min_coverage,
            methyl_threshold=methyl_threshold,
            threads=threads, chroms=chroms,
            force_ingest=force_ingest,
            progress_cb=progress_cb, pct_cb=pct_cb,
            status_cb=status_cb,
        )
        return target, None

    # Pair mode: split progress 50/50 between control and target.
    def half_pct(offset: float):
        if pct_cb is None:
            return None

        def cb(f: float, msg: str = "") -> None:
            pct_cb(offset + 0.5 * f, msg)

        return cb

    if progress_cb:
        progress_cb("[pair] running CONTROL sample")
    control = run_one_sample(
        control_spec,
        out_dir=out_dir, fasta=fasta,
        entropy_mode=entropy_mode,
        cpgs_per_bin=cpgs_per_bin,
        min_coverage=min_coverage,
        methyl_threshold=methyl_threshold,
        threads=threads, chroms=chroms,
        force_ingest=force_ingest,
        progress_cb=progress_cb, pct_cb=half_pct(0.0),
        status_cb=status_cb,
    )
    if progress_cb:
        progress_cb("[pair] running TARGET sample")
    target = run_one_sample(
        target_spec,
        out_dir=out_dir, fasta=fasta,
        entropy_mode=entropy_mode,
        cpgs_per_bin=cpgs_per_bin,
        min_coverage=min_coverage,
        methyl_threshold=methyl_threshold,
        threads=threads, chroms=chroms,
        force_ingest=force_ingest,
        progress_cb=progress_cb, pct_cb=half_pct(0.5),
        status_cb=status_cb,
    )
    return target, control
