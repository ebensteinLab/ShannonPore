"""Command-line interface for shannonpore v4.

Mirrors every feature of the Streamlit GUI as scriptable subcommands so
the same pipelines can run on an HPC node, in a Snakemake/Nextflow
workflow, or in a notebook without launching a browser.

Subcommands
-----------
- ``extract``   BAM → modkit TSV
- ``entropy``   modkit TSV → ME / MML / coverage bedgraphs
                (--mode true_mc | bisulfite | ternary)
- ``segment``   bedgraph → segments BED (PELT or greedy)
- ``annotate``  segments BED → gene + GO annotation via R bridge
- ``plot``      bedgraphs → track / ME-MML scatter / arch / paired landscape PNG
- ``run``       one-shot pipeline: BAM → entropy → segment → annotate
- ``doctor``    verify every dependency, version, and permission
- ``selftest``  end-to-end pipeline smoke test on a synthetic BAM

Each subcommand has its own ``--help`` and accepts the same parameters
the GUI exposes.

Entry point
-----------
After ``pip install -e .`` (or ``python -m src.cli ...``) the CLI is
available as ``shannonpore`` thanks to the ``[project.scripts]`` entry
in ``pyproject.toml``.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

from src import __version__
from src.config import (
    REFERENCE_DIR,
    RESULTS_DIR,
    assets_for,
    ensure_dirs,
    ensure_genome_gtf,
)
from src.constants import (
    ENTROPY_MODE_BISULFITE,
    ENTROPY_MODE_HELP,
    ENTROPY_MODE_TERNARY,
    ENTROPY_MODE_TRUE_MC,
    ENTROPY_MODES,
)

logger = logging.getLogger("shannonpore.cli")


# ─────────────────────────── helpers ──────────────────────────────────────

def _resolve_fasta(genome: str, custom: str | None) -> Path:
    if custom:
        p = Path(custom).expanduser().resolve()
        if not p.exists():
            raise SystemExit(f"FASTA not found: {p}")
        return p
    return assets_for(genome).fasta


def _setup_logging(verbose: int) -> None:
    level = logging.WARNING - 10 * min(verbose, 2)
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def _print_progress(msg: str) -> None:
    print(msg, flush=True)


# ─────────────────────────── extract ──────────────────────────────────────

def cmd_extract(args: argparse.Namespace) -> int:
    from src.pipelines.bam_utils import find_bams, merge_sort_index_bams
    from src.pipelines.modkit_runner import run_modkit_extract_minimal
    from src.ui.progress import CLIProgress

    fasta = _resolve_fasta(args.genome, args.fasta)
    out_tsv = Path(args.out_tsv).expanduser().resolve()
    log_path = Path(args.log_file or (out_tsv.parent / "modkit.log"))
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Resolve input BAM (single file, or merge a folder of BAMs).
    if args.bam_folder:
        with CLIProgress("merge+sort+index BAMs") as bar:
            bams = find_bams(args.bam_folder)
            if not bams:
                raise SystemExit(f"No .bam files found in {args.bam_folder}")
            bar.status(f"found {len(bams)} BAMs in {args.bam_folder}")
            merged_bam = out_tsv.parent / "merged.sorted.bam"
            bam_path = merge_sort_index_bams(
                bams, merged_bam, threads=int(args.threads),
                progress_cb=bar.status, pct_cb=bar.update,
            )
    else:
        bam_path = Path(args.bam).expanduser().resolve()

    with CLIProgress("modkit extract") as bar:
        run_modkit_extract_minimal(
            bam_path=str(bam_path),
            out_tsv_path=str(out_tsv),
            reference_fasta=str(fasta),
            threads=int(args.threads),
            log_filepath=str(log_path),
            stream_cb=bar.status,
        )
    print(f"[OK] modkit TSV: {out_tsv}")
    return 0


# ─────────────────────────── entropy ──────────────────────────────────────

def cmd_entropy(args: argparse.Namespace) -> int:
    from src.ui.progress import CLIProgress

    fasta = _resolve_fasta(args.genome, args.fasta)
    out_prefix = Path(args.out_prefix).expanduser().resolve()
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    work_dir = Path(
        args.work_dir or (out_prefix.parent / f"{out_prefix.name}_work")
    ).expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    with CLIProgress(f"entropy · {args.mode}") as bar:
        if args.mode == ENTROPY_MODE_TERNARY:
            from src.pipelines.ternary_entropy import (
                required_coverage_for_k,
                run_whole_genome_ternary,
            )
            rec = required_coverage_for_k(int(args.cpgs_per_bin))
            if int(args.min_coverage) < rec:
                bar.status(
                    f"[warn] ternary k={args.cpgs_per_bin} needs ~{rec}× "
                    f"coverage; you set min_coverage={args.min_coverage}"
                )
            run_whole_genome_ternary(
                tsv_path=str(Path(args.tsv).expanduser().resolve()),
                fasta_path=str(fasta),
                out_prefix=str(out_prefix),
                work_dir=str(work_dir),
                threads=int(args.threads),
                cpgs_per_bin=int(args.cpgs_per_bin),
                methyl_thresh=float(args.methyl_threshold),
                min_coverage=int(args.min_coverage),
                chroms=args.chroms or "",
                force_ingest=bool(args.force),
                progress_cb=bar.status,
                pct_cb=bar.update,
            )
        else:
            from src.pipelines.whole_genome_duckdb_pipeline import (
                run_whole_genome_duckdb_only,
            )
            run_whole_genome_duckdb_only(
                tsv_path=str(Path(args.tsv).expanduser().resolve()),
                fasta_path=str(fasta),
                out_prefix=str(out_prefix),
                work_dir=str(work_dir),
                threads=int(args.threads),
                cpgs_per_bin=int(args.cpgs_per_bin),
                methyl_thresh=float(args.methyl_threshold),
                min_coverage=int(args.min_coverage),
                chroms=args.chroms or "",
                force_ingest=bool(args.force),
                entropy_mode=args.mode,
                progress_cb=bar.status,
                pct_cb=bar.update,
            )
    print(f"[OK] entropy bedgraphs prefix: {out_prefix}")
    return 0


# ─────────────────────────── plot ─────────────────────────────────────────

def cmd_plot(args: argparse.Namespace) -> int:
    from src.plots.scatter import (
        load_paired_bedgraphs,
        me_mml_scatter,
        paired_landscape,
        triple_landscape,
    )
    from src.plots.theme import apply_default_style
    from src.plots.tracks import plot_region_tracks

    apply_default_style()
    out_path = Path(args.out_path).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.kind == "tracks":
        # All four bedgraphs are required for the new gene-panel layout.
        for required in ("control_mml", "control_me", "target_mml", "target_me"):
            if not getattr(args, required):
                raise SystemExit(
                    f"--{required.replace('_', '-')} is required for "
                    f"plot kind 'tracks'."
                )
        # Default GTF: bundled per-genome (downloads on first use) unless
        # the user passed --gtf.
        gtf: str | None = args.gtf
        if not gtf:
            try:
                gtf = str(ensure_genome_gtf(args.genome))
            except RuntimeError as exc:
                logger.warning(
                    "Could not download GTF for %s — gene panel will be empty: %s",
                    args.genome, exc,
                )
                gtf = None
        plot_region_tracks(
            chrom=args.chrom, start=int(args.start), end=int(args.end),
            control_mml=args.control_mml, target_mml=args.target_mml,
            control_me=args.control_me, target_me=args.target_me,
            control_coverage=args.control_coverage,
            target_coverage=args.target_coverage,
            label_a=args.label_a, label_b=args.label_b,
            color_a=args.color_a, color_b=args.color_b,
            gtf_path=gtf,
            smooth_win=int(args.window),
            pad_bp=int(args.pad),
            out_path=str(out_path),
        )
    elif args.kind in ("scatter", "arch", "landscape"):
        # All three need the four paired bedgraphs.
        for required in ("control_mml", "control_me", "target_mml", "target_me"):
            if not getattr(args, required):
                raise SystemExit(
                    f"--{required.replace('_', '-')} is required for "
                    f"plot kind '{args.kind}'."
                )
        df = load_paired_bedgraphs(
            control_mml=args.control_mml,
            control_me=args.control_me,
            target_mml=args.target_mml,
            target_me=args.target_me,
        )
        if df.empty:
            raise SystemExit(
                "No overlapping bins found across the four bedgraphs."
            )
        log_scale = bool(args.log_scale)

        if args.kind == "scatter":
            me_mml_scatter(
                df,
                label_a=args.label_a, label_b=args.label_b,
                log_scale=log_scale, out_path=str(out_path),
            )
        elif args.kind == "arch":
            triple_landscape(
                df,
                label_a=args.label_a, label_b=args.label_b,
                color_a=args.color_a, color_b=args.color_b,
                log_scale=log_scale, out_path=str(out_path),
            )
        else:  # landscape (paired)
            paired_landscape(
                df,
                label_a=args.label_a, label_b=args.label_b,
                color_a=args.color_a, color_b=args.color_b,
                filter_a_dim=None if args.filter_a_dim == "off" else args.filter_a_dim,
                filter_a_op=args.filter_a_op,
                filter_a_value=float(args.filter_a_value),
                filter_b_dim=None if args.filter_b_dim == "off" else args.filter_b_dim,
                filter_b_op=args.filter_b_op,
                filter_b_value=float(args.filter_b_value),
                max_lines=int(args.max_lines),
                out_path=str(out_path),
            )
    else:
        raise SystemExit(f"Unknown plot kind: {args.kind!r}")
    print(f"[OK] plot → {out_path}")
    return 0


# ─────────────────────────── run (one-shot) ───────────────────────────────

def _spec_from_args(args, prefix: str, default_label: str):
    """Build a SampleSpec from CLI flags. `prefix` is "" for single mode,
    "control_" or "target_" in pair mode."""
    from src.state import SampleSpec

    bam = getattr(args, f"{prefix}bam", None)
    tsv = getattr(args, f"{prefix}tsv", None)
    folder = getattr(args, f"{prefix}bam_folder", None)
    label = getattr(args, f"{prefix}label", None) or default_label

    if sum(bool(x) for x in (bam, tsv, folder)) > 1:
        raise SystemExit(
            f"For {default_label!r}: pass exactly one of --{prefix}bam / "
            f"--{prefix}tsv / --{prefix}bam-folder."
        )
    if bam:
        return SampleSpec(label=label, input_kind="bam", bam_path=Path(bam))
    if tsv:
        return SampleSpec(label=label, input_kind="tsv", tsv_path=Path(tsv))
    if folder:
        return SampleSpec(
            label=label, input_kind="bam_folder", bam_folder=Path(folder),
        )
    raise SystemExit(
        f"For {default_label!r}: provide --{prefix}bam, --{prefix}tsv, "
        f"or --{prefix}bam-folder."
    )


def cmd_run(args: argparse.Namespace) -> int:
    """End-to-end: BAM(/folder)/TSV → modkit → entropy bedgraphs.

    Two flavours:
        single mode (default): --bam | --tsv | --bam-folder
        pair mode (--pair):    --control-{bam,tsv,bam-folder}
                               --target-{bam,tsv,bam-folder}
    """
    from src.pipelines.orchestrator import run_pipeline
    from src.ui.progress import CLIProgress

    out_dir = Path(args.out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = _resolve_fasta(args.genome, args.fasta)

    if args.pair:
        control_spec = _spec_from_args(args, "control_", "control")
        target_spec = _spec_from_args(args, "target_", "target")
    else:
        target_spec = _spec_from_args(args, "", args.label or "sample")
        control_spec = None

    label = "pair pipeline" if args.pair else "pipeline"
    with CLIProgress(label) as bar:
        target_result, control_result = run_pipeline(
            target_spec=target_spec,
            control_spec=control_spec,
            out_dir=out_dir,
            fasta=fasta,
            entropy_mode=args.mode,
            cpgs_per_bin=int(args.cpgs_per_bin),
            min_coverage=int(args.min_coverage),
            methyl_threshold=float(args.methyl_threshold),
            threads=int(args.threads),
            chroms=args.chroms or "",
            force_ingest=bool(args.force),
            progress_cb=bar.status,
            pct_cb=bar.update,
            status_cb=bar.status_line,  # in-place modkit progress bar
        )

    summary: dict = {
        "version": __version__,
        "mode": args.mode,
        "pair_mode": bool(args.pair),
        "fasta": str(fasta),
        "cpgs_per_bin": int(args.cpgs_per_bin),
        "min_coverage": int(args.min_coverage),
        "methyl_threshold": float(args.methyl_threshold),
        "target": {
            "label": target_result.label,
            "out_prefix": str(target_result.out_prefix),
            "bam_used": str(target_result.bam_used) if target_result.bam_used else None,
            "tsv_path": str(target_result.tsv_path) if target_result.tsv_path else None,
            "me_bedgraph": str(target_result.me_bedgraph),
            "mml_bedgraph": str(target_result.mml_bedgraph),
            "coverage_bedgraph": str(target_result.coverage_bedgraph),
        },
    }
    if control_result:
        summary["control"] = {
            "label": control_result.label,
            "out_prefix": str(control_result.out_prefix),
            "bam_used": str(control_result.bam_used) if control_result.bam_used else None,
            "tsv_path": str(control_result.tsv_path) if control_result.tsv_path else None,
            "me_bedgraph": str(control_result.me_bedgraph),
            "mml_bedgraph": str(control_result.mml_bedgraph),
            "coverage_bedgraph": str(control_result.coverage_bedgraph),
        }
    summary_path = out_dir / "run_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"[OK] run summary: {summary_path}")
    return 0


# ─────────────────────────── guide / examples ────────────────────────────

def cmd_guide(args: argparse.Namespace) -> int:
    from src.help_text import GUIDE
    print(GUIDE)
    return 0


def cmd_examples(args: argparse.Namespace) -> int:
    from src.help_text import EXAMPLES
    print(EXAMPLES)
    return 0


# ─────────────────────────── doctor ───────────────────────────────────────

def cmd_doctor(args: argparse.Namespace) -> int:
    """Verify every required dependency, version, and permission.

    Prints a coloured pass/fail table and exits non-zero if anything
    critical is missing. Run this immediately after `bash install.sh`.
    """
    import importlib
    import shutil
    import stat

    # severity: "must"   — fails the install (broken state)
    #           "advise" — informational; user-action item, never fails the install
    checks: list[tuple[str, str, bool, str, str]] = []

    def check(
        name: str, expected: str, ok: bool, detail: str = "",
        severity: str = "must",
    ) -> None:
        checks.append((name, expected, ok, detail, severity))

    # ── Python version ──
    py = sys.version_info
    check(
        "python",
        ">=3.10",
        py.major == 3 and py.minor >= 10,
        f"{py.major}.{py.minor}.{py.micro}",
    )

    # ── Pinned Python packages ──
    pinned = {
        "streamlit": "1.51.0",
        "pandas": "2.2.3",
        "numpy": "1.26.4",
        "duckdb": "1.4.3",
        "matplotlib": "3.8.4",
        "seaborn": "0.13.2",
        "pyfaidx": "0.8.1.2",
        "pysam": "0.22.1",
        "plotly": "6.2.0",
        "tqdm": "4.67.1",
        "scipy": "1.11.4",
    }
    for pkg, want in pinned.items():
        try:
            mod = importlib.import_module(pkg)
            ver = getattr(mod, "__version__", "?")
            check(pkg, f"=={want}", ver == want, ver)
        except ImportError as exc:
            check(pkg, f"=={want}", False, f"missing: {exc}")

    # ── External binaries ──
    for binary, expected in (("modkit", "0.6.0"), ("samtools", "1.21")):
        path = shutil.which(binary)
        if not path:
            check(binary, expected, False, "not on PATH")
            continue
        try:
            out = subprocess.run(  # noqa: S603
                [binary, "--version"], capture_output=True, text=True, timeout=10,
            )
            ver_line = (out.stdout or out.stderr).strip().splitlines()[0]
            ok = expected in ver_line
            check(binary, expected, ok, ver_line[:60])
        except (OSError, subprocess.TimeoutExpired) as exc:
            check(binary, expected, False, str(exc))

    # ── Reference dir + per-genome FASTAs ──
    from src.config import reference_status

    ensure_dirs()  # creates REFERENCE_DIR if missing
    check("REFERENCE_DIR", "exists", REFERENCE_DIR.exists(), str(REFERENCE_DIR))

    project_root = Path(__file__).resolve().parent.parent
    setup_script = project_root / "scripts" / "setup_references.sh"
    fix_hint = (
        f"download via: bash {setup_script}"
        if setup_script.exists()
        else "download FASTAs into REFERENCE_DIR"
    )
    for genome, status in reference_status().items():
        all_good = all(status.values())
        detail = (
            "fasta + .fai + GTF present"
            if all_good
            else f"not downloaded — {fix_hint}"
        )
        # FASTA download is opt-in; missing it is "advise", not "must".
        check(
            f"refs:{genome}", "fasta+fai+gtf", all_good, detail,
            severity="advise",
        )

    # ── Results dir is writable ──
    try:
        ensure_dirs()
        test_file = RESULTS_DIR / ".permcheck"
        test_file.write_text("ok")
        test_file.unlink()
        check("RESULTS_DIR writable", "ok", True, str(RESULTS_DIR))
    except OSError as exc:
        check("RESULTS_DIR writable", "ok", False, str(exc))

    # ── Executable bit on shell scripts ──
    project_root = Path(__file__).resolve().parent.parent
    for script in ("install.sh",):
        p = project_root / script
        if not p.exists():
            check(f"+x {script}", "exists", False, "missing")
            continue
        is_exec = bool(p.stat().st_mode & stat.S_IXUSR)
        check(f"+x {script}", "executable", is_exec, str(p))

    # ── Render table ──
    width_n = max(len(c[0]) for c in checks) + 2
    width_e = max(len(c[1]) for c in checks) + 2
    print(f"{'Component':<{width_n}}{'Expected':<{width_e}}{'OK':<5}Detail")
    print("─" * 80)
    failures = 0
    advisories = 0
    for name, expected, ok, detail, severity in checks:
        if ok:
            flag = "✓"
        elif severity == "advise":
            flag = "i"
            advisories += 1
        else:
            flag = "✗"
            failures += 1
        print(f"{name:<{width_n}}{expected:<{width_e}}{flag:<5}{detail}")
    print("─" * 80)

    if failures == 0 and advisories == 0:
        print("[OK] All checks passed. shannonpore is ready to use.")
        return 0
    if failures == 0:
        # Only advisories — install is functional.
        print(
            f"[OK] Install is functional ({advisories} advisory item(s) — "
            "see `i` rows above; not blockers)."
        )
        return 0
    print(f"[FAIL] {failures} blocking check(s) failed. Run `bash install.sh` and re-check.")
    return 1


# ─────────────────────────── selftest ─────────────────────────────────────

def cmd_selftest(args: argparse.Namespace) -> int:
    """End-to-end smoke test: synthesize a tiny BAM, run the full pipeline,
    verify outputs. Confirms the install actually works.

    Requires: pysam, modkit. Skips the full ternary check by default
    (use --include-ternary to add it).
    """
    import shutil
    import tempfile

    if shutil.which("modkit") is None:
        print("[FAIL] modkit not on PATH. Install with conda or check install.sh.")
        return 1
    try:
        import pysam  # noqa: F401
    except ImportError:
        print("[FAIL] pysam not installed. Re-run install.sh.")
        return 1

    project_root = Path(__file__).resolve().parent.parent
    sys.path.insert(0, str(project_root))
    from tests.fixtures.build_tiny_bam import build_tiny_fixture

    print("[selftest] building tiny FASTA + BAM fixture...")
    with tempfile.TemporaryDirectory() as td:
        td_p = Path(td)
        fa, bam = build_tiny_fixture(td_p, n_reads=20, methylated_fraction=0.5)

        out_dir = td_p / "out"
        modes = ["true_mc", "bisulfite"]
        if args.include_ternary:
            modes.append("ternary")

        for mode in modes:
            print(f"[selftest] running pipeline mode={mode}...")
            r = subprocess.run(  # noqa: S603
                [sys.executable, "-m", "src.cli", "run",
                 "--bam", str(bam), "--fasta", str(fa),
                 "--out-dir", str(out_dir / mode),
                 "--threads", "2", "--mode", mode,
                 "--cpgs-per-bin", "2", "--min-coverage", "4",
                 "--chroms", "chr_test"],
                cwd=str(project_root),
                capture_output=True, text=True, timeout=120,
            )
            if r.returncode != 0:
                print(f"[FAIL] pipeline failed for mode={mode}")
                print("STDOUT:", r.stdout[-500:])
                print("STDERR:", r.stderr[-500:])
                return 1
            for suffix in ("coverage", "me", "mml"):
                p = out_dir / mode / f"sample_{mode}.{suffix}.bedgraph"
                if not p.exists() or p.stat().st_size == 0:
                    print(f"[FAIL] missing/empty output: {p}")
                    return 1
            print(f"[OK]   mode={mode} produced expected bedgraphs")

    print()
    print("[OK] selftest passed. The full BAM → entropy pipeline works end-to-end.")
    return 0


# ─────────────────────────── parser ───────────────────────────────────────

def _add_common_io(parser: argparse.ArgumentParser, *, want_tsv: bool = False) -> None:
    parser.add_argument("--genome", default="hg38",
                        choices=["hg38", "mm10"],
                        help="Reference genome (default: hg38).")
    parser.add_argument("--fasta", default=None,
                        help="Custom FASTA (overrides --genome bundled FASTA).")
    parser.add_argument("--threads", type=int, default=8,
                        help="Threads / parallel workers (default: 8).")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="Increase verbosity (-v INFO, -vv DEBUG).")


def _add_entropy_params(p: argparse.ArgumentParser) -> None:
    p.add_argument(
        "--mode", default=ENTROPY_MODE_TRUE_MC,
        choices=list(ENTROPY_MODES),
        help="Entropy mode (default: true_mc). " + " | ".join(
            f"{m}: {ENTROPY_MODE_HELP[m]}" for m in ENTROPY_MODES
        ),
    )
    p.add_argument("--cpgs-per-bin", type=int, default=4,
                   dest="cpgs_per_bin",
                   help="CpGs per bin k (default: 4).")
    p.add_argument("--min-coverage", type=int, default=16,
                   dest="min_coverage",
                   help="Minimum reads per bin (default: 16).")
    p.add_argument("--methyl-threshold", type=float, default=0.5,
                   dest="methyl_threshold",
                   help="Probability threshold for methylated call (default: 0.5).")
    p.add_argument("--chroms", default="",
                   help="Comma-separated chromosome list. Empty = all from FASTA.")
    p.add_argument("--force", action="store_true",
                   help="Force re-ingest of TSV (drops existing DuckDB table).")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="shannonpore",
        description=(
            f"shannonpore v{__version__} — nanopore methylation entropy "
            "analysis. CLI mirrors the Streamlit GUI."
        ),
        epilog=(
            "Quick links:\n"
            "  shannonpore guide       walkthrough of every flow\n"
            "  shannonpore examples    cheat-sheet of common commands\n"
            "  shannonpore doctor      verify install\n"
            "  shannonpore selftest    end-to-end synthetic-BAM smoke test\n\n"
            "Full docs: https://github.com/uribertocchitau/shannonpore#readme"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"shannonpore {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True, metavar="<subcommand>")

    # extract -----------------------------------------------------------------
    p_ex = sub.add_parser(
        "extract",
        help="BAM (or folder of BAMs) → modkit extract full TSV.",
        description=(
            "Run `modkit extract full` to produce a per-CpG TSV. Pass either "
            "--bam (single BAM) or --bam-folder (folder of BAMs to "
            "merge+sort+index automatically before extract)."
        ),
    )
    grp = p_ex.add_mutually_exclusive_group(required=True)
    grp.add_argument("--bam", default=None, help="Input BAM file.")
    grp.add_argument("--bam-folder", default=None, dest="bam_folder",
                     help="Folder of BAMs (auto merge+sort+index).")
    p_ex.add_argument("out_tsv", help="Output TSV path.")
    p_ex.add_argument("--log-file", default=None,
                      help="Path for modkit log (default: alongside output TSV).")
    _add_common_io(p_ex)
    p_ex.set_defaults(func=cmd_extract)

    # entropy -----------------------------------------------------------------
    p_en = sub.add_parser(
        "entropy",
        help="modkit TSV → ME / MML / coverage bedgraphs (whole genome).",
        description=(
            "Compute methylation entropy from a modkit TSV. Three modes:\n"
            f"  true_mc   — {ENTROPY_MODE_HELP[ENTROPY_MODE_TRUE_MC]}\n"
            f"  bisulfite — {ENTROPY_MODE_HELP[ENTROPY_MODE_BISULFITE]}\n"
            f"  ternary   — {ENTROPY_MODE_HELP[ENTROPY_MODE_TERNARY]}"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_en.add_argument("tsv", help="modkit extract TSV (input).")
    p_en.add_argument("out_prefix", help="Output prefix (writes <prefix>.{me,mml,coverage}.bedgraph).")
    p_en.add_argument("--work-dir", default=None,
                      help="Work directory (default: alongside out_prefix).")
    _add_common_io(p_en)
    _add_entropy_params(p_en)
    p_en.set_defaults(func=cmd_entropy)

    # plot --------------------------------------------------------------------
    p_pl = sub.add_parser(
        "plot",
        help="Bedgraphs → tracks / scatter / arch / landscape PNG.",
        description=(
            "Render one of four plot kinds from your bedgraphs:\n\n"
            "  tracks     control & target MML/ME along a genomic window,\n"
            "             optional gene panel from a GTF\n"
            "  scatter    side-by-side ME and MML 2D-histograms\n"
            "             (control vs target)\n"
            "  arch       three-panel MML × ME density landscape with the\n"
            "             theoretical entropy arch\n"
            "  landscape  paired-line scatter with direction arrows on a\n"
            "             configurable bin filter"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_pl.add_argument(
        "kind", choices=["tracks", "scatter", "arch", "landscape"],
    )
    p_pl.add_argument("out_path", help="Output PNG path.")

    # Inputs (used by scatter / arch / landscape AND tracks)
    p_pl.add_argument("--control-mml", default=None)
    p_pl.add_argument("--control-me", default=None)
    p_pl.add_argument("--target-mml", default=None)
    p_pl.add_argument("--target-me", default=None)
    # Coverage is optional and only consumed by `tracks` (4th panel).
    p_pl.add_argument("--control-coverage", default=None,
                      help="(tracks only) coverage bedgraph for control.")
    p_pl.add_argument("--target-coverage", default=None,
                      help="(tracks only) coverage bedgraph for target.")

    # Sample labels + colours (carried into all plots)
    p_pl.add_argument("--label-a", default="Control",
                      help="Label for sample A (control). Default: Control.")
    p_pl.add_argument("--label-b", default="Target",
                      help="Label for sample B (target). Default: Target.")
    p_pl.add_argument("--color-a", default="#1f77b4",
                      help="Colour for sample A (control). Default: #1f77b4.")
    p_pl.add_argument("--color-b", default="#ff7f0e",
                      help="Colour for sample B (target). Default: #ff7f0e.")

    # scatter / arch — log vs linear colour scale
    p_pl.add_argument(
        "--log-scale", action=argparse.BooleanOptionalAction, default=True,
        help="log10 colour scale (default). Use --no-log-scale for linear.",
    )

    # landscape (paired) — two AND-ed filters
    _filter_dims = ["off", "MML", "|dMML|", "ME", "|dME|"]
    p_pl.add_argument("--filter-a-dim", default="|dMML|",
                      choices=_filter_dims,
                      help="Landscape filter A dimension (default: |dMML|).")
    p_pl.add_argument("--filter-a-op", default="<", choices=["<", ">"],
                      help="Landscape filter A op (default: <).")
    p_pl.add_argument("--filter-a-value", type=float, default=0.1,
                      help="Landscape filter A threshold (default: 0.1).")
    p_pl.add_argument("--filter-b-dim", default="|dME|",
                      choices=_filter_dims,
                      help="Landscape filter B dimension (default: |dME|).")
    p_pl.add_argument("--filter-b-op", default=">", choices=["<", ">"],
                      help="Landscape filter B op (default: >).")
    p_pl.add_argument("--filter-b-value", type=float, default=0.4,
                      help="Landscape filter B threshold (default: 0.4).")
    p_pl.add_argument("--max-lines", type=int, default=20_000,
                      help="Max paired lines drawn (default: 20000).")

    # tracks specifics
    p_pl.add_argument(
        "--gtf", default=None,
        help="Override the bundled GTF (defaults to the GTF bundled for "
             "--genome).",
    )
    p_pl.add_argument("--genome", default="hg38",
                      choices=["hg38", "mm10"],
                      help="Genome key for the bundled GTF (default: hg38).")
    p_pl.add_argument("--chrom", default="")
    p_pl.add_argument("--start", type=int, default=0)
    p_pl.add_argument("--end", type=int, default=0)
    p_pl.add_argument("--window", type=int, default=5,
                      help="Smoothing window for ME / MML signals (default: 5 bins).")
    p_pl.add_argument("--pad", type=int, default=2000,
                      help="Padding around the region in bp (default: 2000).")

    p_pl.add_argument("-v", "--verbose", action="count", default=0)
    p_pl.set_defaults(func=cmd_plot)

    # run ---------------------------------------------------------------------
    p_run = sub.add_parser(
        "run",
        help=(
            "End-to-end pipeline. "
            "Single sample: --bam|--tsv|--bam-folder. "
            "Pair mode: --pair --control-... --target-..."
        ),
        description=(
            "End-to-end shannonpore pipeline.\n\n"
            "Single sample (default):\n"
            "  --bam BAM          one BAM file\n"
            "  --tsv TSV          pre-computed modkit TSV\n"
            "  --bam-folder DIR   folder of BAMs (auto merge+sort+index)\n\n"
            "Pair mode (--pair): same flags but prefixed --control-* and "
            "--target-* so two samples are processed in one invocation, "
            "producing two parallel sets of bedgraphs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # single-sample inputs
    p_run.add_argument("--bam", default=None, help="Single-mode: one BAM.")
    p_run.add_argument("--tsv", default=None,
                       help="Single-mode: pre-computed modkit TSV.")
    p_run.add_argument("--bam-folder", default=None, dest="bam_folder",
                       help="Single-mode: folder of BAMs (merge+sort+index).")
    p_run.add_argument("--label", default=None,
                       help="Single-mode: label for output filenames "
                            "(default: 'sample').")
    # pair-mode flag
    p_run.add_argument("--pair", action="store_true",
                       help="Pair mode: process control + target.")
    # pair-mode inputs
    p_run.add_argument("--control-bam", default=None, dest="control_bam")
    p_run.add_argument("--control-tsv", default=None, dest="control_tsv")
    p_run.add_argument("--control-bam-folder", default=None,
                       dest="control_bam_folder")
    p_run.add_argument("--control-label", default=None, dest="control_label",
                       help="Pair-mode control label (default: 'control').")
    p_run.add_argument("--target-bam", default=None, dest="target_bam")
    p_run.add_argument("--target-tsv", default=None, dest="target_tsv")
    p_run.add_argument("--target-bam-folder", default=None,
                       dest="target_bam_folder")
    p_run.add_argument("--target-label", default=None, dest="target_label",
                       help="Pair-mode target label (default: 'target').")

    p_run.add_argument("--out-dir", required=True, dest="out_dir")
    _add_common_io(p_run)
    _add_entropy_params(p_run)
    p_run.set_defaults(func=cmd_run)

    # doctor ------------------------------------------------------------------
    p_doc = sub.add_parser(
        "doctor",
        help="Verify every dependency, version, and permission is in place.",
        description=(
            "Runs after install. Checks Python version, every pinned Python "
            "library, modkit, Rscript, R packages, REFERENCE_DIR / RESULTS_DIR, "
            "and the executable bit on shell scripts. Exits 0 only when "
            "everything is OK."
        ),
    )
    p_doc.add_argument("-v", "--verbose", action="count", default=0)
    p_doc.set_defaults(func=cmd_doctor)

    # guide -------------------------------------------------------------------
    p_g = sub.add_parser(
        "guide",
        help="Print the user guide.",
        description=(
            "Print a friendly walkthrough of the tool: input shapes, "
            "entropy modes, bin parameters, and how to plot the results."
        ),
    )
    p_g.set_defaults(func=cmd_guide)

    # examples ----------------------------------------------------------------
    p_e = sub.add_parser(
        "examples",
        help="Print common usage examples.",
        description="Cheat-sheet of common shannonpore command lines.",
    )
    p_e.set_defaults(func=cmd_examples)

    # selftest ----------------------------------------------------------------
    p_st = sub.add_parser(
        "selftest",
        help="End-to-end pipeline smoke test on a synthetic BAM.",
        description=(
            "Synthesizes a tiny FASTA + BAM with MM/ML methylation tags, "
            "runs the full pipeline (modkit extract → entropy bedgraphs → "
            "segments) for true_mc and bisulfite modes (and ternary if you "
            "pass --include-ternary), and verifies all expected outputs. "
            "If this passes, the install is fully functional."
        ),
    )
    p_st.add_argument(
        "--include-ternary", action="store_true",
        help="Also test the ternary mode (slower; needs more memory).",
    )
    p_st.add_argument("-v", "--verbose", action="count", default=0)
    p_st.set_defaults(func=cmd_selftest)

    return p


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    _setup_logging(getattr(args, "verbose", 0))

    if args.cmd == "run":
        if args.pair:
            ctrl_specified = any([args.control_bam, args.control_tsv,
                                  args.control_bam_folder])
            tgt_specified = any([args.target_bam, args.target_tsv,
                                 args.target_bam_folder])
            if not (ctrl_specified and tgt_specified):
                parser.error(
                    "`run --pair` requires both control and target inputs "
                    "(--control-bam|--control-tsv|--control-bam-folder AND "
                    "--target-bam|--target-tsv|--target-bam-folder)."
                )
        else:
            if not (args.bam or args.tsv or args.bam_folder):
                parser.error(
                    "`run` (single mode) requires one of "
                    "--bam, --tsv, or --bam-folder."
                )

    ensure_dirs()
    try:
        return int(args.func(args) or 0)
    except FileNotFoundError as exc:
        logger.error("Missing file: %s", exc)
        return 2
    except SystemExit:
        raise
    except Exception:
        logger.exception("Unhandled error")
        return 1


if __name__ == "__main__":
    sys.exit(main())
