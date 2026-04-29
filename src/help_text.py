"""Shared user-facing help text.

Single source of truth for prose used by:
- the GUI (`Getting started` expander, widget tooltips, glossary)
- the CLI (`shannonpore guide`, `shannonpore examples`)

Keep these strings short and copy-edited; both surfaces render them
verbatim.
"""

from __future__ import annotations

GUIDE = """\
# shannonpore — quick start

`shannonpore` turns a nanopore BAM (with MM/ML methylation tags) into
per-CpG methylation entropy bedgraphs you can plot, segment, and
compare across samples — through either a polished Streamlit GUI or a
Bash-friendly CLI.

## 1.  Pick your input shape

Every run produces one or two **samples**. A sample is one of:

  * **single BAM**      — `--bam sample.bam`
  * **single TSV**      — `--tsv sample_modkit_extract.tsv`
                         (skips modkit extract; you already have a TSV)
  * **folder of BAMs**  — `--bam-folder /path/to/bams_dir`
                         (auto-merge → sort → index → modkit extract)

Pair mode (`--pair`) processes a control AND a target sample in one
invocation, producing two parallel sets of bedgraphs.

## 2.  Pick your entropy mode

`--mode` controls how 5hmC calls are interpreted:

  * **true_mc**    (default)  — drop 5hmC; only true 5mC drives entropy
  * **bisulfite**             — count 5hmC as 5mC (matches BS-seq)
  * **ternary**               — 3-state (0=C, 1=5mC, 2=5hmC).
                                 Requires ≥ 3^k× coverage.

## 3.  Pick bin size + coverage

  * `--cpgs-per-bin k`   — k consecutive CpGs per entropy bin (default 4)
  * `--min-coverage N`   — minimum reads per bin to emit a value (default 16)
  * `--methyl-threshold` — probability threshold for "methylated" (default 0.5)

## 4.  Run

    shannonpore run --bam sample.bam --out-dir results/ \\
        --genome hg38 --mode true_mc --threads 16

Outputs land at `<out-dir>/<label>_<mode>.{coverage,me,mml}.bedgraph`.
A `run_summary.json` file in the output dir lists every produced file.

## 5.  Plot

CLI:

    shannonpore plot scatter scatter.png \\
        --x-bedgraph control.me.bedgraph \\
        --y-bedgraph target.me.bedgraph \\
        --x-label "control ME" --y-label "target ME"

GUI:
The Plotting tab is auto-populated with the bedgraphs from your last run.
Just open it and click "render".

## Verification & troubleshooting

  * `shannonpore doctor`    — verify every dependency, version, permission
  * `shannonpore selftest`  — run a synthetic BAM through the full pipeline
  * `shannonpore --help`    — list all subcommands
"""


EXAMPLES = """\
# shannonpore — common recipes

# Single sample, BAM input
shannonpore run \\
    --bam sample.bam \\
    --out-dir results/sample \\
    --genome hg38 --mode true_mc \\
    --cpgs-per-bin 4 --min-coverage 16 --threads 16

# Single sample, label-prefixed outputs
shannonpore run \\
    --bam sample.bam --label patient42 \\
    --out-dir results/ \\
    --genome hg38 --mode true_mc

# Folder of BAMs (auto merge+sort+index)
shannonpore run \\
    --bam-folder /lab/runs/2026-04/bams \\
    --label batch_2026_04 \\
    --out-dir results/ \\
    --genome hg38 --mode true_mc --threads 32

# Pair mode (control vs target)
shannonpore run --pair \\
    --control-bam ctrl.bam --control-label control \\
    --target-bam tgt.bam   --target-label treated \\
    --out-dir results/pair --genome hg38 --mode true_mc

# Pair mode, both sides are folders of BAMs
shannonpore run --pair \\
    --control-bam-folder /lab/runs/control \\
    --target-bam-folder  /lab/runs/treated \\
    --out-dir results/pair --genome hg38

# Ternary mode (5hmC-aware) — needs high coverage (3^k)
shannonpore run \\
    --bam sample.bam --out-dir results/ \\
    --genome hg38 --mode ternary \\
    --cpgs-per-bin 3 --min-coverage 27

# Just modkit extract (BAM → TSV)
shannonpore extract --bam sample.bam sample.tsv --genome hg38 --threads 8

# Just entropy from a pre-computed TSV
shannonpore entropy sample.tsv out/sample \\
    --genome hg38 --mode bisulfite --cpgs-per-bin 4

# Plot density scatter of two samples
shannonpore plot hexbin scatter.png \\
    --x-bedgraph control.me.bedgraph \\
    --y-bedgraph treated.me.bedgraph \\
    --x-label "control ME" --y-label "treated ME"

# Distribution of ME values
shannonpore plot distribution dist.png \\
    --x-bedgraph sample.me.bedgraph --dist-kind violin

# Verify install
shannonpore doctor
shannonpore selftest
"""


# ─── GUI-only help snippets ───────────────────────────────────────────────

GUI_GETTING_STARTED = """\
**shannonpore** computes per-CpG methylation entropy from nanopore reads.
You feed it a BAM (or a folder of BAMs, or a pre-computed modkit TSV) and
it produces three bedgraphs you can plot in the next tab:

| File | Meaning |
|---|---|
| `*.coverage.bedgraph` | reads per bin |
| `*.mml.bedgraph`      | mean methylation level (0–1) per bin |
| `*.me.bedgraph`       | normalised methylation entropy (0–1) per bin |

**Three steps:**

1. **File Preparation** *(this tab)* — pick inputs, choose a 5hmC mode, run.
2. **Graph Preparation** *(next tab)* — render tracks, scatter, distributions.
   Bedgraph paths are auto-loaded after step 1.

You can run a **single sample**, or **pair mode** (control + target) which
produces two parallel sets of bedgraphs ready for side-by-side plots.
"""

GUI_GLOSSARY = """\
| term | meaning |
|---|---|
| **CpG** | a cytosine immediately followed by a guanine in the reference |
| **k**   | number of consecutive CpGs that form one entropy bin |
| **MML** | *mean methylation level* — fraction of reads × CpGs in the bin that are methylated, ∈ [0,1] |
| **ME**  | *methylation entropy* — Shannon entropy over per-read methylation patterns in the bin, normalised to [0,1] |
| **MhML** | mean **5h**ydroxy-methyl level (ternary mode only) |
| **methyl threshold** | probability ≥ threshold ⇒ "methylated"; below ⇒ "unmodified" |
| **5mC**  | 5-methylcytosine (canonical methylation) |
| **5hmC** | 5-hydroxymethylcytosine (oxidised intermediate) |
"""

GUI_HELP_INPUT_KIND = (
    "single BAM = one file. "
    "pre-computed TSV = skip modkit extract (you already ran it). "
    "folder of BAMs = auto merge + sort + index, then extract."
)

GUI_HELP_PAIR_MODE = (
    "On: process two samples (control + target) in one run. "
    "Off: process one sample. "
    "Both sides can independently be a BAM, a TSV, or a folder of BAMs."
)

GUI_HELP_K = (
    "How many consecutive CpGs form one entropy bin. "
    "Larger k → richer pattern space (2^k binary, 3^k ternary) but "
    "needs more coverage."
)

GUI_HELP_MIN_COV = (
    "Minimum reads covering a bin to emit a metric. "
    "Bins below this are written with no value (empty). "
    "Default 16 — works well for most nanopore datasets."
)


# ── Per-section descriptions for the File Preparation tab ────────────────
SEC_FP_SAMPLES = (
    "Pick what to feed in: a single BAM, a pre-computed modkit TSV, or a "
    "folder of BAMs that the tool will merge → sort → index for you. "
    "Toggle pair mode to process control + target side-by-side."
)
SEC_FP_REFERENCE = (
    "Choose hg38 or mm10 — modkit aligns CpG positions against this FASTA. "
    "Override with a custom FASTA if your samples were called against a "
    "different reference."
)
SEC_FP_ENTROPY_MODE = (
    "Decides how 5-hydroxymethylcytosine (5hmC) is handled. True-mC is the "
    "safe default; pick bisulfite-equivalent only when comparing to BS-seq, "
    "and ternary only when you have ≥ 3^k× coverage."
)
SEC_FP_BIN_PARAMS = (
    "k controls how many CpGs make one bin (richer pattern space, more "
    "coverage needed). min coverage filters bins with too few reads. "
    "Threads parallelise the per-chromosome pass."
)
SEC_FP_OUTPUT = (
    "Where the run writes its bedgraphs, the modkit TSV, the merge log, "
    "and `run_summary.json`. Each run is idempotent — re-running with the "
    "same inputs reuses cached intermediates."
)


# ── Per-section descriptions for the Graph Preparation tab ───────────────
SEC_GP_SAMPLES = (
    "Each side (control / target) takes the three bedgraphs produced by "
    "File Prep. Paths auto-fill after a successful run; you can also paste "
    "them in by hand to plot prior runs."
)
SEC_GP_REGION = (
    "Pick a single contig + interval to draw track plots. Smaller windows "
    "(<1 Mb) render fastest. The optional GTF turns on a gene panel under "
    "the methylation tracks."
)
SEC_GP_SCATTER = (
    "Compare control vs target for one metric (MML or ME) across the whole "
    "genome. Hexbin shows density; scatter shows individual bins (subsample "
    "for speed on dense data)."
)
SEC_GP_DISTRIBUTION = (
    "Distribution of values from one bedgraph across all bins. Useful for "
    "spotting bimodal MML or wide ME spreads at a glance."
)

GUI_HELP_METHYL_THRESHOLD = (
    "Per-call probability threshold. mod_qual > threshold ⇒ 'methylated'. "
    "0.5 is the modkit default."
)
