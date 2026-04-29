<img width="1024" height="506" alt="Gemini_Generated_Image_kcm5chkcm5chkcm5" src="https://github.com/user-attachments/assets/12957920-fdb1-4bdf-b704-989ddac66d00" />
<div align="center">

# shannonpore

**Nanopore methylation entropy analysis — GUI + CLI**

Per-CpG methylation entropy (ME), mean methylation level (MML), and coverage from
modkit-extracted nanopore reads. Three modes for handling 5hmC: True-mC,
Bisulfite-equivalent, and Ternary (3-state).

[![CI](https://github.com/uribertocchitau/shannonpore/actions/workflows/ci.yml/badge.svg)](https://github.com/uribertocchitau/shannonpore/actions/workflows/ci.yml)
[![Streamlit smoke](https://github.com/uribertocchitau/shannonpore/actions/workflows/streamlit-smoke.yml/badge.svg)](https://github.com/uribertocchitau/shannonpore/actions/workflows/streamlit-smoke.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-31012/)
[![Streamlit 1.51](https://img.shields.io/badge/streamlit-1.51-ff4b4b.svg)](https://streamlit.io)
[![modkit 0.6](https://img.shields.io/badge/modkit-0.6.0-orange.svg)](https://github.com/nanoporetech/modkit)
[![DOI](https://img.shields.io/badge/cite-CITATION.cff-blue.svg)](CITATION.cff)

</div>

---

## What it does

`shannonpore` turns a nanopore BAM with methylation tags (`MM`/`ML`) into
per-CpG entropy bedgraphs you can plot, segment, and compare across
samples — through either a polished Streamlit GUI or a Bash-friendly CLI.

```
              ┌──────────────┐  modkit extract   ┌──────────────┐
   sample.bam │     BAM      │ ───────────────►  │   per-CpG    │
              │  (MM/ML tags)│                   │     TSV      │
              └──────────────┘                   └──────┬───────┘
                                                        │
                                                        ▼
              ┌──────────────────────────────────────────────────┐
              │  whole-genome DuckDB pipeline (k CpGs / bin)     │
              │  ├─ true_mc    : 5hmC → unmodified C            │
              │  ├─ bisulfite  : 5hmC counted as 5mC            │
              │  └─ ternary    : 3-state (0=C, 1=5mC, 2=5hmC)   │
              └──────────────────────────────────────────────────┘
                                        │
                                        ▼
              ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐
              │ coverage │ │   MML    │ │    ME    │ │  + MhML  │ ◄ ternary
              │.bedgraph │ │.bedgraph │ │.bedgraph │ │.bedgraph │
              └──────────┘ └──────────┘ └──────────┘ └──────────┘
```

## Highlights

- **Two interfaces, one core** — Streamlit GUI for interactive exploration,
  Bash-friendly argparse CLI for HPC, Snakemake, Nextflow, or any shell.
- **Three input shapes** — single BAM, pre-computed modkit TSV, or a
  **folder of BAMs** that the tool auto-merges, sorts, and indexes
  before running modkit extract.
- **Single sample or pair mode** — process one sample, or control + target
  in one invocation. After a run, bedgraph paths auto-load into the
  Graph Preparation tab.
- **Three 5hmC modes** — True-mC (default), Bisulfite-equivalent, Ternary.
  Mode selection is a single radio button / `--mode` flag.



## Install

Requires conda/mamba/micromamba. The installer auto-bootstraps micromamba
into `~/.local/bin` if none is present.

```bash
git clone https://github.com/uribertocchitau/shannonpore.git
cd shannonpore
bash install.sh
```

The installer:

1. Resolves a conda CLI (or installs micromamba — no admin needed).
2. Creates the `shannonpore` environment from `conda-lock.yml`.
3. `pip install -e .` so `shannonpore` lands on `PATH`.
4. Installs Bash completion to `~/.local/share/bash-completion/completions/`.
5. Runs `shannonpore doctor` (every pinned version, every permission).
6. Runs `shannonpore selftest` (synthetic-BAM end-to-end pipeline).

If both pass, the install is fully functional. See [docs/INSTALL.md](docs/INSTALL.md)
for the full verification flow.

```bash
# After install:
conda activate shannonpore

# Run the GUI
streamlit run app.py
# or the CLI
shannonpore --help
```

**Reference data**: by default the tool looks in `<repo>/reference_files/`
for FASTAs + GTFs. Either populate that directory (run
`bash scripts/setup_references.sh --genome hg38`), pass `--fasta` /
"custom FASTA" per invocation, or set `SHANNONPORE_REF_DIR` to a folder
you already have. `selftest` doesn't need any of this — it builds its
own synthetic FASTA.

## Usage

Two interfaces, three input shapes, two sample designs.

| | single sample | pair (control + target) |
|---|---|---|
| **single BAM** | `--bam sample.bam` | `--pair --control-bam ... --target-bam ...` |
| **pre-computed TSV** | `--tsv sample.tsv` | `--pair --control-tsv ... --target-tsv ...` |
| **folder of BAMs** *(auto merge → sort → index)* | `--bam-folder /dir` | `--pair --control-bam-folder ... --target-bam-folder ...` |

After a successful run on the **File Preparation** tab, output bedgraph
paths are auto-loaded into the **Graph Preparation** tab — open it and
plot.

### CLI

```bash
# Print the walkthrough or cheat-sheet
shannonpore guide
shannonpore examples

# Single sample (BAM)
shannonpore run \
    --bam sample.bam --label patient42 \
    --out-dir results/ \
    --genome hg38 --mode true_mc --threads 16

# Folder of BAMs — merge + sort + index automatically
shannonpore run \
    --bam-folder /lab/runs/2026-04/bams \
    --label batch_2026_04 \
    --out-dir results/ \
    --genome hg38 --mode true_mc --threads 32

# Pair mode (control vs target) — BAM folders on both sides
shannonpore run --pair \
    --control-bam-folder /lab/runs/control \
    --target-bam-folder  /lab/runs/treated \
    --out-dir results/pair --genome hg38

# Just the modkit extract step (BAM or folder → TSV)
shannonpore extract --bam       sample.bam   sample.tsv  --genome hg38
shannonpore extract --bam-folder /lab/runs   merged.tsv  --genome hg38

# Plot
shannonpore plot hexbin scatter.png \
    --x-bedgraph control.me.bedgraph --y-bedgraph treated.me.bedgraph \
    --x-label "control ME" --y-label "treated ME"

# Verify install / debug
shannonpore doctor
shannonpore selftest --include-ternary
```

Tab completion (subcommands · options · mode/genome enums · file
paths) is installed automatically by `install.sh`. Live progress is
shown via tqdm bars; rerun with `-v` for streaming worker output.

### GUI

```bash
streamlit run app.py
```

Two tabs — `01 · file preparation` and `02 · graph preparation`. Both
ship inline `📖 How to use this tab` panels (auto-expanded for first-time
users), contextual `?` tooltips on every parameter, and a glossary
(CpG · k · ME · MML · 5mC · 5hmC). Long-running steps get a Streamlit
progress bar with live status.

## Three entropy modes

| Mode | 5hmC handling | When to use |
|---|---|---|
| **`true_mc`** _(default)_ | drop 5hmC calls (treat as unmodified C) | most analyses — isolates true 5mC signal |
| **`bisulfite`** | sum 5mC + 5hmC into one "methylated" signal | comparing to bisulfite-seq data (BS can't tell them apart) |
| **`ternary`** | 3 states per CpG (0=C, 1=5mC, 2=5hmC) | requires ≥ 3<sup>k</sup>× coverage (27× at k=3, 81× at k=4); produces an extra `MhML` bedgraph for 5hmC level |

The CLI surfaces a coverage warning when ternary's pattern-space exceeds
your `--min-coverage`. The GUI does the same in Tab 1.

## Project layout

```
v4/
├── app.py                       # Streamlit entrypoint (<200 lines)
├── bin/
│   ├── shannonpore              # bash launcher (works without pip install)
│   └── shannonpore.bash-completion
├── install.sh                   # idempotent installer
├── conda-lock.yml               # bit-for-bit reproducible env (generated)
├── environment.yml              # high-level conda spec
├── pyproject.toml               # project metadata + pip pins
├── requirements.txt             # pip-only deps
├── src/
│   ├── cli.py                   # argparse CLI: extract|entropy|plot|run|doctor|guide|examples|selftest
│   ├── config.py                # env-var resolved paths
│   ├── constants.py
│   ├── help_text.py             # shared user-facing prose (GUI + CLI)
│   ├── state.py                 # AppState · SampleSpec · TrackPlot
│   ├── tabs/                    # Streamlit tab modules
│   │   ├── tab_file_prep.py
│   │   └── tab_graph_prep.py
│   ├── pipelines/               # modkit · roi · whole_genome · ternary · bam_utils · orchestrator
│   ├── plots/                   # matplotlib/plotly: tracks, scatter, distributions
│   ├── io/                      # bedgraph, GTF, ROI, filesystem helpers
│   └── ui/                      # error_handler, progress, widgets, style
├── tests/
│   ├── unit/                    # 65+ unit tests
│   ├── integration/             # CLI, doctor/selftest, full-pipeline-on-BAM, pair, folder-merge
│   ├── e2e/                     # Streamlit AppTest smoke
│   └── fixtures/                # synthetic FASTA + BAM-with-MM/ML-tags builder
└── docs/
    ├── ARCHITECTURE.md
    ├── INSTALL.md
    ├── MIGRATION_FROM_V3.md
    └── SESSION_STATE_SCHEMA.md
```


## Citation

If you use `shannonpore` in published work, please cite via the
[`CITATION.cff`](CITATION.cff) file (GitHub renders a "Cite this repository"
button on the right sidebar).

## License

[MIT](LICENSE) © 2026 Uri Bertocchi, Ebenstein Lab, Tel Aviv University.

## Acknowledgements

- [modkit](https://github.com/nanoporetech/modkit) (Oxford Nanopore) —
  per-read CpG methylation calls.
- [Streamlit](https://streamlit.io) — interactive frontend.
- [DuckDB](https://duckdb.org) — out-of-core SQL for whole-genome ingest.
- [pyfaidx](https://github.com/mdshw5/pyfaidx) — FASTA indexing.
- The Ebenstein Lab for the science.
