<p align="center">
  <a href="https://github.com/uribertocchitau/shannonpore">
    <img src="https://github.com/user-attachments/assets/12957920-fdb1-4bdf-b704-989ddac66d00"
         alt="shannonpore — methylation entropy from nanopore reads"
         width="100%">
  </a>
</p>

<p align="center">
  <sub>
    Per-CpG methylation entropy (ME), mean methylation level (MML), and coverage from modkit-extracted nanopore reads.<br/>
    Three 5hmC modes: <b>True-mC</b> · <b>Bisulfite-equivalent</b> · <b>Ternary</b>.
  </sub>
</p>

<p align="center">

[![CI](https://github.com/uribertocchitau/shannonpore/actions/workflows/ci.yml/badge.svg)](https://github.com/uribertocchitau/shannonpore/actions/workflows/ci.yml)
[![Streamlit smoke](https://github.com/uribertocchitau/shannonpore/actions/workflows/streamlit-smoke.yml/badge.svg)](https://github.com/uribertocchitau/shannonpore/actions/workflows/streamlit-smoke.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-31012/)
[![Streamlit 1.51](https://img.shields.io/badge/streamlit-1.51-ff4b4b.svg)](https://streamlit.io)
[![modkit 0.6](https://img.shields.io/badge/modkit-0.6.0-orange.svg)](https://github.com/nanoporetech/modkit)
[![DOI](https://img.shields.io/badge/cite-CITATION.cff-blue.svg)](CITATION.cff)

</p>

<br/>

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
# After install — call from any shell, no activation needed:
shannonpore-gui              # launch the Streamlit GUI
shannonpore --help           # the CLI
```

`install.sh` symlinks both wrappers into `~/.local/bin/`, so they're
global commands. They use `mamba run` under the hood, which means you
**never** need to `mamba activate shannonpore` — the env is resolved
automatically.

**Reference data**: bundled hg38 + mm10 FASTAs + RefSeq GTFs are
downloaded automatically by `install.sh` (multi-GiB; 20–60 min). Anything
still missing is fetched lazily on first use — type a region or hit
**Run** in the GUI and the GTF / FASTA will pull from UCSC on demand
(~40 MiB GTF, ~1 GiB FASTA → ~3 GiB unzipped + `samtools faidx`).

To pre-fetch ahead of time: `shannonpore prefetch hg38 mm10`. To skip
the install-time fetch: `SHANNONPORE_SKIP_REFERENCES=1 bash install.sh`.
You can also pass `--fasta` (CLI) or "custom FASTA" (GUI) per invocation,
or set `SHANNONPORE_REF_DIR` to point at an existing shared directory.
`selftest` doesn't need any of this — it builds its own synthetic FASTA.

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

# Plot — four kinds. All take the same control + target bedgraphs.
shannonpore plot scatter   scatter.png   --control-mml ctrl.mml.bedgraph --control-me ctrl.me.bedgraph --target-mml tgt.mml.bedgraph --target-me tgt.me.bedgraph --label-a Control --label-b Treated
shannonpore plot arch      arch.png      --control-mml ...               --control-me ...              --target-mml ...              --target-me ...
shannonpore plot landscape paired.png    --control-mml ...               --control-me ...              --target-mml ...              --target-me ...               --filter-a-dim '|dMML|' --filter-a-op '<' --filter-a-value 0.1 --filter-b-dim '|dME|' --filter-b-op '>' --filter-b-value 0.4
shannonpore plot tracks    tracks.png    --control-mml ...               --control-me ...              --target-mml ...              --target-me ...               --control-coverage ctrl.coverage.bedgraph --target-coverage tgt.coverage.bedgraph --genome hg38 --chrom chr3 --start 10141778 --end 10153676

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

Two tabs:

* **`01 · file preparation`** — pick BAM(s) / TSV / folder, choose
  entropy mode (true_mc / bisulfite / ternary) and bin parameters,
  run the pipeline. Outputs land at `<results>/file_prep/`.
* **`02 · graph preparation`** — five sections in order:
  `01 · samples → 02 · ME / MML scatter → 03 · arch landscape →
  04 · paired landscape → 05 · region (track plot)`. Bedgraph paths
  are auto-populated from the last File Prep run. The track section
  loads the bundled hg38 / mm10 RefSeq GTF on first use and lets you
  search by gene name to auto-fill `chrom / start / end`.

Both tabs ship inline `📖 How to use this tab` panels (auto-expanded
for first-time users), contextual `?` tooltips, and a glossary
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
