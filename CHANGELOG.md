# Changelog

## Unreleased

### Fixed
- **bioconda package name** — `environment.yml` pinned `modkit=0.6.0`
  but bioconda's recipe is `ont-modkit` (the on-disk binary is still
  `modkit`). Mamba aborted with `modkit =0.6.0 * does not exist
  (perhaps a typo or a missing channel)` and the entire install
  cascaded: env was never created, every subsequent `mamba run -n
  nanoentropy_v4 …` failed with "Environment must first be created…".
  Renamed to `ont-modkit=0.6.0`.
- **install.sh hard-fails on env-create error** — previously errors
  silently passed through and downstream steps emitted confusing
  "Environment must first be created" messages. Now we abort with
  a clear message at the source.
- `install.sh` no longer auto-downloads ~6 GiB of FASTAs by default —
  this hung most fresh installs. Reference downloads are now opt-in
  via `NANOENTROPY_DOWNLOAD_REFERENCES=1|hg38|mm10`, or by running
  `scripts/setup_references.sh` manually. See `docs/INSTALL.md`.
- Added `pysam==0.22.1` to `environment.yml`, `requirements.txt`,
  `pyproject.toml`, and the doctor's pinned-package check. `selftest`
  imports pysam to build its synthetic BAM, so without this it would
  fail immediately after a clean install.

## v0.1.0 — 2026-04-29

### Added

- **Three input shapes** for both GUI and CLI: single BAM, pre-computed
  modkit TSV, or a folder of BAMs that is auto-merged → sorted → indexed
  via `samtools` before modkit extract.
- **Pair mode** — process control + target samples in one invocation,
  producing two parallel sets of bedgraphs.
- **Auto-population** of the Graph Preparation tab with the bedgraph
  paths produced by the most-recent run on File Preparation.
- **Live progress reporting** — tqdm bars in the CLI, `st.progress` in
  the GUI; the per-chromosome stage streams "X/Y chroms" updates via
  `imap_unordered`.
- **In-app instructions** — `nanoentropy guide` (walkthrough) and
  `nanoentropy examples` (cheat-sheet) CLI subcommands;
  `📖 How to use this tab` expanders, contextual tooltips, and a
  glossary in the GUI.
- Reproducible install via `conda-lock` (Python 3.10.12, modkit 0.6.0,
  samtools 1.21).
- Centralized session state (`src/state.py` `AppState` dataclass).
- Modular tab architecture (`src/tabs/`).
- Plot library decoupled from Streamlit callbacks (`src/plots/`).
- Full QA test pyramid (83 tests: unit + integration + Streamlit E2E +
  full BAM pipeline).
- Environment-variable-resolved paths (no hardcoded lab paths).
- `@show_error` decorator replaces all bare `except Exception:` clauses.
- **CLI** — `nanoentropy` console script with 6 subcommands:
  `extract`, `entropy`, `plot`, `run`, `doctor`, `selftest`. Mirrors
  every pipeline feature of the Streamlit GUI.
- **3 entropy modes** in the GUI and CLI: True-mC (default; 5hmC → C),
  Bisulfite (5hmC counted as 5mC), Ternary (3-state, requires 3<sup>k</sup>×
  coverage).
- **`nanoentropy doctor`** — post-install audit of every dependency,
  version, permission, and reference-data path.
- **`nanoentropy selftest`** — synthesises a tiny FASTA + BAM with
  MM/ML tags and runs the full pipeline end-to-end across all modes.
- **Bash completion** for the CLI (subcommands, options, mode/genome
  enums, file/dir paths). Installed automatically by `install.sh`.
- **`bin/nanoentropy`** — bash launcher that works without
  `pip install -e .` (auto-detects the conda env).
- **Custom Streamlit theme** — `.streamlit/config.toml` + `src/ui/style.py`
  for a restrained scientific aesthetic (IBM Plex fonts, hairline rules,
  flat buttons, hidden Streamlit chrome).
- **GitHub project files** — issue/PR templates, CI workflows, security
  policy, code of conduct, citation file, dependabot config.
- **Adversarial review pipeline** (optional) — 5-agent ClawTeam workflow
  for code/UX review under `clawteam/`.
- **MIT License**.

### Removed

- Tab 2 "Data Analysis" feature (segmentation + GO/GeneHancer enrichment)
  and its R bridge — out of scope for v4.

### Notes

- v3 remains untouched as fallback during transition.
