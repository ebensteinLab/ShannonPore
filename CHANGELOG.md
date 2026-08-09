# Changelog

## Unreleased

### Package rename: `src` → `shannonpore` (BREAKING for imports)

- The top-level Python package is now `shannonpore` (was the generic
  `src`, which collides with any other package that ships a `src`
  module and blocks bioconda acceptance). All imports change
  accordingly: `from shannonpore.plots.scatter import …`; the module
  CLI is `python -m shannonpore.cli`. The installed `shannonpore`
  console script, GUI, and all pipeline behaviour are unchanged.
  Re-run `pip install -e .` in existing dev checkouts.

### Multi-format plot export + paired-landscape line toggle

**Features**
- Every plot (ME/MML scatter, arch landscape, paired landscape, region
  tracks) can now be exported as high-resolution **PNG / JPG / SVG /
  PDF** in one render. New shared helper `shannonpore.plots.export.save_figure`
  centralises saving: raster formats honour a configurable DPI
  (bounded 30–1200), JPG is flattened onto white at quality 95, and
  SVG/PDF are true vector output.
- GUI: new **export settings** row on the Graph Preparation tab —
  format multiselect + DPI input (default PNG @ 300) applied to all
  four render buttons; every saved file is listed with its size.
- GUI: new **connecting lines** toggle on the paired landscape hides
  the black control → target `LineCollection` (`show_lines=False`),
  keeping only the per-sample dots.
- CLI: `shannonpore plot` gains `--formats png,jpg,svg,pdf` (validated,
  alias-deduped), `--dpi N` (bounded, friendly argparse errors), and
  `--lines/--no-lines` for the paired landscape.

**Back-compat**
- `out_path`-only calls (no `formats=`) keep the old behaviour exactly,
  including pass-through of any matplotlib-supported suffix (`.eps`,
  `.tif`, …) and the dpi=200 default. `plot_region_tracks` still
  returns a single `Path`.

**QC**
- +22 tests (export helper unit tests, line-toggle assertions,
  multi-format integration through the CLI). Security review closed a
  resource-exhaustion gap by bounding DPI at both the library and
  argparse layers.

### Repo-wide coherence sweep (3 parallel reviewers)

**Code coherence**
- `src/cli.py` module docstring no longer lists removed `segment` /
  `annotate` subcommands; surface the actual 9 (`extract`, `entropy`,
  `plot`, `run`, `prefetch`, `doctor`, `selftest`, `guide`, `examples`).
- `plot landscape` CLI help text dropped the stale "direction arrows"
  claim; describes the caption-based bin counts (Δ ME ↑ / ↓) instead.
- `plot tracks` CLI description updated to match the 4-panel layout.
- Color palette centralized via `src.constants.PALETTE_CONTROL` /
  `PALETTE_TARGET` (`#2980b9` / `#e67e22`). Single source of truth used
  by `state.TrackPlot` defaults, `GraphPrepState` factories,
  `tab_file_prep` widget initializers, `cli.py --color-{a,b}` defaults,
  and `scatter.py` plot fn defaults. Old `#1f77b4` / `#ff7f0e` literals
  removed everywhere.
- `scatter.py` module docstring updated for the caption-based paired
  landscape.

**Docs drift**
- `docs/INSTALL.md`: corrected the "FASTAs not downloaded by default"
  claim (install.sh now downloads by default); reordered the section
  to lead with the default behaviour, listed all four reference-data
  paths (default, skip, single-genome, post-install prefetch), removed
  the unsupported `hg19` from the bash-completion enum list.
- `README.md`: flipped the reference-fetching narrative — install.sh
  default is auto-download, lazy fetch is the always-on fallback,
  `shannonpore prefetch` is the explicit pre-warm.
- Deleted vestigial `.R-version` (R was removed in v4).

**CI + packaging**
- `.github/workflows/ci.yml` now runs `black --check src tests app.py`
  alongside `ruff check`. CONTRIBUTING.md's claimed style policy is
  finally enforced in CI.
- Applied `black` over the whole codebase (33 files reformatted) so
  the new check passes on day one.
- `environment.yml` historical comment about `modkit=0.6.0` recipe
  trimmed to the load-bearing fact (the bioconda recipe is
  `ont-modkit`, the binary is `modkit`).
- `streamlit-smoke.yml` already exercises `AppTest.from_file("app.py")`
  — no change needed.

**Dependency bumps (rolled in from Dependabot PR #9)**
- `streamlit` 1.51.0 → 1.57.0
- `duckdb` 1.4.3 → 1.5.2
- `matplotlib` 3.8.4 → 3.10.9
- `pyfaidx` 0.8.1.2 → 0.9.0.4
- `pysam` 0.22.1 → 0.24.0
- `plotly` 6.2.0 → 6.7.0
- `tqdm` 4.67.1 → 4.67.3
- `scipy` 1.11.4 → 1.15.3
- `ruff` 0.7.4 → 0.15.12
- `mypy` 1.13.0 → 1.20.2
- Doctor's pinned-version map, environment.yml, requirements.txt, and
  pyproject.toml all aligned. `pytest` (132/132), `ruff check`, and
  `black --check` all pass on the bumped stack.

### Added — graph prep & references
- **Region track plot** rewritten to a 4-panel layout: gene structure
  (exons, 1 kb promoter, strand arrows) → smoothed ME → smoothed MML →
  optional smoothed coverage. Multiple genes in a window stack onto
  separate rows automatically.
- **Bundled GTFs auto-load.** hg38 and mm10 RefSeq GTFs are downloaded
  on first use (~30–40 MB each) via `ensure_genome_gtf()`; the gene
  panel and gene-name search "just work" without `setup_references.sh`.
- **Gene-name search** in the Graph Prep track section: type a symbol
  (e.g. `VHL`), click FIND, and chrom/start/end auto-fill. Near-misses
  surface as "did you mean…" hints.
- **Coverage panel** in the track plot when control/target coverage
  bedgraphs are supplied. Autoscales (no 0–1 clamp).
- **Saved-plot paths** — every render shows the on-disk PNG path with
  file size; the tab header lists the base output directory.
- **ME / MML scatter, arch landscape, paired landscape** plot families
  (replacing the old hexbin / scatter / distribution UI). Paired
  landscape draws bin counts (Δ ME ↑ / Δ ME ↓) in a caption block
  below the axes — no in-plot arrows.
- CLI `plot tracks` learns `--control-coverage` / `--target-coverage`,
  `--genome`, `--smooth`/`--window`, `--pad`. `plot scatter|arch|landscape`
  accept matched control/target bedgraphs and per-plot knobs.

### Fixed
- **GUI state refresh** — the Graph Prep tab was using a one-cycle-stale
  AppState after `update_section`, so genome / GTF changes only took
  effect on the *next* widget interaction. Now refetches via `get_state()`.
- **paired_landscape performance** — replaced the per-row `ax.plot` loop
  with a single `LineCollection`, dropping render time on 20 k bins
  from seconds to ~100 ms.
- **vmax_pct ignored in linear colour scale** for me_mml_scatter and
  triple_landscape — both modes now respect the parameter.
- **track plot edge-clipped bins** — `_read_region` now uses half-open
  overlap so a bedgraph bin that straddles the window edge isn't
  silently dropped.
- **Path-traversal hardening** — sanitize `chrom` (track-plot filename)
  and `SampleSpec.label` (output-path component) before use.
- **modkit zombies** — terminate path now `wait()`s after `terminate()`,
  with `kill()` fallback.
- **DuckDB SQL injection-resistant TSV path** — quote single-quotes in
  TSV paths in `read_csv_auto('…')` so paths like `/data/O'Brien/x.tsv`
  parse correctly.
- **GTF download stalls** — global socket timeout for the read loop so a
  hung mirror can't block forever; failed downloads now also clean up
  non-zero partials.
- **Orchestrator** — `force_ingest=True` from `--force` is now actually
  forwarded to the entropy pipelines; `run_pipeline`'s `min_coverage`
  default aligned with the GUI/CLI default of 16.
- **`open(p, "w").close()` FD leaks** — replaced with `Path.touch()` in
  the empty-output paths of both whole-genome and ternary pipelines.
- **`bioconda package name`** — `environment.yml` pinned `modkit=0.6.0`
  but bioconda's recipe is `ont-modkit`. Renamed to `ont-modkit=0.6.0`.
- **`install.sh` micromamba bootstrap** — added `--fail --retry 3` to
  the curl call so a non-2xx response aborts instead of feeding HTML
  into tar.
- **`install.sh` hard-fails on env-create error** instead of silently
  passing through to a confusing "Environment must first be created"
  downstream.
- Reference FASTAs are downloaded **by default** by `install.sh`. Skip
  with `SHANNONPORE_SKIP_REFERENCES=1`; pin to one genome with
  `SHANNONPORE_DOWNLOAD_REFERENCES=hg38|mm10`. (GTFs are separately
  lazy-downloaded on first GUI use.)
- Added `pysam==0.22.1` to `environment.yml`, `requirements.txt`,
  `pyproject.toml`, and the doctor's pinned-package check.

### Tests
- New unit tests for `GeneStructure`, `load_gene_structures` (caching
  + promoter strand logic), `find_gene_by_name`, `search_gene_names`,
  `ensure_genome_gtf` (mocked download success / failure cleanup),
  `me_mml_scatter` / `triple_landscape` / `paired_landscape` (including
  no-bins-pass-filter and subsample-aware caption), and
  `plot_region_tracks` 3-panel + 4-panel branches.
- New integration tests for `plot tracks` CLI — both with and without
  the coverage panel.

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
- **In-app instructions** — `shannonpore guide` (walkthrough) and
  `shannonpore examples` (cheat-sheet) CLI subcommands;
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
- **CLI** — `shannonpore` console script with 6 subcommands:
  `extract`, `entropy`, `plot`, `run`, `doctor`, `selftest`. Mirrors
  every pipeline feature of the Streamlit GUI.
- **3 entropy modes** in the GUI and CLI: True-mC (default; 5hmC → C),
  Bisulfite (5hmC counted as 5mC), Ternary (3-state, requires 3<sup>k</sup>×
  coverage).
- **`shannonpore doctor`** — post-install audit of every dependency,
  version, permission, and reference-data path.
- **`shannonpore selftest`** — synthesises a tiny FASTA + BAM with
  MM/ML tags and runs the full pipeline end-to-end across all modes.
- **Bash completion** for the CLI (subcommands, options, mode/genome
  enums, file/dir paths). Installed automatically by `install.sh`.
- **`bin/shannonpore`** — bash launcher that works without
  `pip install -e .` (auto-detects the conda env).
- **Custom Streamlit theme** — `.streamlit/config.toml` + `src/ui/style.py`
  for a restrained scientific aesthetic (IBM Plex fonts, hairline rules,
  flat buttons, hidden Streamlit chrome).
- **GitHub project files** — issue/PR templates, CI workflows, security
  policy, code of conduct, citation file, dependabot config.
- **MIT License**.

### Removed

- Tab 2 "Data Analysis" feature (segmentation + GO/GeneHancer enrichment)
  and its R bridge — out of scope for v4.

### Notes

- v3 remains untouched as fallback during transition.
