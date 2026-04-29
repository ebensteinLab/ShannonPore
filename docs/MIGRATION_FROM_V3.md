# Migrating from v3 to v4

v3 lives at `../v3/` and is **untouched** by this refactor. Use it as a
fallback while you adopt v4. Only delete v3 (or rename to
`v3-archived/`) once you have validated v4 against your full workflow.

## What stays the same

- **Reference files** — v4 reads from `NANOENTROPY_REF_DIR` which by
  default points at the v3 `reference_files/` directory. No multi-GB
  FASTAs are duplicated.
- **Output layout** — v4 creates `results/` with the same per-run
  subdirectory pattern (`bedgraphs/`, `wg_work/`).
- **Bedgraph schema** — still `chrom \t start \t end \t value`.
- **modkit invocation** — same `--cpg --threads N --force` flags.

## What changed

### 1. Install is reproducible

v3 had no `requirements.txt`. v4 ships:

- `pyproject.toml` (pinned Python deps with `==`)
- `environment.yml` (conda spec; Python + modkit + system libs)
- `conda-lock.yml` (bit-for-bit lockfile)
- `install.sh` (one-command installer)
- `nanoentropy doctor` + `nanoentropy selftest` for verification

### 2. Code is modular

- `app.py` shrunk from 7,450 lines to <200.
- 2 tabs in `src/tabs/` (was 3 in v3).
- Plotting in `src/plots/` (one file per plot family).
- Pipelines in `src/pipelines/`.

### 3. State is typed

107+ `st.session_state[...]` accesses → one `AppState` dataclass.
See [SESSION_STATE_SCHEMA.md](SESSION_STATE_SCHEMA.md).

### 4. Errors are handled centrally

44 bare `except Exception:` → `@show_error` decorator that logs the full
traceback and surfaces a friendly message in the UI.

### 5. Three entropy modes in Tab 1 / CLI

| Mode | Behaviour |
|---|---|
| **True-mC** _(default)_ | 5hmC dropped; only 5mC drives entropy. |
| **Bisulfite** | 5hmC counted as 5mC (matches bisulfite seq). |
| **Ternary** | 3 states (C/5mC/5hmC); requires ≥ 3<sup>k</sup>× coverage. |

The ternary pipeline is `src/pipelines/ternary_entropy.py`, ported from
`5hmc_confound_analysis/analysis/ternary_entropy.py`.

### 6. CLI — `nanoentropy` console script

v3 had no CLI; the Streamlit GUI was the only interface. v4 ships a
Bash-friendly CLI with subcommands `extract`, `entropy`, `plot`, `run`,
`doctor`, `selftest`. Tab completion supplied by the bash-completion
file in `bin/`.

### 7. GUI redesigned

Generic Streamlit chrome replaced with a restrained scientific theme:
custom `.streamlit/config.toml`, IBM Plex fonts, hairline rules, flat
buttons, hidden MainMenu/footer, monospace metrics. Section headers
follow a `01 · input source / 02 · reference / ...` numbered convention.

### 8. Removed: Data Analysis tab

The v3 "Data Analysis" tab (block segmentation + ChIPseeker /
GeneHancer / gprofiler2 enrichment) is out of scope for v4. The 1,456-line
embedded-R-string is gone, and the R bridge has been removed. If you
need block segmentation or GO enrichment, run them separately on the
exported bedgraphs.

### 9. Bug fixes

- `build_parquet_from_tsv_stream_duckdb` (v3) referenced undefined
  `force` and `stream_cb` — fixed in v4 by adding them as parameters.
- scipy bumped 1.8.0 → 1.11.4 to match numpy 1.26.4.
- Whole-genome ingest no longer raises `SystemExit` on missing
  duckdb/pyfaidx.

### 10. Tests exist

83 tests: unit (60), integration (CLI, doctor, selftest, full BAM
pipeline), e2e (Streamlit smoke). The synthetic-BAM fixture
(`tests/fixtures/build_tiny_bam.py`) generates real `MM`/`ML`-tagged
BAMs that modkit extracts from.

## Migration checklist

1. `cd v4 && bash install.sh` — produces `nanoentropy` conda env.
2. `cp .env.example .env` and edit if your reference paths differ.
3. `source .env && conda activate nanoentropy`.
4. `nanoentropy doctor` (must exit 0).
5. `nanoentropy selftest` (must exit 0).
6. `streamlit run app.py` — verify both tabs render.
7. Reproduce one of your v3 runs: same BAM + same genome → same
   bedgraphs (numerical tolerance ~1e-6).
8. `pytest -v` — must be green.
9. Once happy: rename `../v3/` to `../v3-archived/`.

## Rolling back

If v4 misbehaves:

- v3 is fully intact. `cd ../v3 && streamlit run app.py` resumes the v3
  workflow.
- v4 results live in `v4/results/`; they are independent of v3.
- The conda env `nanoentropy` can be removed without affecting v3:
  `conda env remove -n nanoentropy`.
