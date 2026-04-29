# Install & verification — nanoentropy v4

This page documents the full install flow. Goal: a user who downloads
v4 today and runs `bash install.sh` should have a fully working GUI +
CLI with no extra steps. The same install procedure run in 2027 should
reproduce the 2026 stack bit-for-bit.

## TL;DR

```bash
git clone https://github.com/uribertocchitau/nanoentropy.git && cd nanoentropy
bash install.sh
conda activate nanoentropy
nanoentropy doctor      # auto-checks every dep + permission
nanoentropy selftest    # synthetic-BAM end-to-end smoke test
```

## Reference data (multi-GiB)

The default `bash install.sh` does **not** download reference FASTAs —
they're multi-GiB and hang most installs. The tool needs a populated
`REFERENCE_DIR` only when you point it at real BAM/TSV data; `selftest`
ships its own synthetic FASTA so you can verify the install end-to-end
without fetching anything multi-GiB.

`REFERENCE_DIR` defaults to `<repo>/reference_files/`. You only need
to set `NANOENTROPY_REF_DIR` if you want it to point somewhere else
(e.g. a shared lab volume). Three ways to populate it:

```bash
# 1. Have install.sh fetch into the default location
NANOENTROPY_DOWNLOAD_REFERENCES=hg38 bash install.sh   # ~3.2 GiB
NANOENTROPY_DOWNLOAD_REFERENCES=mm10 bash install.sh   # ~2.8 GiB
NANOENTROPY_DOWNLOAD_REFERENCES=1    bash install.sh   # both

# 2. Run install.sh now (skips refs), fetch later
bash install.sh
bash scripts/setup_references.sh --genome hg38

# 3. Use FASTAs you already have at a different path
NANOENTROPY_REF_DIR=/path/to/your/refs bash install.sh
```

You can also bypass `REFERENCE_DIR` entirely and pass `--fasta` to the
CLI or use the "custom FASTA" field in the GUI per invocation.

## What `install.sh` does

1. **Resolves a conda CLI** — uses `micromamba` if present; else `mamba`,
   else `conda`. If none available, downloads micromamba into
   `~/.local/bin` (no admin rights needed).
2. **Creates the env** — from `conda-lock.yml` (preferred; bit-for-bit)
   or `environment.yml` (fallback; resolver may pick newer transitives).
   Env name: `nanoentropy`.
3. **Pip belt-and-braces** — `pip install --no-deps -r requirements.txt`
   into the env to backstop any pip-only deps.
4. **Permissions** — `chmod +x` on `install.sh`, `bin/nanoentropy`,
   `clawteam/orchestrate_review.sh`.
5. **Editable install** — `pip install --no-deps -e .` so the
   `nanoentropy` console script lands on PATH.
6. **Bash completion** — copies `bin/nanoentropy.bash-completion` to
   `~/.local/share/bash-completion/completions/nanoentropy` so tab
   completion works in subsequent shells.
7. **Doctor** — runs `nanoentropy doctor` to verify Python version, all
   pinned Python packages, modkit, samtools, `REFERENCE_DIR`,
   `RESULTS_DIR` writability, and executable bits.
8. **Selftest** — runs `nanoentropy selftest` which synthesises a tiny
   FASTA + BAM with MM/ML methylation tags and pipes it through
   modkit → entropy → bedgraphs for `true_mc` and `bisulfite` modes.

To skip the selftest (e.g. on a slow CI runner), set
`NANOENTROPY_SKIP_SELFTEST=1` before `bash install.sh`. Reference FASTAs
are NOT downloaded by default — see the *Reference data* section above
for the three ways to populate them.

## What `nanoentropy doctor` checks

| Component | Expected | Action if failing |
|---|---|---|
| `python` | `>=3.10` | re-run install.sh |
| 10 pinned Python packages | exact `==` versions | re-run install.sh |
| `modkit` | `0.6.0` | re-run install.sh; check conda env active |
| `REFERENCE_DIR` | exists | set `NANOENTROPY_REF_DIR` env var |
| `RESULTS_DIR` | writable | check filesystem permissions |
| `+x install.sh`, `+x clawteam/orchestrate_review.sh` | executable bit set | re-run install.sh |

Doctor exits `0` only when every check passes; non-zero otherwise.

## What `nanoentropy selftest` checks

Synthesises:

- a 200-bp FASTA with 8 known CpG positions
- a 30-read BAM with `MM:Z:C+m?` and `ML:B:C` methylation tags

Then runs the full CLI pipeline:

```
nanoentropy run --bam tiny.bam --fasta tiny.fa --out-dir out/<mode> ...
```

for `--mode true_mc` and `--mode bisulfite` (and `--mode ternary` if you
pass `--include-ternary`). For each mode, asserts that
`<out>/sample_<mode>.{coverage,me,mml}.bedgraph` exist and are non-empty.

Exits `0` if every mode produced expected outputs; non-zero on any
failure with the offending mode's stdout/stderr printed.

## Reference data

v4 does **not** ship the multi-GB reference FASTAs. Set
`NANOENTROPY_REF_DIR` to a directory containing:

```
hg38.fa, hg38.fa.fai
mm10.fa, mm10.fa.fai
hg38.ncbiRefSeq.gtf.gz
mm10.ncbiRefSeq.gtf.gz
```

`.env.example` ships with the path to a recommended starting reference
directory. Copy to `.env` and source before launching the GUI/CLI:

```bash
cp .env.example .env
# edit if your reference paths differ
source .env
```

The doctor flags a missing `NANOENTROPY_REF_DIR`; selftest doesn't need
it (the synthetic FASTA is self-contained).

## Reproducibility year-over-year

- `conda-lock.yml` pins every conda package to an exact build hash. The
  same lockfile installed in 2027 produces the same package set as in
  2026.
- `pyproject.toml` and `requirements.txt` use `==` for every Python dep.

If a future Python or conda package release breaks the API of a pinned
version, conda installs the pinned version anyway (it does not silently
upgrade).

## Bash completion

`install.sh` copies `bin/nanoentropy.bash-completion` to
`~/.local/share/bash-completion/completions/nanoentropy`. Open a new
shell and tab-complete on `nanoentropy <Tab>`. Completes:

- subcommands (`extract`, `entropy`, `plot`, `run`, `doctor`, `selftest`)
- option flags per subcommand
- mode enums (`true_mc`, `bisulfite`, `ternary`)
- genome enums (`hg38`, `hg19`, `mm10`)
- `--bam`, `--fasta`, `--gtf`, etc. → file completions
- `--out-dir`, `--work-dir` → directory completions

## Rolling back

`conda env remove -n nanoentropy` cleanly removes the v4 env without
touching any other env.

## Troubleshooting

```bash
nanoentropy doctor                       # something broke after install
nanoentropy selftest --include-ternary   # exercise full pipeline
streamlit run app.py 2>&1 | tee app.log  # GUI crashes
which nanoentropy                        # confirm console script on PATH
python -m src.cli --help                 # equivalent direct invocation
```
