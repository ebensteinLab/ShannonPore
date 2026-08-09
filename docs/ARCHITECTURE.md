# Architecture — shannonpore v4

## High level

```
┌──────────────────────────────────────────────────────────────────┐
│                  app.py (Streamlit entrypoint)                   │
│        page config · sidebar · tab dispatch · style              │
└──────────────────────────────────────────────────────────────────┘
                                │
        ┌───────────────────────┼───────────────────────────┐
        ▼                       ▼                           ▼
┌──────────────────┐   ┌──────────────────┐    ┌─────────────────────┐
│ shannonpore/tabs/        │   │ shannonpore/state.py     │    │ shannonpore/config.py       │
│ ├ tab_file_prep  │   │ AppState (DC)    │    │ env-resolved paths      │
│ └ tab_graph_prep │   │ FilePrepState    │    │ SHANNONPORE_REF_DIR     │
└──────────────────┘   │ GraphPrepState   │    │ SHANNONPORE_RESULTS_DIR │
        │              └──────────────────┘    │ HG38 / MM10 assets      │
        │                                      └─────────────────────────┘
        │
        ├─► shannonpore/pipelines/  (orchestrator · modkit_runner · bam_utils ·
        │                    whole_genome_duckdb · ternary · roi_entropy)
        ├─► shannonpore/plots/      (tracks · scatter [ME/MML + arch + paired] ·
        │                    export [multi-format save] · theme)
        ├─► shannonpore/io/         (bedgraph · gtf_utils · roi_utils · utils_io)
        └─► shannonpore/ui/         (error_handler · progress · widgets · style)
```

`shannonpore/cli.py` exposes the same pipelines as a Bash-friendly argparse CLI
with subcommands `extract`, `entropy`, `plot`, `run`, `doctor`, `selftest`,
plus `guide` and `examples` for inline documentation. The CLI shares 100%
of its compute path with the GUI — both call into `shannonpore/pipelines/*`.

## File budget

| Layer | Lines | Notes |
|---|---:|---|
| `app.py` | <200 | page config, sidebar, tab dispatch, style injection |
| `shannonpore/cli.py` | ~950 | argparse + 8 subcommands; the largest module |
| `shannonpore/tabs/tab_file_prep.py` | ~450 | input + entropy-mode + bins + output |
| `shannonpore/tabs/tab_graph_prep.py` | ~600 | samples · ME/MML scatter · arch · paired · region track |
| `shannonpore/state.py` | ~150 | AppState dataclass + accessors |
| `shannonpore/pipelines/*` | 6 modules | orchestrator + extract/merge + 3 entropy strategies |
| `shannonpore/plots/*` | 4 modules | scatter (3 plot families) · tracks · export (png/jpg/svg/pdf + DPI) · theme |
| `shannonpore/ui/style.py` | ~260 | CSS injection + header band |

Plot logic is decoupled from Streamlit; every plot module is importable
without a Streamlit context. The current outlier on size is `cli.py`,
which is mostly argparse boilerplate — splitting it is a tracked nice-to-have.

## Session state

`shannonpore/state.py` defines `AppState` as a frozen-equivalent (mutated via
`dataclasses.replace`) dataclass. The Streamlit `session_state` holds a
single key, `_app_state`, pointing to the current instance. Reads via
`get_state()`, writes via `update_section('file_prep'|'graph_prep', …)`.

See [SESSION_STATE_SCHEMA.md](SESSION_STATE_SCHEMA.md) for the full schema.

## Visual identity

The GUI ships its own theme via:

- `.streamlit/config.toml` — base palette (light, prussian primary).
- `shannonpore/ui/style.py` — IBM Plex font stack, hairline rules, flat
  buttons, monospaced metrics, and a `shannonpore · v4 · GUI` band at
  the top. Hides Streamlit's MainMenu / footer chrome.

Replaces the generic Streamlit prototype look with something that reads
as laboratory software: restrained, high-contrast, monospaced where
data is shown.

## Error handling

`shannonpore/ui/error_handler.show_error` decorates every Streamlit callback
that touches the filesystem or external binaries. It logs the full
traceback (stdlib `logging`) and renders a friendly `st.error` plus a
collapsible traceback expander.

## Entropy modes (Tab 1 default = True-mC)

| Mode | mod_code filter | Math behaviour |
|---|---|---|
| `true_mc` (default) | `WHERE mod_code = 'm'` | 5hmC ignored; only true 5mC drives entropy. |
| `bisulfite` | `WHERE mod_code IN ('m','h')` + SUM | 5hmC counted as 5mC. |
| `ternary` | pivots m/h, assigns state ∈ {0,1,2} | 3-state entropy normalised by k·log2(3). Requires ≥ 3^k coverage. Emits extra `mhml` bedgraph. |

Ternary uses `shannonpore/pipelines/ternary_entropy.py`; the 2-state modes share
`shannonpore/pipelines/whole_genome_duckdb_pipeline.py`.

## Reference data

Reference FASTAs / GTFs are NOT bundled. Path resolved from
`SHANNONPORE_REF_DIR` (default `v4/reference_files/`), set via `.env` or
the user's shell. v4 ships `.env.example` for continuity.

## Reproducibility chain

```
bash install.sh
   ├── conda-lock.yml     (Python · modkit · system libs at exact versions)
   ├── requirements.txt   (== pinned pip deps)
   ├── pip install -e .   (shannonpore console script on PATH)
   ├── shannonpore doctor (every pinned version, every permission)
   └── shannonpore selftest (synthetic-BAM end-to-end pipeline)
```

If a 2027 user runs `bash install.sh`, they get the 2026 stack.
