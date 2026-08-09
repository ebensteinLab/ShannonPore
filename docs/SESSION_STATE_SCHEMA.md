# Session State Schema

`shannonpore/state.py` defines the complete session state for the v4 Streamlit
app. Every key is typed; reads/writes go through `get_state()` /
`update_section()`. Streamlit's `session_state` dict holds one entry,
`_app_state`, that points to an `AppState` instance.

## Top level — `AppState`

| Field | Type | Default | Notes |
|---|---|---|---|
| `file_prep` | `FilePrepState` | dataclass default | Tab 1 |
| `graph_prep` | `GraphPrepState` | dataclass default | Tab 2 |
| `sidebar_collapsed` | `bool` | `False` | UI |
| `last_error` | `str \| None` | `None` | last unhandled exception |

## `FilePrepState` (Tab 1)

| Field | Type | Default |
|---|---|---|
| `target` | `SampleSpec` | `SampleSpec(label="target")` |
| `control` | `SampleSpec` | `SampleSpec(label="control")` |
| `pair_mode` | `bool` | `False` |
| `genome` | `"hg38" \| "mm10"` | `"hg38"` |
| `custom_fasta` | `Path \| None` | `None` |
| `roi_text` | `str` | `""` |
| `roi_bed` | `Path \| None` | `None` |
| `use_whole_genome` | `bool` | `True` |
| `modkit_threads` | `int` | `8` |
| `modkit_queue_size` | `int` | `1000` |
| `modkit_extra_args` | `str` | `""` |
| `include_bed` | `Path \| None` | `None` |
| `output_dir` | `Path \| None` | `None` |
| `cpgs_per_bin` | `int` | `4` |
| `min_coverage` | `int` | `16` |
| `methyl_threshold` | `float` | `0.5` |
| `entropy_mode` | `"true_mc" \| "bisulfite" \| "ternary"` | `"true_mc"` |
| `last_run_id` | `str \| None` | `None` |
| `last_run_summary` | `dict[str, Any]` | `{}` |

### `SampleSpec`

| Field | Type | Default |
|---|---|---|
| `label` | `str` | `"sample"` (sanitized to `[A-Za-z0-9._-]` before becoming a path component) |
| `input_kind` | `"bam" \| "tsv" \| "bam_folder"` | `"bam"` |
| `bam_path` | `Path \| None` | `None` |
| `tsv_path` | `Path \| None` | `None` |
| `bam_folder` | `Path \| None` | `None` |

## `GraphPrepState` (Tab 2)

| Field | Type | Default |
|---|---|---|
| `control` | `TrackPlot` | name="Control", color="#2980b9" |
| `target` | `TrackPlot` | name="Target", color="#e67e22" |
| `genome` | `"hg38" \| "mm10"` | `"hg38"` |
| `gtf_path` | `Path \| None` | `None` (custom override; `None` ⇒ use bundled GTF for `genome`) |
| `region_chrom` | `str` | `""` |
| `region_start` | `int` | `0` |
| `region_end` | `int` | `0` |
| `plots_dir` | `Path \| None` | `None` |
| `scatter_subsample` | `int` | `50_000` (legacy; current plots don't use it) |
| `scatter_color_by_density` | `bool` | `True` (legacy; current plots don't use it) |

### `TrackPlot`

| Field | Type | Default |
|---|---|---|
| `name` | `str` | `""` |
| `color` | `str` | `"#1f77b4"` (overridden per side by `GraphPrepState` factories) |
| `mml_path` | `Path \| None` | `None` |
| `me_path` | `Path \| None` | `None` |
| `coverage_path` | `Path \| None` | `None` |

## Invariants

- All `Path` fields hold absolute or `~`-expanded paths.
- `entropy_mode` is one of the values in `shannonpore.constants.ENTROPY_MODES`.
- `cpgs_per_bin >= 2`.
- `region_end >= region_start`.

These are enforced at widget level (Streamlit `min_value` / radio
options); pipeline functions also re-validate.
