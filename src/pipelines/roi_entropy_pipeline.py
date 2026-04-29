"""ROI-based methylation entropy pipeline.

Streams modkit-extracted CpG calls (TSV), bins per region, computes per-bin
methylation level (MML) and per-bin Shannon entropy normalised by the
number of CpGs in the bin (ME).

v4 fix: in v3 `build_parquet_from_tsv_stream_duckdb` referenced undefined
`force` and `stream_cb` parameters. The signature now exposes both.
"""

from __future__ import annotations

import logging
import os
import re
from collections.abc import Callable, Iterator
from typing import Any

import numpy as np
import pandas as pd
from pyfaidx import Fasta

from src.io.utils_io import safe_mkdir

logger = logging.getLogger(__name__)


# ─── TSV streaming ────────────────────────────────────────────────────────

_MODKIT_KEEP_COLS = [
    "read_id", "forward_read_position", "ref_position",
    "chrom", "ref_strand", "mod_qual",
]


def stream_grouped_modkit_table_tsv(
    tsv_path: str,
    chunk_size: int,
) -> Iterator[pd.DataFrame]:
    compression = "gzip" if tsv_path.lower().endswith(".gz") else None
    reader = pd.read_csv(
        tsv_path,
        sep="\t",
        usecols=_MODKIT_KEEP_COLS,
        chunksize=chunk_size,
        iterator=True,
        comment="#",
        compression=compression,
    )
    for chunk in reader:
        grouped = chunk.groupby(
            ["read_id", "forward_read_position", "ref_position", "chrom", "ref_strand"],
            sort=False,
            as_index=False,
        )["mod_qual"].sum()
        yield grouped


def can_use_duckdb() -> bool:
    try:
        import duckdb  # noqa: F401
    except ImportError:
        return False
    return True


def build_parquet_from_tsv_stream_duckdb(
    input_path: str,
    parquet_out_path: str,
    chrom_whitelist: list[str],
    chunk_size: int = 1_000_000,
    force: bool = False,
    stream_cb: Callable[[str], None] | None = None,
) -> str:
    """Stream a modkit TSV, group rows by (read_id, ref_position, chrom,
    strand), and write a parquet cache. Idempotent: skips when the parquet
    already exists unless `force=True`.
    """
    import duckdb

    safe_mkdir(os.path.dirname(parquet_out_path))

    if (
        os.path.exists(parquet_out_path)
        and os.path.getsize(parquet_out_path) > 0
        and not force
    ):
        msg = f"[INFO] Parquet already exists, skipping TSV import: {parquet_out_path}"
        if stream_cb:
            stream_cb(msg)
        else:
            logger.info(msg)
        return parquet_out_path
    if os.path.exists(parquet_out_path) and force:
        os.remove(parquet_out_path)

    con = duckdb.connect(database=":memory:")
    con.execute(
        """
        CREATE TABLE t (
            read_id VARCHAR,
            forward_read_position BIGINT,
            ref_position BIGINT,
            chrom VARCHAR,
            ref_strand VARCHAR,
            mod_qual DOUBLE
        );
        """
    )

    chrom_set = set(chrom_whitelist)
    for g in stream_grouped_modkit_table_tsv(input_path, chunk_size=chunk_size):
        if not g.empty:
            g = g[g["chrom"].astype(str).isin(chrom_set)]
            if not g.empty:
                con.register("chunk_df", g)
                con.execute("INSERT INTO t SELECT * FROM chunk_df;")
                con.unregister("chunk_df")

    con.execute(f"COPY t TO '{parquet_out_path}' (FORMAT PARQUET);")
    con.close()

    if not os.path.exists(parquet_out_path) or os.path.getsize(parquet_out_path) == 0:
        raise RuntimeError("DuckDB parquet export failed or produced empty file.")
    return parquet_out_path


def load_grouped_table_to_df_from_tsv(tsv_path: str, chunk_size: int) -> pd.DataFrame:
    chunks = []
    for g in stream_grouped_modkit_table_tsv(tsv_path, chunk_size=chunk_size):
        chunks.append(g)
    if not chunks:
        return pd.DataFrame(columns=_MODKIT_KEEP_COLS)
    df = pd.concat(chunks, ignore_index=True)
    df.loc[df["ref_strand"] == "-", "ref_position"] -= 1
    return df


# ─── Region/CpG handling ──────────────────────────────────────────────────

def find_cpg_locations_from_atlas(
    fasta_path: str, regions_df: pd.DataFrame
) -> pd.DataFrame:
    df = regions_df.copy()
    if "target" not in df.columns:
        df["target"] = df["chr"]
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    genome = Fasta(fasta_path)

    def extract_sequence(row: pd.Series) -> str | None:
        try:
            return genome[row["chr"]][row["start"]: row["end"]].seq
        except (KeyError, ValueError) as exc:
            logger.warning("FASTA lookup failed for %s:%s-%s: %s",
                           row["chr"], row["start"], row["end"], exc)
            return None

    df["sequence"] = df.apply(extract_sequence, axis=1)

    def extract_cpg(row: pd.Series) -> list[int]:
        seq = row["sequence"]
        if seq is None or pd.isna(seq):
            return []
        return [row["start"] + m.start() for m in re.finditer(r"CG", seq, re.IGNORECASE)]

    df["cpgs"] = df.apply(extract_cpg, axis=1)
    df = df.drop(columns=["sequence"])
    return df


def bin_cpg_positions_per_region(
    df_regions: pd.DataFrame, cpg_per_bin_count: int
) -> pd.DataFrame:
    results = []
    for _, region in df_regions.iterrows():
        region_chr = region["chr"]
        r_start, r_end = int(region["start"]), int(region["end"])
        target = region.get("target", region_chr)
        cpg_pos = sorted(region["cpgs"])

        bin_edges: list[tuple[int, int]] = []
        cpg_per_bin: list[list[int]] = []

        if cpg_pos and cpg_per_bin_count > 0:
            num_bins = int(np.ceil(len(cpg_pos) / cpg_per_bin_count))
            bins_array = np.array_split(np.array(cpg_pos), num_bins)
            for i, bin_arr in enumerate(bins_array):
                if bin_arr.size == 0:
                    continue
                start_edge = r_start if i == 0 else bin_edges[-1][1] + 1
                bin_edges.append((start_edge, int(bin_arr[-1])))
                cpg_per_bin.append(bin_arr.tolist())
            if bin_edges:
                bin_edges[-1] = (bin_edges[-1][0], r_end)

        results.append({
            "chr": region_chr,
            "start": r_start,
            "end": r_end,
            "target": target,
            "cpg_per_bin": cpg_per_bin,
            "bin_edges": bin_edges,
        })
    return pd.DataFrame(results)


# ─── Entropy math ─────────────────────────────────────────────────────────

def vect_to_num(bits: np.ndarray) -> int:
    return int(sum(2 ** (len(bits) - 1 - i) * int(bits[i]) for i in range(len(bits))))


def entropy_vec(values: np.ndarray, config_options_num: int) -> float:
    """Shannon entropy over methylation patterns, normalised by log2(2^k)=k."""
    if values.size == 0 or config_options_num <= 0:
        return 0.0
    _unique, counts = np.unique(values, return_counts=True)
    total = int(np.sum(counts))
    if total == 0:
        return 0.0
    probs = counts / total
    probs = probs[probs > 0]
    return float((-np.sum(probs * np.log2(probs))) / config_options_num)


def process_bin(
    bin_events: pd.DataFrame, cpg_rel: np.ndarray, methyl_score_thresh: float
) -> dict[str, np.ndarray]:
    if bin_events.empty or cpg_rel.size == 0:
        return {}
    event_pos = bin_events["rel_pos"].values
    match_matrix = (event_pos[:, None] == cpg_rel[None, :])
    if not match_matrix.any():
        return {}

    matched_rows = np.where(match_matrix.any(axis=1))[0]
    matched_cpg_idx = np.argmax(match_matrix[matched_rows], axis=1)

    matched = bin_events.iloc[matched_rows].copy()
    matched["cpg_idx"] = matched_cpg_idx
    matched["status"] = (matched["mod_qual"] > methyl_score_thresh).astype(int)

    pivot = matched.pivot_table(
        index="read_id", columns="cpg_idx", values="status"
    ).dropna()
    return {rid: row.values.astype(int) for rid, row in pivot.iterrows()}


def merge_asof_by_chrom(
    df_left: pd.DataFrame,
    df_right: pd.DataFrame,
    left_on: str,
    right_on: str,
    by: str,
) -> pd.DataFrame:
    merged_groups = []
    right_groups = set(df_right[by].astype(str).unique().tolist())
    for group, left_group in df_left.groupby(by):
        if str(group) not in right_groups:
            continue
        right_group = df_right[df_right[by] == group].sort_values(right_on)
        merged = pd.merge_asof(
            left_group.sort_values(left_on),
            right_group,
            left_on=left_on,
            right_on=right_on,
            direction="backward",
        )
        merged_groups.append(merged)
    if not merged_groups:
        return pd.DataFrame()
    return pd.concat(merged_groups).reset_index(drop=True)


def fill_cpg_methylation_optimized_df(
    df_grouped: pd.DataFrame,
    df_cpgs_binned: pd.DataFrame,
    methyl_score_thresh: float,
) -> pd.DataFrame:
    df_regions = df_cpgs_binned.sort_values(["chr", "start"]).reset_index(drop=True)
    df_regions["region_id"] = df_regions.index
    df_regions = df_regions.rename(columns={"chr": "chrom"})

    merged = merge_asof_by_chrom(
        df_grouped.sort_values(["chrom", "ref_position"]),
        df_regions,
        left_on="ref_position",
        right_on="start",
        by="chrom",
    )

    if merged.empty:
        df_regions["meth_data"] = [[] for _ in range(len(df_regions))]
        return df_regions

    merged = merged[merged["ref_position"] <= merged["end"]].copy()
    grouped = merged.groupby("region_id")

    def process_region(region_id: int) -> list[dict]:
        region = df_regions.loc[region_id]
        if region_id not in grouped.groups:
            return [{} for _ in region["bin_edges"]]

        r_start = int(region["start"])
        bin_edges = region["bin_edges"]
        cpg_per_bin = region["cpg_per_bin"]

        events = grouped.get_group(region_id).copy()
        events["rel_pos"] = events["ref_position"] - r_start
        rel_cpg_per_bin = [np.array(cpg_list) - r_start for cpg_list in cpg_per_bin]

        meth_data: list[dict] = [{} for _ in range(len(bin_edges))]
        for i in range(len(bin_edges)):
            b_start = bin_edges[i][0] - r_start
            b_end = bin_edges[i][1] - r_start
            mask = (events["rel_pos"] >= b_start) & (events["rel_pos"] <= b_end)
            meth_data[i] = process_bin(
                events.loc[mask], rel_cpg_per_bin[i], methyl_score_thresh
            )
        return meth_data

    df_regions["meth_data"] = df_regions["region_id"].apply(process_region)
    return df_regions


def calculate_methylation_entropy_df(df_result: pd.DataFrame) -> pd.DataFrame:
    def calc_one(meth_data: list[dict]) -> list[list[float]]:
        num_bins = len(meth_data)
        mml = np.full(num_bins, -1.0)
        ent = np.full(num_bins, -1.0)
        var = np.full(num_bins, -1.0)
        cov = np.zeros(num_bins, dtype=int)
        for i, bin_data in enumerate(meth_data):
            if bin_data:
                cov[i] = len(bin_data)
                all_meth = np.array(list(bin_data.values()))
                if all_meth.size > 0:
                    mml[i] = float(np.mean(all_meth))
                    if all_meth.shape[0] > 1:
                        ent[i] = entropy_vec(
                            np.apply_along_axis(vect_to_num, 1, all_meth),
                            all_meth.shape[1],
                        )
                        var[i] = float(np.var(all_meth))
        return [mml.tolist(), ent.tolist(), var.tolist(), cov.tolist()]

    df_result[
        ["mml_by_bin", "entropy_by_bin", "variance_by_bin", "coverage_by_bin"]
    ] = pd.DataFrame(
        df_result["meth_data"].apply(calc_one).tolist(),
        index=df_result.index,
    )
    return df_result


def unroll_for_plotting(
    df_result: pd.DataFrame, coverage_threshold: int
) -> pd.DataFrame:
    rows = []
    for _, region in df_result.iterrows():
        chrom = region["chrom"]
        target_label = region.get(
            "target", f"{chrom}:{region['start']}-{region['end']}"
        )
        for edges, mml, ent, cov in zip(
            region["bin_edges"],
            region["mml_by_bin"],
            region["entropy_by_bin"],
            region["coverage_by_bin"],
        ):
            if int(cov) >= int(coverage_threshold):
                rows.append({
                    "chrom": chrom,
                    "start": int(edges[0]),
                    "end": int(edges[1]),
                    "mml": float(mml),
                    "entropy": float(ent),
                    "coverage": int(cov),
                    "region_target": target_label,
                })
    if not rows:
        return pd.DataFrame(
            columns=[
                "chrom", "start", "end", "mml", "entropy", "coverage", "region_target",
            ]
        )
    return pd.DataFrame(rows)


def create_bedgraph_from_df(df: pd.DataFrame, value_col: str, out_path: str) -> None:
    safe_mkdir(os.path.dirname(out_path))
    if df.empty:
        open(out_path, "w").close()
        return
    df2 = df[["chrom", "start", "end", value_col]].copy()
    df2.to_csv(out_path, sep="\t", header=False, index=False)


# ─── Top-level orchestration ──────────────────────────────────────────────

def run_roi_pipeline(
    control_tsv_path: str,
    target_tsv_path: str,
    reference_fasta: str,
    regions_df: pd.DataFrame,
    cpgs_per_bin: int,
    coverage_threshold: int,
    ingest_chunk_rows: int,
    out_dir: str,
    methyl_score_thresh: float = 0.5,
    prefer_duckdb: bool = True,
    progress_cb: Callable[[str], None] | None = None,
) -> dict[str, Any]:
    tables_dir = safe_mkdir(os.path.join(out_dir, "tables"))
    bedgraph_dir = safe_mkdir(os.path.join(out_dir, "bedgraphs"))
    parquet_dir = safe_mkdir(os.path.join(out_dir, "grouped_parquet_cache"))

    if progress_cb:
        progress_cb("ROI: finding CpGs in regions...")
    df_regions_cpg = find_cpg_locations_from_atlas(reference_fasta, regions_df)
    df_binned = bin_cpg_positions_per_region(df_regions_cpg, int(cpgs_per_bin))
    chrom_whitelist = sorted(set(df_binned["chr"].astype(str).tolist()))

    use_duck = bool(prefer_duckdb) and can_use_duckdb()
    if use_duck:
        if progress_cb:
            progress_cb("ROI: DuckDB available, building grouped parquet caches...")
        target_parquet = os.path.join(parquet_dir, "target_grouped.parquet")
        control_parquet = os.path.join(parquet_dir, "control_grouped.parquet")

        build_parquet_from_tsv_stream_duckdb(
            input_path=target_tsv_path,
            parquet_out_path=target_parquet,
            chrom_whitelist=chrom_whitelist,
            chunk_size=int(ingest_chunk_rows),
            stream_cb=progress_cb,
        )
        build_parquet_from_tsv_stream_duckdb(
            input_path=control_tsv_path,
            parquet_out_path=control_parquet,
            chrom_whitelist=chrom_whitelist,
            chunk_size=int(ingest_chunk_rows),
            stream_cb=progress_cb,
        )

        df_target_raw = pd.read_parquet(target_parquet)
        df_control_raw = pd.read_parquet(control_parquet)

        df_target_raw.loc[df_target_raw["ref_strand"] == "-", "ref_position"] -= 1
        df_control_raw.loc[df_control_raw["ref_strand"] == "-", "ref_position"] -= 1
    else:
        if progress_cb:
            progress_cb("ROI: DuckDB not available, falling back to pandas chunk streaming...")
        df_target_raw = load_grouped_table_to_df_from_tsv(
            target_tsv_path, chunk_size=int(ingest_chunk_rows)
        )
        df_control_raw = load_grouped_table_to_df_from_tsv(
            control_tsv_path, chunk_size=int(ingest_chunk_rows)
        )

    if progress_cb:
        progress_cb("ROI: filling CpG methylation per region for TARGET...")
    df_target_regions = fill_cpg_methylation_optimized_df(
        df_target_raw, df_binned, float(methyl_score_thresh)
    )
    df_target_regions = calculate_methylation_entropy_df(df_target_regions)

    if progress_cb:
        progress_cb("ROI: filling CpG methylation per region for CONTROL...")
    df_control_regions = fill_cpg_methylation_optimized_df(
        df_control_raw, df_binned, float(methyl_score_thresh)
    )
    df_control_regions = calculate_methylation_entropy_df(df_control_regions)

    df_target_regions.to_pickle(os.path.join(tables_dir, "target_regions.pkl"))
    df_control_regions.to_pickle(os.path.join(tables_dir, "control_regions.pkl"))

    plot_target = unroll_for_plotting(df_target_regions, int(coverage_threshold))
    plot_control = unroll_for_plotting(df_control_regions, int(coverage_threshold))

    paths = {
        "control_mml": os.path.join(bedgraph_dir, "control_mml.bedgraph"),
        "target_mml": os.path.join(bedgraph_dir, "target_mml.bedgraph"),
        "control_me": os.path.join(bedgraph_dir, "control_me.bedgraph"),
        "target_me": os.path.join(bedgraph_dir, "target_me.bedgraph"),
    }
    create_bedgraph_from_df(plot_control, "mml", paths["control_mml"])
    create_bedgraph_from_df(plot_target, "mml", paths["target_mml"])
    create_bedgraph_from_df(plot_control, "entropy", paths["control_me"])
    create_bedgraph_from_df(plot_target, "entropy", paths["target_me"])

    return {
        "df_control_regions": df_control_regions,
        "df_target_regions": df_target_regions,
        "plot_control": plot_control,
        "plot_target": plot_target,
        "control_bedgraphs": {"mml": paths["control_mml"], "me": paths["control_me"]},
        "target_bedgraphs": {"mml": paths["target_mml"], "me": paths["target_me"]},
    }
