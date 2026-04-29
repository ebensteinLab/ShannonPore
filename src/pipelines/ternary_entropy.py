"""Whole-genome TERNARY methylation entropy (0=C, 1=5mC, 2=5hmC).

Adapted from `entropy/5hmc_confound_analysis/analysis/ternary_entropy.py`
to be importable (callable from Streamlit) rather than CLI-only. The
mathematics and DuckDB SQL match the original; only the I/O orchestration
is wrapped in a function instead of `main()`.

Coverage warning: ternary entropy partitions reads into 3^k pattern bins,
so per-bin coverage requirements grow as 3^k:
    k=3 CpGs/bin → ~27× coverage minimum
    k=4 CpGs/bin → ~81× coverage minimum
The Streamlit Tab 1 surfaces this warning to the user.
"""

from __future__ import annotations

import logging
import multiprocessing as mp
import os
import re
from collections.abc import Callable

import duckdb
import numpy as np
from pyfaidx import Fasta

logger = logging.getLogger(__name__)

_BASE3_CACHE: dict[int, np.ndarray] = {}


def _ensure_dir(path: str) -> None:
    if path:
        os.makedirs(path, exist_ok=True)


def _list_fasta_chroms(fasta_path: str) -> list[str]:
    return list(Fasta(fasta_path).keys())


def _factorize(values: np.ndarray) -> tuple[np.ndarray, int]:
    uniques, inv = np.unique(values, return_inverse=True)
    return inv.astype(np.int64, copy=False), int(uniques.shape[0])


def _find_cpg_positions(
    fasta_path: str, chrom: str, chunk_size: int = 10_000_000
) -> np.ndarray:
    genome = Fasta(fasta_path)
    chrom_len = len(genome[chrom])
    positions: list[int] = []
    carry = ""
    pos = 0
    while pos < chrom_len:
        end = min(pos + chunk_size, chrom_len)
        seq = carry + genome[chrom][pos:end].seq.upper()
        offset = pos - len(carry)
        for m in re.finditer("CG", seq):
            positions.append(offset + m.start())
        carry = seq[-1:] if seq else ""
        pos = end
    return np.array(positions, dtype=np.int64) if positions else np.array([], dtype=np.int64)


def _shannon(counts: np.ndarray) -> float:
    total = int(counts.sum())
    if total <= 0:
        return 0.0
    p = counts.astype(np.float64) / float(total)
    p = p[p > 0]
    return float(-np.sum(p * np.log2(p)))


def _base3_powers(k: int) -> np.ndarray:
    if k not in _BASE3_CACHE:
        _BASE3_CACHE[k] = np.array(
            [3 ** (k - 1 - i) for i in range(k)], dtype=np.int64
        )
    return _BASE3_CACHE[k]


def stage_a_build_duckdb_ternary(
    tsv_path: str,
    db_path: str,
    table_name: str,
    tmp_dir: str,
    threads: int,
    methyl_thresh: float,
    force: bool,
) -> None:
    """Stage A: ingest modkit TSV, pivot m/h per CpG, assign ternary state.

    State assignment:
        h_qual > thresh AND h_qual >= m_qual → 2 (5hmC)
        m_qual > thresh                      → 1 (5mC)
        else                                 → 0 (unmodified)
    """
    _ensure_dir(os.path.dirname(db_path))
    _ensure_dir(tmp_dir)

    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"TSV not found: {tsv_path}")
    with open(tsv_path, "rb") as fh:
        head = fh.read(4096)
    if b"\t" not in head:
        raise ValueError("Input TSV does not look tab-delimited.")

    con = duckdb.connect(db_path)
    con.execute(f"PRAGMA temp_directory='{tmp_dir}';")
    con.execute(f"PRAGMA threads={int(threads)};")

    exists = con.execute(
        "SELECT COUNT(*) FROM information_schema.tables WHERE table_name = ?",
        [table_name],
    ).fetchone()[0] > 0
    if exists and not force:
        con.close()
        return
    if exists and force:
        con.execute(f"DROP TABLE {table_name};")

    thresh = float(methyl_thresh)
    con.execute(
        f"""
        CREATE TABLE {table_name} AS
        WITH raw AS (
            SELECT
                read_id::VARCHAR  AS read_id,
                chrom::VARCHAR    AS chrom,
                CASE
                    WHEN ref_strand = '-' THEN (ref_position::BIGINT - 1)
                    ELSE ref_position::BIGINT
                END               AS ref_position,
                mod_qual::DOUBLE  AS mod_qual,
                mod_code::VARCHAR AS mod_code
            FROM read_csv_auto(
                '{tsv_path}',
                delim='\t',
                header=true,
                ignore_errors=true,
                null_padding=true
            )
            WHERE read_id IS NOT NULL
              AND chrom IS NOT NULL
              AND ref_position IS NOT NULL
              AND mod_qual IS NOT NULL
              AND ref_strand IS NOT NULL
              AND mod_code IN ('m', 'h')
        ),
        pivoted AS (
            SELECT
                read_id, chrom, ref_position,
                MAX(CASE WHEN mod_code = 'm' THEN mod_qual ELSE 0.0 END) AS m_qual,
                MAX(CASE WHEN mod_code = 'h' THEN mod_qual ELSE 0.0 END) AS h_qual
            FROM raw
            GROUP BY read_id, chrom, ref_position
        )
        SELECT
            read_id, chrom, ref_position, m_qual, h_qual,
            CASE
                WHEN h_qual > {thresh} AND h_qual >= m_qual THEN 2
                WHEN m_qual > {thresh} THEN 1
                ELSE 0
            END AS state
        FROM pivoted;
        """
    )
    con.execute(f"CREATE INDEX idx_{table_name}_chrom_pos ON {table_name}(chrom, ref_position);")
    con.execute(f"CREATE INDEX idx_{table_name}_chrom ON {table_name}(chrom);")
    con.close()


def compute_chrom_ternary(
    chrom: str,
    fasta_path: str,
    db_path: str,
    table_name: str,
    out_chrom_dir: str,
    cpgs_per_bin: int,
    min_coverage: int,
    fasta_chunk: int,
) -> str:
    """Worker: per-chrom ternary ME / MML / MhML / coverage / cov_5mc / cov_5hmc."""
    _ensure_dir(out_chrom_dir)
    suffixes = ["coverage", "me", "mml", "mhml", "cov_5mc", "cov_5hmc"]
    out_paths = {s: os.path.join(out_chrom_dir, f"{chrom}.{s}.bedgraph") for s in suffixes}
    if all(os.path.exists(p) for p in out_paths.values()):
        return f"[INFO] {chrom}: outputs exist, skipping"

    k = int(cpgs_per_bin)
    cpg_positions = _find_cpg_positions(fasta_path, chrom, chunk_size=int(fasta_chunk))
    if cpg_positions.size < k:
        for p in out_paths.values():
            open(p, "w").close()
        return f"[INFO] {chrom}: not enough CpGs ({cpg_positions.size})"

    n_full = int(cpg_positions.size // k)
    starts = cpg_positions[0: n_full * k: k]
    ends = cpg_positions[(k - 1): n_full * k: k] + 2

    con = duckdb.connect(db_path, read_only=True)
    rows = con.execute(
        f"SELECT read_id, ref_position, state FROM {table_name} WHERE chrom = ?",
        [chrom],
    ).fetchall()
    con.close()

    def _empty_outputs(reason: str) -> str:
        with open(out_paths["coverage"], "w") as f:
            for i in range(n_full):
                f.write(f"{chrom}\t{int(starts[i])}\t{int(ends[i])}\t0\n")
        for s in suffixes[1:]:
            open(out_paths[s], "w").close()
        return f"[INFO] {chrom}: {reason}"

    if not rows:
        return _empty_outputs("no events")

    read_ids = np.array([r[0] for r in rows], dtype=object)
    pos = np.array([r[1] for r in rows], dtype=np.int64)
    state = np.array([r[2] for r in rows], dtype=np.int8)
    read_code, _ = _factorize(read_ids)

    idx = np.searchsorted(cpg_positions, pos)
    valid = (idx >= 0) & (idx < cpg_positions.size)
    idx2 = idx[valid]
    matches = cpg_positions[idx2] == pos[valid]
    if not np.any(matches):
        return _empty_outputs("no CpG-matching events")

    ordinal = idx2[matches].astype(np.int64)
    state = state[valid][matches]
    read_code = read_code[valid][matches].astype(np.int64)
    bin_id = (ordinal // k).astype(np.int64)
    in_full = bin_id < n_full
    if not np.any(in_full):
        return _empty_outputs("events but no full bins")

    ordinal = ordinal[in_full]
    bin_id = bin_id[in_full]
    state = state[in_full]
    read_code = read_code[in_full]
    cpg_idx = (ordinal % k).astype(np.int8)

    keys = np.rec.fromarrays([read_code, bin_id, cpg_idx], names="r,b,c")
    order = np.argsort(keys, kind="mergesort")
    keys = keys[order]
    state = state[order]

    uniq_mask = np.ones(keys.shape[0], dtype=bool)
    uniq_mask[1:] = (keys[1:] != keys[:-1])
    uniq_idx = np.where(uniq_mask)[0]
    group_ends = np.r_[uniq_idx[1:], keys.shape[0]]

    dedup_read = read_code[order][uniq_idx]
    dedup_bin = bin_id[order][uniq_idx]
    dedup_cpg_idx = cpg_idx[order][uniq_idx]
    dedup_state = np.zeros(uniq_idx.shape[0], dtype=np.int8)
    for i in range(uniq_idx.shape[0]):
        dedup_state[i] = state[uniq_idx[i]:group_ends[i]].max()

    powers = _base3_powers(k)
    weighted = (dedup_state.astype(np.int64) * powers[dedup_cpg_idx]).astype(np.int64)
    obs_bit = (1 << dedup_cpg_idx.astype(np.int64)).astype(np.int64)
    is_5mc = (dedup_state == 1).astype(np.int64)
    is_5hmc = (dedup_state == 2).astype(np.int64)

    rb_keys = np.rec.fromarrays([dedup_read, dedup_bin], names="r,b")
    rb_order = np.argsort(rb_keys, kind="mergesort")
    rb_keys = rb_keys[rb_order]
    weighted = weighted[rb_order]
    obs_bit = obs_bit[rb_order]
    is_5mc = is_5mc[rb_order]
    is_5hmc = is_5hmc[rb_order]

    rb_uniq = np.ones(rb_keys.shape[0], dtype=bool)
    rb_uniq[1:] = (rb_keys[1:] != rb_keys[:-1])
    rb_idx = np.where(rb_uniq)[0]
    rb_end = np.r_[rb_idx[1:], rb_keys.shape[0]]

    rb_bin = dedup_bin[rb_order][rb_idx]
    n_rb = rb_idx.shape[0]
    pattern = np.zeros(n_rb, dtype=np.int64)
    obsmask = np.zeros(n_rb, dtype=np.int64)
    mc_count = np.zeros(n_rb, dtype=np.int64)
    hmc_count = np.zeros(n_rb, dtype=np.int64)
    for i in range(n_rb):
        s0, s1 = rb_idx[i], rb_end[i]
        pattern[i] = weighted[s0:s1].sum()
        obsmask[i] = obs_bit[s0:s1].sum()
        mc_count[i] = is_5mc[s0:s1].sum()
        hmc_count[i] = is_5hmc[s0:s1].sum()

    full_mask_val = (1 << k) - 1
    full = (obsmask == full_mask_val)
    if not np.any(full):
        return _empty_outputs(f"no reads cover all {k} CpGs in any bin")

    rb_bin = rb_bin[full]
    pattern = pattern[full]
    mc_count = mc_count[full]
    hmc_count = hmc_count[full]

    pair = np.rec.fromarrays([rb_bin, pattern], names="b,p")
    pair_order = np.argsort(pair, kind="mergesort")
    rb_bin = rb_bin[pair_order]
    pattern = pattern[pair_order]
    mc_count = mc_count[pair_order]
    hmc_count = hmc_count[pair_order]

    uniq = np.ones(rb_bin.shape[0], dtype=bool)
    uniq[1:] = (rb_bin[1:] != rb_bin[:-1]) | (pattern[1:] != pattern[:-1])
    uidx = np.where(uniq)[0]
    uend = np.r_[uidx[1:], rb_bin.shape[0]]
    bins_u = rb_bin[uidx]
    cnt_u = (uend - uidx).astype(np.int64)
    mc_sum_u = np.zeros(uidx.shape[0], dtype=np.int64)
    hmc_sum_u = np.zeros(uidx.shape[0], dtype=np.int64)
    for i in range(uidx.shape[0]):
        mc_sum_u[i] = mc_count[uidx[i]:uend[i]].sum()
        hmc_sum_u[i] = hmc_count[uidx[i]:uend[i]].sum()

    buniq = np.ones(bins_u.shape[0], dtype=bool)
    buniq[1:] = (bins_u[1:] != bins_u[:-1])
    bidx = np.where(buniq)[0]
    bend = np.r_[bidx[1:], bins_u.shape[0]]

    cov_arr = np.zeros(n_full, dtype=np.int64)
    me_arr = np.full(n_full, np.nan, dtype=np.float64)
    mml_arr = np.full(n_full, np.nan, dtype=np.float64)
    mhml_arr = np.full(n_full, np.nan, dtype=np.float64)
    cov5mc = np.zeros(n_full, dtype=np.int64)
    cov5hmc = np.zeros(n_full, dtype=np.int64)
    max_H = float(k) * np.log2(3.0)

    for i in range(bidx.shape[0]):
        s0, s1 = bidx[i], bend[i]
        b = int(bins_u[s0])
        counts = cnt_u[s0:s1]
        coverage = int(counts.sum())
        cov_arr[b] = coverage
        total_mc = int(mc_sum_u[s0:s1].sum())
        total_hmc = int(hmc_sum_u[s0:s1].sum())
        cov5mc[b] = total_mc
        cov5hmc[b] = total_hmc
        if coverage < int(min_coverage):
            continue
        mml_arr[b] = total_mc / (coverage * float(k))
        mhml_arr[b] = total_hmc / (coverage * float(k))
        me_arr[b] = _shannon(counts) / max_H

    def _write(path: str, arr: np.ndarray, is_int: bool = False) -> None:
        with open(path, "w") as f:
            for i in range(n_full):
                if is_int:
                    f.write(f"{chrom}\t{int(starts[i])}\t{int(ends[i])}\t{int(arr[i])}\n")
                elif not np.isnan(arr[i]):
                    f.write(f"{chrom}\t{int(starts[i])}\t{int(ends[i])}\t{arr[i]:.6f}\n")

    _write(out_paths["coverage"], cov_arr, is_int=True)
    _write(out_paths["me"], me_arr)
    _write(out_paths["mml"], mml_arr)
    _write(out_paths["mhml"], mhml_arr)
    _write(out_paths["cov_5mc"], cov5mc, is_int=True)
    _write(out_paths["cov_5hmc"], cov5hmc, is_int=True)
    return f"[INFO] {chrom}: done ({len(rows):,} rows, {int(cov_arr.sum()):,} coverage)"


def concat_ternary_outputs(
    chroms: list[str], out_chrom_dir: str, out_prefix: str
) -> None:
    suffixes = ["coverage", "me", "mml", "mhml", "cov_5mc", "cov_5hmc"]
    _ensure_dir(os.path.dirname(out_prefix))
    for suffix in suffixes:
        out_path = f"{out_prefix}.{suffix}.bedgraph"
        with open(out_path, "w") as fout:
            for chrom in chroms:
                p = os.path.join(out_chrom_dir, f"{chrom}.{suffix}.bedgraph")
                if os.path.exists(p):
                    with open(p) as fin:
                        fout.write(fin.read())


def run_whole_genome_ternary(
    tsv_path: str,
    fasta_path: str,
    out_prefix: str,
    work_dir: str,
    threads: int,
    cpgs_per_bin: int,
    methyl_thresh: float,
    min_coverage: int,
    chroms: str = "",
    fasta_chunk: int = 10_000_000,
    force_ingest: bool = False,
    duckdb_table: str = "reads_ternary",
    duckdb_path: str = "",
    tmp_dir: str = "",
    progress_cb: Callable[[str], None] | None = None,
    pct_cb: Callable[[float, str], None] | None = None,
) -> None:
    """Streamlit-friendly wrapper for the full ternary pipeline."""
    tsv_path = os.path.abspath(tsv_path)
    fasta_path = os.path.abspath(fasta_path)
    out_prefix = os.path.abspath(out_prefix)
    work_dir = os.path.abspath(work_dir)

    _ensure_dir(os.path.dirname(out_prefix))
    _ensure_dir(work_dir)
    db_path = os.path.abspath(duckdb_path) if duckdb_path else os.path.join(work_dir, "reads_ternary.duckdb")
    tmp = os.path.abspath(tmp_dir) if tmp_dir else os.path.join(work_dir, "duckdb_tmp")
    out_chrom_dir = os.path.join(work_dir, "chrom_bedgraphs_ternary")
    _ensure_dir(tmp)
    _ensure_dir(out_chrom_dir)

    suffixes = ["coverage", "me", "mml", "mhml", "cov_5mc", "cov_5hmc"]
    final_paths = [f"{out_prefix}.{s}.bedgraph" for s in suffixes]
    if all(os.path.exists(p) and os.path.getsize(p) > 0 for p in final_paths):
        if progress_cb:
            progress_cb("[INFO] Ternary outputs already exist; skipping.")
        return

    chrom_list = (
        [c.strip() for c in chroms.split(",") if c.strip()]
        if chroms.strip()
        else _list_fasta_chroms(fasta_path)
    )

    if progress_cb:
        progress_cb("[ternary] Stage A: ingest TSV with m/h pivot")
    stage_a_build_duckdb_ternary(
        tsv_path=tsv_path, db_path=db_path, table_name=duckdb_table,
        tmp_dir=tmp, threads=int(threads),
        methyl_thresh=float(methyl_thresh), force=bool(force_ingest),
    )

    nproc = max(1, min(int(threads), len(chrom_list)))
    if progress_cb:
        progress_cb(f"[ternary] Stage B: per-chrom metrics with {nproc} workers")
    args = [
        (chrom, fasta_path, db_path, duckdb_table, out_chrom_dir,
         int(cpgs_per_bin), int(min_coverage), int(fasta_chunk))
        for chrom in chrom_list
    ]
    total = len(args)
    done = 0
    if pct_cb:
        pct_cb(0.0, f"0/{total} chroms")
    with mp.get_context("spawn").Pool(processes=nproc) as pool:
        for msg in pool.imap_unordered(_ternary_worker_unpack, args):
            done += 1
            if progress_cb:
                progress_cb(str(msg))
            if pct_cb:
                pct_cb(done / total, f"{done}/{total} chroms")

    if progress_cb:
        progress_cb("[ternary] Stage C: concatenating per-chrom outputs")
    concat_ternary_outputs(chroms=chrom_list, out_chrom_dir=out_chrom_dir, out_prefix=out_prefix)


def _ternary_worker_unpack(args: tuple) -> str:
    return compute_chrom_ternary(*args)


def required_coverage_for_k(cpgs_per_bin: int) -> int:
    """Minimum recommended coverage for ternary entropy at k CpGs/bin.

    3^k is the cap on distinct ternary patterns; we need at least that
    much coverage to estimate the entropy without severe finite-sample
    bias.
    """
    return int(3 ** int(cpgs_per_bin))
