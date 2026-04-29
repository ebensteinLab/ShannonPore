"""Whole-genome methylation entropy via DuckDB + per-chrom multiprocessing.

Stage A: ingest modkit TSV into DuckDB grouped by (read_id, chrom, ref_position).
Stage B: per chromosome, find CpG positions, intersect with reads, compute
per-bin coverage / MML / ME, write per-chrom bedgraphs.
Stage C: concatenate per-chrom bedgraphs into a single output prefix.

v4 changes from v3:
- Removed `SystemExit` raised at import time on missing duckdb/pyfaidx.
  v4 expects these via the conda lockfile; if missing the user gets a
  normal ImportError, not a SystemExit that would crash Streamlit.
- Replaced bare except clauses with specific exceptions.
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


def ensure_dir(path: str) -> None:
    if path == "":
        return
    os.makedirs(path, exist_ok=True)


def list_fasta_chroms(fasta_path: str) -> list[str]:
    return list(Fasta(fasta_path).keys())


def shannon_entropy_from_counts(counts: np.ndarray) -> float:
    total = int(counts.sum())
    if total <= 0:
        return 0.0
    probs = counts.astype(np.float64) / float(total)
    probs = probs[probs > 0]
    return float(-np.sum(probs * np.log2(probs)))


def find_cpg_positions_chrom_chunked(
    fasta_path: str, chrom: str, chunk_size: int = 10_000_000
) -> np.ndarray:
    genome = Fasta(fasta_path)
    chrom_len = len(genome[chrom])

    cpg_positions: list[int] = []
    carry = ""
    pos = 0
    while pos < chrom_len:
        end = min(pos + chunk_size, chrom_len)
        seq = carry + genome[chrom][pos:end].seq.upper()
        offset = pos - len(carry)
        for m in re.finditer("CG", seq):
            cpg_positions.append(offset + m.start())
        carry = seq[-1:] if len(seq) > 0 else ""
        pos = end

    if not cpg_positions:
        return np.array([], dtype=np.int64)
    return np.array(cpg_positions, dtype=np.int64)


def pd_factorize(values: np.ndarray) -> tuple[np.ndarray, int]:
    uniques, inv = np.unique(values, return_inverse=True)
    return inv.astype(np.int64, copy=False), int(uniques.shape[0])


def stage_a_build_duckdb_table(
    tsv_path: str,
    db_path: str,
    table_name: str,
    tmp_dir: str,
    threads: int,
    force: bool,
    entropy_mode: str = "true_mc",
) -> None:
    """Stage A: ingest modkit TSV into DuckDB.

    `entropy_mode` controls how 5hmC calls are interpreted:
      - "true_mc"   : drop 'h' rows; treat 5hmC sites as unmodified C.
                      mod_qual reflects 5mC probability only.
      - "bisulfite" : keep both 'm' and 'h' rows; AVG(mod_qual) sums their
                      contribution so 5hmC counts as 5mC (BS-equivalent).
                      This is what v3 did.
      - "ternary"   : NOT handled here — use src/pipelines/ternary_entropy.py.
                      Raises ValueError if requested.
    """
    if entropy_mode == "ternary":
        raise ValueError(
            "Use src/pipelines/ternary_entropy.run_whole_genome_ternary "
            "for ternary mode; this binary stage handles 2-state modes only."
        )

    ensure_dir(os.path.dirname(db_path))
    ensure_dir(tmp_dir)

    if not os.path.exists(tsv_path):
        raise FileNotFoundError(f"TSV not found: {tsv_path}")
    with open(tsv_path, "rb") as fh:
        head = fh.read(4096)
    if b"\t" not in head:
        raise ValueError(
            "Input TSV does not look tab-delimited (no TAB found in first 4KB). "
            "Your modkit extract TSV should be TAB-separated."
        )

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

    # mode-dependent WHERE clause
    if entropy_mode == "true_mc":
        mode_filter = "AND mod_code = 'm'"
    elif entropy_mode == "bisulfite":
        mode_filter = "AND mod_code IN ('m', 'h')"
    else:
        raise ValueError(f"Unknown entropy_mode: {entropy_mode!r}")

    con.execute(
        f"""
        CREATE TABLE {table_name} AS
        WITH raw AS (
            SELECT
                read_id::VARCHAR AS read_id,
                chrom::VARCHAR AS chrom,
                CASE
                    WHEN ref_strand = '-' THEN (ref_position::BIGINT - 1)
                    ELSE ref_position::BIGINT
                END AS ref_position,
                mod_qual::DOUBLE AS mod_qual,
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
              {mode_filter}
        )
        SELECT
            read_id,
            chrom,
            ref_position,
            SUM(mod_qual) AS mod_qual
        FROM raw
        GROUP BY 1,2,3;
        """
    )
    con.execute(f"CREATE INDEX idx_{table_name}_chrom_pos ON {table_name}(chrom, ref_position);")
    con.execute(f"CREATE INDEX idx_{table_name}_chrom ON {table_name}(chrom);")
    con.close()


def compute_chrom_metrics_from_db(
    chrom: str,
    fasta_path: str,
    db_path: str,
    table_name: str,
    out_chrom_dir: str,
    cpgs_per_bin: int,
    methyl_thresh: float,
    min_coverage: int,
    fasta_chunk: int,
) -> str:
    ensure_dir(out_chrom_dir)
    out_cov = os.path.join(out_chrom_dir, f"{chrom}.coverage.bedgraph")
    out_mml = os.path.join(out_chrom_dir, f"{chrom}.mml.bedgraph")
    out_me = os.path.join(out_chrom_dir, f"{chrom}.me.bedgraph")

    if os.path.exists(out_cov) and os.path.exists(out_mml) and os.path.exists(out_me):
        return f"[INFO] {chrom}: outputs exist, skipping"

    cpg_positions = find_cpg_positions_chrom_chunked(
        fasta_path, chrom, chunk_size=int(fasta_chunk)
    )
    k = int(cpgs_per_bin)
    if cpg_positions.size < k:
        for p in (out_cov, out_mml, out_me):
            open(p, "w").close()
        return f"[INFO] {chrom}: not enough CpGs, wrote empty outputs"

    n_full_bins = int(cpg_positions.size // k)
    starts = cpg_positions[0: n_full_bins * k: k]
    ends = cpg_positions[(k - 1): (n_full_bins * k): k] + 2

    con = duckdb.connect(db_path, read_only=True)
    rows = con.execute(
        f"SELECT read_id, ref_position, mod_qual FROM {table_name} WHERE chrom = ?",
        [chrom],
    ).fetchall()
    con.close()

    def _write_zero_coverage(reason: str) -> str:
        with open(out_cov, "w") as f_cov:
            for i in range(n_full_bins):
                f_cov.write(f"{chrom}\t{int(starts[i])}\t{int(ends[i])}\t0\n")
        open(out_mml, "w").close()
        open(out_me, "w").close()
        return f"[INFO] {chrom}: {reason}"

    if not rows:
        return _write_zero_coverage("no events, wrote zero-coverage outputs")

    read_id = np.array([r[0] for r in rows], dtype=object)
    pos = np.array([r[1] for r in rows], dtype=np.int64)
    mod = np.array([r[2] for r in rows], dtype=np.float32)
    read_code, _ = pd_factorize(read_id)

    idx = np.searchsorted(cpg_positions, pos)
    valid = (idx >= 0) & (idx < cpg_positions.size)
    idx2 = idx[valid]
    matches = cpg_positions[idx2] == pos[valid]

    if not np.any(matches):
        return _write_zero_coverage("no CpG-matching events")

    ordinal = idx2[matches].astype(np.int64)
    mod = mod[valid][matches]
    read_code = read_code[valid][matches].astype(np.int64)

    bin_id = (ordinal // k).astype(np.int64)
    in_full_bins = bin_id < n_full_bins
    if not np.any(in_full_bins):
        return _write_zero_coverage("CpG events but no full bins")

    ordinal = ordinal[in_full_bins]
    bin_id = bin_id[in_full_bins]
    mod = mod[in_full_bins]
    read_code = read_code[in_full_bins]

    cpg_idx = (ordinal % k).astype(np.int8)
    shift = (k - 1 - cpg_idx).astype(np.int8)
    status = (mod >= float(methyl_thresh)).astype(np.int8)

    keys = np.rec.fromarrays([read_code, bin_id, shift], names="r,b,s")
    order = np.argsort(keys, kind="mergesort")
    keys = keys[order]
    status = status[order]

    uniq_mask = np.ones(keys.shape[0], dtype=bool)
    uniq_mask[1:] = (keys[1:] != keys[:-1])
    uniq_idx = np.where(uniq_mask)[0]
    group_ends = np.r_[uniq_idx[1:], keys.shape[0]]

    collapsed_read = read_code[order][uniq_idx]
    collapsed_bin = bin_id[order][uniq_idx]
    collapsed_shift = shift[order][uniq_idx]

    group_max = np.zeros(uniq_idx.shape[0], dtype=np.int8)
    for i in range(uniq_idx.shape[0]):
        s0, s1 = uniq_idx[i], group_ends[i]
        group_max[i] = 1 if np.any(status[s0:s1]) else 0

    full_mask = (1 << k) - 1
    bit = (group_max.astype(np.int64) << collapsed_shift.astype(np.int64)).astype(np.int64)
    obs = (1 << collapsed_shift.astype(np.int64)).astype(np.int64)

    rb_keys = np.rec.fromarrays([collapsed_read, collapsed_bin], names="r,b")
    rb_order = np.argsort(rb_keys, kind="mergesort")
    rb_keys = rb_keys[rb_order]
    bit = bit[rb_order]
    obs = obs[rb_order]

    rb_uniq = np.ones(rb_keys.shape[0], dtype=bool)
    rb_uniq[1:] = (rb_keys[1:] != rb_keys[:-1])
    rb_idx = np.where(rb_uniq)[0]
    rb_end = np.r_[rb_idx[1:], rb_keys.shape[0]]

    rb_bin = collapsed_bin[rb_order][rb_idx]
    pattern = np.zeros(rb_idx.shape[0], dtype=np.int64)
    obsmask = np.zeros(rb_idx.shape[0], dtype=np.int64)
    for i in range(rb_idx.shape[0]):
        s0, s1 = rb_idx[i], rb_end[i]
        pattern[i] = bit[s0:s1].sum()
        obsmask[i] = obs[s0:s1].sum()

    full = (obsmask == full_mask)
    if not np.any(full):
        return _write_zero_coverage("no full k-CpG patterns")

    rb_bin = rb_bin[full]
    pattern = pattern[full]

    pair = np.rec.fromarrays([rb_bin, pattern], names="b,p")
    pair_order = np.argsort(pair, kind="mergesort")
    rb_bin = rb_bin[pair_order]
    pattern = pattern[pair_order]

    uniq = np.ones(rb_bin.shape[0], dtype=bool)
    uniq[1:] = (rb_bin[1:] != rb_bin[:-1]) | (pattern[1:] != pattern[:-1])
    uidx = np.where(uniq)[0]
    uend = np.r_[uidx[1:], rb_bin.shape[0]]

    bins_u = rb_bin[uidx]
    pat_u = pattern[uidx]
    cnt_u = (uend - uidx).astype(np.int64)

    buniq = np.ones(bins_u.shape[0], dtype=bool)
    buniq[1:] = (bins_u[1:] != bins_u[:-1])
    bidx = np.where(buniq)[0]
    bend = np.r_[bidx[1:], bins_u.shape[0]]

    cov_arr = np.zeros(n_full_bins, dtype=np.int64)
    mml_arr = np.full(n_full_bins, np.nan, dtype=np.float64)
    me_arr = np.full(n_full_bins, np.nan, dtype=np.float64)

    max_pat = 1 << k
    ones = np.array([int(bin(x).count("1")) for x in range(max_pat)], dtype=np.int64)

    for i in range(bidx.shape[0]):
        s0, s1 = bidx[i], bend[i]
        b = int(bins_u[s0])
        counts = cnt_u[s0:s1]
        pats = pat_u[s0:s1]

        coverage = int(counts.sum())
        cov_arr[b] = coverage
        if coverage < int(min_coverage):
            continue

        ones_weighted = int((ones[pats] * counts).sum())
        mml_arr[b] = ones_weighted / (coverage * float(k))
        me_arr[b] = shannon_entropy_from_counts(counts) / float(k)

    with open(out_cov, "w") as f_cov:
        for i in range(n_full_bins):
            f_cov.write(f"{chrom}\t{int(starts[i])}\t{int(ends[i])}\t{int(cov_arr[i])}\n")
    with open(out_mml, "w") as f_mml:
        for i in range(n_full_bins):
            if not np.isnan(mml_arr[i]):
                f_mml.write(f"{chrom}\t{int(starts[i])}\t{int(ends[i])}\t{float(mml_arr[i]):.6f}\n")
    with open(out_me, "w") as f_me:
        for i in range(n_full_bins):
            if not np.isnan(me_arr[i]):
                f_me.write(f"{chrom}\t{int(starts[i])}\t{int(ends[i])}\t{float(me_arr[i]):.6f}\n")

    return f"[INFO] {chrom}: done (rows fetched: {len(rows)})"


def concat_chrom_outputs(chroms: list[str], out_chrom_dir: str, out_prefix: str) -> None:
    out_cov = f"{out_prefix}.coverage.bedgraph"
    out_mml = f"{out_prefix}.mml.bedgraph"
    out_me = f"{out_prefix}.me.bedgraph"
    ensure_dir(os.path.dirname(out_cov))

    def _concat(suffix: str, dest_file) -> None:
        for chrom in chroms:
            p = os.path.join(out_chrom_dir, f"{chrom}.{suffix}.bedgraph")
            if os.path.exists(p):
                with open(p, "r", encoding="utf-8", errors="replace") as r:
                    dest_file.write(r.read())

    with open(out_cov, "w") as f_cov, open(out_mml, "w") as f_mml, open(out_me, "w") as f_me:
        _concat("coverage", f_cov)
        _concat("mml", f_mml)
        _concat("me", f_me)


def run_whole_genome_duckdb_only(
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
    duckdb_table: str = "reads_agg",
    duckdb_path: str = "",
    tmp_dir: str = "",
    progress_cb: Callable[[str], None] | None = None,
    pct_cb: Callable[[float, str], None] | None = None,
    entropy_mode: str = "true_mc",
) -> None:
    tsv_path = os.path.abspath(tsv_path)
    fasta_path = os.path.abspath(fasta_path)
    out_prefix = os.path.abspath(out_prefix)
    work_dir = os.path.abspath(work_dir)

    ensure_dir(os.path.dirname(out_prefix))
    ensure_dir(work_dir)

    db_path = os.path.abspath(duckdb_path) if duckdb_path else os.path.join(work_dir, "reads.duckdb")
    tmp_dir2 = os.path.abspath(tmp_dir) if tmp_dir else os.path.join(work_dir, "duckdb_tmp")
    out_chrom_dir = os.path.join(work_dir, "chrom_bedgraphs")
    ensure_dir(tmp_dir2)
    ensure_dir(out_chrom_dir)

    final_cov = f"{out_prefix}.coverage.bedgraph"
    final_mml = f"{out_prefix}.mml.bedgraph"
    final_me = f"{out_prefix}.me.bedgraph"
    if all(os.path.exists(p) and os.path.getsize(p) > 0 for p in (final_cov, final_mml, final_me)):
        if progress_cb:
            progress_cb("[INFO] Whole-genome final bedGraphs already exist, skipping computation.")
        return

    chrom_list = (
        [c.strip() for c in chroms.split(",") if c.strip()]
        if chroms.strip()
        else list_fasta_chroms(fasta_path)
    )

    if progress_cb:
        progress_cb(f"[WG] Stage A: ingest TSV into DuckDB (mode={entropy_mode})")
    stage_a_build_duckdb_table(
        tsv_path=tsv_path,
        db_path=db_path,
        table_name=duckdb_table,
        tmp_dir=tmp_dir2,
        threads=int(threads),
        force=bool(force_ingest),
        entropy_mode=entropy_mode,
    )

    nproc = max(1, min(int(threads), len(chrom_list)))
    if progress_cb:
        progress_cb(f"[WG] Stage B: per-chrom metrics with {nproc} processes")

    worker_args = [
        (chrom, fasta_path, db_path, duckdb_table, out_chrom_dir,
         int(cpgs_per_bin), float(methyl_thresh), int(min_coverage), int(fasta_chunk))
        for chrom in chrom_list
    ]

    total = len(worker_args)
    done = 0
    if pct_cb:
        pct_cb(0.0, f"0/{total} chroms")
    with mp.get_context("spawn").Pool(processes=nproc) as pool:
        for msg in pool.imap_unordered(_chrom_worker_unpack, worker_args):
            done += 1
            if progress_cb:
                progress_cb(str(msg))
            if pct_cb:
                pct_cb(done / total, f"{done}/{total} chroms")

    if progress_cb:
        progress_cb("[WG] Concatenating per-chrom outputs")
    concat_chrom_outputs(chroms=chrom_list, out_chrom_dir=out_chrom_dir, out_prefix=out_prefix)


def _chrom_worker_unpack(args: tuple) -> str:
    """Adapter for `imap_unordered` — accepts a tuple instead of *args."""
    return compute_chrom_metrics_from_db(*args)
