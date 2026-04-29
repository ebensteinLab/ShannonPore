"""Wrapper around `modkit extract full` for nanopore CpG methylation calls."""

from __future__ import annotations

import logging
import os
import subprocess
from collections.abc import Callable

from src.io.utils_io import safe_mkdir, validate_existing_path

logger = logging.getLogger(__name__)


def run_modkit_extract_minimal(
    bam_path: str,
    out_tsv_path: str,
    reference_fasta: str,
    threads: int,
    log_filepath: str,
    stream_cb: Callable[[str], None] | None = None,
    status_cb: Callable[[str], None] | None = None,
) -> None:
    """Run `modkit extract full` with minimal args:
    `--reference <FASTA> --cpg --threads N --force --log-filepath LOG`.

    modkit emits two kinds of output on stdout/stderr:

    * **Final lines** end with ``\\n`` — log records like
      "found BAM index, processing reads in 100000 base pair chunks".
      These are sent to ``stream_cb`` and are intended to be appended
      to a scrolling log.
    * **Progress frames** end with ``\\r`` — modkit's `indicatif`
      progress bar redraws itself in place. These are sent to
      ``status_cb`` and should *replace* the previous status line in
      the UI (otherwise the user sees a stack of near-identical
      frames). If ``status_cb`` is None, progress frames are
      forwarded to ``stream_cb`` instead so callers that don't
      distinguish them still see the bar.

    Raises RuntimeError if modkit returns non-zero or produces empty
    output.
    """
    if status_cb is None:
        status_cb = stream_cb
    safe_mkdir(os.path.dirname(out_tsv_path))
    safe_mkdir(os.path.dirname(log_filepath))

    if not out_tsv_path.lower().endswith(".tsv"):
        raise ValueError(f"Output path must end with .tsv: {out_tsv_path}")

    ok = validate_existing_path(reference_fasta, "Reference FASTA")
    if ok is not True:
        raise RuntimeError(str(ok))

    cmd = [
        "modkit", "extract", "full",
        "--reference", reference_fasta,
        "--cpg",
        "--threads", str(int(threads)),
        "--force",
        "--log-filepath", log_filepath,
        bam_path,
        out_tsv_path,
    ]

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        bufsize=0,
    )

    screen_lines: list[str] = []
    current_line = ""

    try:
        assert proc.stdout is not None
        while True:
            chunk = proc.stdout.read(4096)
            if not chunk:
                if proc.poll() is not None:
                    break
                continue

            text = chunk.decode("utf-8", errors="replace")
            for ch in text:
                if ch == "\r":
                    if current_line.strip() and status_cb:
                        status_cb(current_line)
                    current_line = ""
                elif ch == "\n":
                    if current_line.strip():
                        screen_lines.append(current_line)
                        if stream_cb:
                            stream_cb(current_line)
                    current_line = ""
                else:
                    current_line += ch

            if len(screen_lines) > 500:
                screen_lines = screen_lines[-500:]

        rc = proc.wait()
        if rc != 0:
            tail = ""
            if os.path.exists(log_filepath):
                try:
                    with open(log_filepath, errors="replace", encoding="utf-8") as f:
                        tail = "".join(f.readlines()[-120:])
                except OSError as exc:
                    logger.warning("Could not read modkit log %s: %s", log_filepath, exc)
            raise RuntimeError(
                f"modkit extract failed (rc={rc}).\n"
                f"modkit log: {log_filepath}\n"
                f"Log tail:\n{tail}"
            )

        if not os.path.exists(out_tsv_path) or os.path.getsize(out_tsv_path) == 0:
            raise RuntimeError(
                f"modkit finished but TSV output missing or empty: {out_tsv_path}"
            )
    finally:
        # Terminate cleanly and reap the child so we never leak a zombie
        # (Streamlit sessions are long-lived and accumulate them otherwise).
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
