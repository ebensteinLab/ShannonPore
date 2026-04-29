"""BAM helpers: discover, merge, sort, index — via the `samtools` CLI.

Pure-subprocess implementation so the project doesn't add `pysam` as a
runtime dep (pysam is only used in test fixtures).
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from collections.abc import Callable
from pathlib import Path

logger = logging.getLogger(__name__)

ProgressCB = Callable[[str], None]
PctCB = Callable[[float, str], None]


def _samtools() -> str:
    p = shutil.which("samtools")
    if p is None:
        raise RuntimeError(
            "samtools is required but not on PATH. Install via the v4 "
            "conda env (it is pinned in environment.yml at 1.21)."
        )
    return p


def find_bams(folder: str | Path, *, recursive: bool = True) -> list[Path]:
    """Return a sorted list of `*.bam` files under `folder`.

    Skips index files (`.bai`) and SAM/CRAM. Stable across runs (sorted).
    """
    base = Path(folder).expanduser().resolve()
    if not base.is_dir():
        raise NotADirectoryError(f"Not a directory: {base}")

    pattern = "**/*.bam" if recursive else "*.bam"
    bams = sorted(p for p in base.glob(pattern) if p.is_file())
    return bams


def _run(cmd: list[str], *, stream: ProgressCB | None = None) -> None:
    """Run a subprocess, streaming stderr line-by-line through `stream`."""
    logger.info("$ %s", " ".join(cmd))
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        bufsize=1, text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        line = line.rstrip("\n")
        if stream:
            stream(line)
        else:
            logger.debug("%s", line)
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f"samtools failed (rc={rc}): {' '.join(cmd)}")


def merge_sort_index_bams(
    bams: list[str | Path],
    out_bam: str | Path,
    *,
    threads: int = 4,
    overwrite: bool = False,
    progress_cb: ProgressCB | None = None,
    pct_cb: PctCB | None = None,
) -> Path:
    """Merge → sort → index a list of BAM files.

    If `len(bams) == 1`, the merge step is skipped and the single BAM is
    sorted+indexed in place to `out_bam`. Idempotent when
    `overwrite=False` and the indexed BAM already exists.

    Returns the path to the indexed sorted BAM.
    """
    if not bams:
        raise ValueError("merge_sort_index_bams: empty BAM list")

    out_bam = Path(out_bam).expanduser().resolve()
    out_bam.parent.mkdir(parents=True, exist_ok=True)
    bai = out_bam.with_suffix(out_bam.suffix + ".bai")

    if bai.exists() and out_bam.exists() and not overwrite:
        if progress_cb:
            progress_cb(f"[bam] reusing existing indexed BAM: {out_bam}")
        return out_bam

    samtools = _samtools()
    n_steps = 3 if len(bams) > 1 else 2
    step = 0

    def advance(label: str) -> None:
        nonlocal step
        step += 1
        if pct_cb:
            pct_cb(step / n_steps, label)
        if progress_cb:
            progress_cb(f"[bam] {step}/{n_steps} {label}")

    work_dir = out_bam.parent
    merged = out_bam.with_name(out_bam.stem + ".unsorted.bam")

    if len(bams) > 1:
        advance(f"merging {len(bams)} BAMs")
        merge_cmd = [
            samtools, "merge", "-f",
            "-@", str(int(threads)),
            str(merged),
            *[str(p) for p in bams],
        ]
        _run(merge_cmd, stream=progress_cb)
        sort_in = merged
    else:
        sort_in = Path(bams[0])

    advance("sorting")
    sort_cmd = [
        samtools, "sort",
        "-@", str(int(threads)),
        "-o", str(out_bam),
        str(sort_in),
    ]
    _run(sort_cmd, stream=progress_cb)

    if merged.exists():
        try:
            os.remove(merged)
        except OSError:
            pass

    advance("indexing")
    _run([samtools, "index", "-@", str(int(threads)), str(out_bam)],
         stream=progress_cb)

    if not bai.exists():
        raise RuntimeError(f"samtools index produced no .bai for {out_bam}")
    return out_bam
