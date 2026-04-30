"""Unit tests for `src/pipelines/bam_utils.py`."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from src.pipelines.bam_utils import find_bams


@pytest.mark.unit
def test_find_bams_returns_only_bams(tmp_path: Path) -> None:
    (tmp_path / "a.bam").write_bytes(b"x")
    (tmp_path / "b.bam").write_bytes(b"x")
    (tmp_path / "c.sam").write_text("not a bam")
    (tmp_path / "d.bam.bai").write_bytes(b"index")
    found = find_bams(tmp_path)
    names = [p.name for p in found]
    assert "a.bam" in names and "b.bam" in names
    assert "c.sam" not in names
    assert all(not n.endswith(".bai") for n in names)


@pytest.mark.unit
def test_find_bams_recursive(tmp_path: Path) -> None:
    sub = tmp_path / "nested"
    sub.mkdir()
    (sub / "x.bam").write_bytes(b"x")
    (tmp_path / "y.bam").write_bytes(b"x")
    found = find_bams(tmp_path, recursive=True)
    assert len(found) == 2
    found_flat = find_bams(tmp_path, recursive=False)
    assert len(found_flat) == 1


@pytest.mark.unit
def test_find_bams_results_sorted(tmp_path: Path) -> None:
    for name in ("z.bam", "a.bam", "m.bam"):
        (tmp_path / name).write_bytes(b"x")
    names = [p.name for p in find_bams(tmp_path)]
    assert names == sorted(names)


@pytest.mark.unit
def test_find_bams_empty_folder_returns_empty_list(tmp_path: Path) -> None:
    assert find_bams(tmp_path) == []


@pytest.mark.unit
def test_find_bams_missing_folder_raises(tmp_path: Path) -> None:
    with pytest.raises(NotADirectoryError):
        find_bams(tmp_path / "nope")


@pytest.mark.unit
@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")
def test_merge_sort_index_single_bam(tmp_path: Path) -> None:
    """One BAM: no merge step, just sort+index."""
    pysam = pytest.importorskip("pysam")
    from src.pipelines.bam_utils import merge_sort_index_bams

    # Build a tiny BAM
    bam = tmp_path / "in.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"LN": 100, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as bf:
        a = pysam.AlignedSegment(bf.header)
        a.query_name = "r1"
        a.query_sequence = "A" * 10
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 0
        a.cigar = [(0, 10)]
        a.mapping_quality = 30
        a.query_qualities = pysam.qualitystring_to_array("I" * 10)
        bf.write(a)

    out = tmp_path / "out.sorted.bam"
    res = merge_sort_index_bams([bam], out, threads=1)
    assert res == out
    assert (tmp_path / "out.sorted.bam.bai").exists()
