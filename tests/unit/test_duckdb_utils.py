"""Tests for shannonpore.pipelines.duckdb_utils: ingest locking and
memory-bounded connection configuration."""

from __future__ import annotations

import threading

import duckdb
import pytest

from shannonpore.pipelines.duckdb_utils import (
    configure_connection,
    ingest_lock,
    ingest_memory_limit_gb,
    worker_memory_limit_gb,
)


class TestIngestLock:
    def test_reacquire_after_release(self, tmp_path):
        db = str(tmp_path / "a.duckdb")
        with ingest_lock(db):
            pass
        with ingest_lock(db):  # must not raise
            pass

    def test_concurrent_acquire_raises(self, tmp_path):
        db = str(tmp_path / "a.duckdb")
        entered = threading.Event()
        release = threading.Event()

        def holder():
            with ingest_lock(db):
                entered.set()
                release.wait(timeout=10)

        t = threading.Thread(target=holder)
        t.start()
        assert entered.wait(timeout=10)
        try:
            with pytest.raises(RuntimeError, match="already ingesting"), ingest_lock(db):
                pass
        finally:
            release.set()
            t.join(timeout=10)

        # After the holder releases, the lock is available again.
        with ingest_lock(db):
            pass

    def test_different_dbs_do_not_contend(self, tmp_path):
        entered = threading.Event()
        release = threading.Event()

        def holder():
            with ingest_lock(str(tmp_path / "a.duckdb")):
                entered.set()
                release.wait(timeout=10)

        t = threading.Thread(target=holder)
        t.start()
        assert entered.wait(timeout=10)
        try:
            with ingest_lock(str(tmp_path / "b.duckdb")):  # must not raise
                pass
        finally:
            release.set()
            t.join(timeout=10)


class TestMemoryLimits:
    def test_positive_limits(self):
        assert ingest_memory_limit_gb() >= 1
        assert worker_memory_limit_gb(8) >= 1
        assert worker_memory_limit_gb(0) >= 1  # guards divide-by-zero

    def test_workers_get_smaller_share_than_ingest(self):
        assert worker_memory_limit_gb(8) <= ingest_memory_limit_gb()

    def test_configure_connection_applies(self, tmp_path):
        con = duckdb.connect(str(tmp_path / "c.duckdb"))
        configure_connection(
            con,
            threads=2,
            memory_limit_gb=1,
            tmp_dir=str(tmp_path / "spill"),
        )
        threads = con.execute("SELECT current_setting('threads')").fetchone()[0]
        mem = con.execute("SELECT current_setting('memory_limit')").fetchone()[0]
        con.close()
        assert int(threads) == 2
        # 1 GB is echoed back by DuckDB as e.g. "953.6 MiB" — just check
        # a limit is set and it's in the right ballpark (< 2 GiB).
        value, unit = str(mem).split()
        assert unit in {"MiB", "GiB", "MB", "GB"}
        assert float(value) <= 2048.0
