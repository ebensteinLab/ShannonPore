"""Shared DuckDB helpers: ingest locking and memory-bounded connections.

Why the lock exists
-------------------
The Streamlit UI runs the pipeline inside the script thread. A page
reconnect / rerun (or a second click on RUN while a long job is still
going) starts a *second* pipeline in the same process. Both threads then
race to ``CREATE TABLE`` in the same ``.duckdb`` file and DuckDB aborts
one of them with::

    TransactionContext Error: Catalog write-write conflict on create
    with "Schema\\0main\\0main\\0Table\\0main\\0reads_agg"

``ingest_lock`` serializes Stage-A ingest per database file. It uses an
``flock`` on a sidecar ``<db>.ingest.lock`` file, which contends both
across threads (separate fds) and across processes, and is released
automatically by the OS if the holder dies — no stale-lock cleanup
needed. On platforms without ``fcntl`` a process-wide registry of
``threading.Lock`` objects still covers the same-process case (which is
the one that actually bites under Streamlit).

Why the memory settings exist
-----------------------------
DuckDB's default ``memory_limit`` is 80% of physical RAM *per database
instance*. Stage A (one big GROUP BY over the whole modkit TSV) plus N
Stage-B worker processes each opening the database meant the pipeline
could try to claim several times the machine's RAM. Every connection
now gets an explicit budget so DuckDB spills to ``temp_directory``
instead of thrashing / getting OOM-killed.
"""

from __future__ import annotations

import contextlib
import os
import threading
from collections.abc import Iterator

import duckdb

try:
    import fcntl
except ImportError:  # pragma: no cover - non-POSIX platforms
    fcntl = None  # type: ignore[assignment]

# Process-wide fallback locks, keyed by absolute db path.
_LOCK_REGISTRY: dict[str, threading.Lock] = {}
_REGISTRY_GUARD = threading.Lock()

# Fraction of physical RAM the whole pipeline may hand to DuckDB.
MEMORY_FRACTION = 0.5
_FALLBACK_TOTAL_RAM_GB = 16.0


def total_ram_gb() -> float:
    """Physical RAM in GiB; conservative fallback when undetectable."""
    try:
        page = os.sysconf("SC_PAGE_SIZE")
        pages = os.sysconf("SC_PHYS_PAGES")
        return page * pages / 1024**3
    except (ValueError, OSError, AttributeError):  # pragma: no cover
        return _FALLBACK_TOTAL_RAM_GB


def ingest_memory_limit_gb() -> int:
    """DuckDB budget for the single Stage-A ingest connection."""
    return max(1, int(total_ram_gb() * MEMORY_FRACTION))


def worker_memory_limit_gb(nproc: int) -> int:
    """Per-worker DuckDB budget so N concurrent Stage-B readers stay
    within the same overall fraction of RAM as Stage A."""
    return max(1, int(total_ram_gb() * MEMORY_FRACTION / max(1, int(nproc))))


def configure_connection(
    con: duckdb.DuckDBPyConnection,
    *,
    threads: int,
    memory_limit_gb: int,
    tmp_dir: str = "",
) -> None:
    """Apply the standard pragmas: temp spill dir, thread count, memory
    cap, and unordered execution (lower memory for bulk CTAS)."""
    if tmp_dir:
        tmp_sql = str(tmp_dir).replace("'", "''")
        con.execute(f"PRAGMA temp_directory='{tmp_sql}';")
    con.execute(f"PRAGMA threads={int(threads)};")
    con.execute(f"PRAGMA memory_limit='{int(memory_limit_gb)}GB';")
    con.execute("SET preserve_insertion_order=false;")


@contextlib.contextmanager
def ingest_lock(db_path: str) -> Iterator[None]:
    """Exclusive, non-blocking lock guarding Stage-A ingest of `db_path`.

    Raises RuntimeError with an actionable message if another ingest
    into the same database is already running.
    """
    abs_db = os.path.abspath(db_path)
    with _REGISTRY_GUARD:
        thread_lock = _LOCK_REGISTRY.setdefault(abs_db, threading.Lock())

    if not thread_lock.acquire(blocking=False):
        raise RuntimeError(_busy_message(abs_db))

    fd: int | None = None
    try:
        if fcntl is not None:
            fd = os.open(abs_db + ".ingest.lock", os.O_CREAT | os.O_RDWR, 0o644)
            try:
                fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError as exc:
                raise RuntimeError(_busy_message(abs_db)) from exc
        yield
    finally:
        if fd is not None:
            os.close(fd)  # closing the fd releases the flock
        thread_lock.release()


def _busy_message(db_path: str) -> str:
    return (
        f"Another ShannonPore run is already ingesting into '{db_path}'. "
        "This usually means a previous pipeline run is still going "
        "(e.g. the browser tab was reloaded or RUN was clicked twice). "
        "Wait for the running job to finish — its progress continues in "
        "the server log — then re-run; the finished ingest will be "
        "reused automatically."
    )
