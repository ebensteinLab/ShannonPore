"""Unified progress reporting for shannonpore.

Two adapters share one abstract interface so pipeline code stays
agnostic to whether it's running under Streamlit or a Bash CLI:

- ``CLIProgress``       wraps tqdm — produces a real terminal bar.
- ``StreamlitProgress`` wraps ``st.progress`` + ``st.status``.
- ``NullProgress``      no-op fallback for tests / library use.

Pipelines accept two optional callbacks:
    progress_cb(msg: str)              — text status line
    pct_cb(fraction: float, msg: str)  — fractional progress 0..1
The Progress objects below expose those as bound methods you can pass
straight in.
"""

from __future__ import annotations

import contextlib
import logging
from contextlib import AbstractContextManager
from typing import Any

logger = logging.getLogger(__name__)


# ─── Public protocol ──────────────────────────────────────────────────────

class Progress(AbstractContextManager):
    """Common interface; concrete subclasses implement these methods."""

    def status(self, msg: str) -> None: ...
    def status_line(self, msg: str) -> None:
        """Replace the most recent status line in place (for ephemeral
        progress-bar frames terminated by ``\\r``). Default: same as
        ``status`` — subclasses override for in-place updates."""
        self.status(msg)
    def update(self, fraction: float, msg: str = "") -> None: ...
    def close(self) -> None: ...

    # Convenience: bound callback factories for pipeline functions.
    @property
    def progress_cb(self):
        return self.status

    @property
    def status_cb(self):
        return self.status_line

    @property
    def pct_cb(self):
        return self.update

    def __exit__(self, *exc: Any) -> None:
        self.close()


# ─── No-op ────────────────────────────────────────────────────────────────

class NullProgress(Progress):
    def status(self, msg: str) -> None:
        logger.debug("%s", msg)

    def update(self, fraction: float, msg: str = "") -> None:
        pass

    def close(self) -> None:
        pass


# ─── tqdm (CLI) ───────────────────────────────────────────────────────────

class CLIProgress(Progress):
    """tqdm-backed progress bar for the CLI.

    `total` of 1.0 means we report fractional progress; the bar shows a
    percentage. Status lines are written ABOVE the bar via tqdm.write().
    """

    def __init__(self, label: str, *, ncols: int = 80) -> None:
        from tqdm import tqdm

        self._label = label
        self._tqdm = tqdm(
            total=100, desc=label, ncols=ncols,
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}% {elapsed}",
            leave=True,
        )
        self._last_pct = 0
        self._tqdm_write = tqdm.write

    def status(self, msg: str) -> None:
        # Write status lines above the bar so they don't disrupt rendering.
        try:
            self._tqdm_write(str(msg))
        except (OSError, ValueError):
            print(msg, flush=True)

    def status_line(self, msg: str) -> None:
        """In-place ephemeral status (e.g. modkit's progress bar frames).

        Stuffs the latest frame into tqdm's right-hand postfix so it
        updates in place rather than scrolling — exactly what the user
        wants for a remote tool's own progress bar.
        """
        # Trim to keep tqdm's line single-row even on wide bars.
        with contextlib.suppress(OSError, ValueError):
            self._tqdm.set_postfix_str(str(msg)[:80], refresh=True)

    def update(self, fraction: float, msg: str = "") -> None:
        pct = max(0, min(100, int(round(fraction * 100))))
        delta = pct - self._last_pct
        if delta > 0:
            self._tqdm.update(delta)
            self._last_pct = pct
        if msg:
            self._tqdm.set_postfix_str(msg, refresh=False)

    def close(self) -> None:
        # Snap to 100% on clean close.
        with contextlib.suppress(OSError, ValueError):
            self._tqdm.update(max(0, 100 - self._last_pct))
        self._tqdm.close()


# ─── Streamlit ────────────────────────────────────────────────────────────

class StreamlitProgress(Progress):
    """Streamlit progress bar + status panel.

    Use as a context manager so the status panel collapses on exit:

        with StreamlitProgress("modkit extract") as p:
            run_pipeline(progress_cb=p.progress_cb, pct_cb=p.pct_cb)
    """

    def __init__(self, label: str) -> None:
        import streamlit as st

        self._st = st
        self._status = st.status(label, expanded=True, state="running")
        self._status.__enter__()
        # Inside the status panel we render a progress bar + a small log.
        self._bar = self._status.progress(0.0, text=label)
        self._log_lines: list[str] = []
        self._log_box = self._status.empty()
        # Dedicated slot for in-place status (modkit's progress bar etc.).
        self._status_line_box = self._status.empty()
        self._label = label

    def status(self, msg: str) -> None:
        if not msg:
            return
        self._log_lines.append(str(msg))
        # Show only the last 12 lines in monospaced code block.
        tail = "\n".join(self._log_lines[-12:])
        self._log_box.code(tail, language="text")

    def status_line(self, msg: str) -> None:
        """In-place status — for ``\\r``-terminated progress-bar frames
        from external tools (modkit). Replaces the previous frame
        rather than appending, so the user sees a live updating
        single line instead of a cascade of near-identical entries."""
        if not msg:
            return
        self._status_line_box.code(str(msg), language="text")

    def update(self, fraction: float, msg: str = "") -> None:
        f = max(0.0, min(1.0, float(fraction)))
        text = f"{self._label} · {msg}" if msg else self._label
        self._bar.progress(f, text=text)

    def close(self) -> None:
        try:
            self._bar.progress(1.0, text=f"{self._label} · done")
            self._status.update(state="complete", expanded=False)
            self._status.__exit__(None, None, None)
        except Exception:  # noqa: BLE001
            pass


# ─── Backwards-compat helpers used by older modules ─────────────────────

def status(label: str, expanded: bool = True):
    """Legacy wrapper around `st.status` kept for callers that don't need
    the full `StreamlitProgress` adapter."""
    import streamlit as st
    return st.status(label, expanded=expanded)


def progress_bar(label: str = "Working...") -> object:
    import streamlit as st
    return st.progress(0.0, text=label)
