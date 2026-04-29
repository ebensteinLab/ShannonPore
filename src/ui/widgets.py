"""Reusable Streamlit widgets (file pickers, ROI selectors, etc.).

Stubs for now — populated as tabs are migrated.
"""

from __future__ import annotations

from pathlib import Path


def file_picker(label: str, key: str, default: Path | None = None) -> Path | None:
    """Render a text input + 'Browse' button for a filesystem path."""
    import streamlit as st

    default_str = str(default) if default else ""
    raw = st.text_input(label, value=default_str, key=key)
    return Path(raw).expanduser() if raw else None
