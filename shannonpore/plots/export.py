"""Multi-format figure export shared by every plot family.

``save_figure`` replaces the per-function ``fig.savefig(out, dpi=200,
bbox_inches="tight")`` calls so all plots agree on how PNG / JPG / SVG /
PDF files are written:

  * raster formats (png, jpg) honour ``dpi``;
  * vector formats (svg, pdf) are resolution-independent — ``dpi`` only
    affects rasterised elements embedded in them;
  * jpg is flattened onto white (JPEG has no alpha) and saved at
    quality 95.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

SUPPORTED_FORMATS: tuple[str, ...] = ("png", "jpg", "svg", "pdf")

# Accepted spellings → canonical extension.
_FORMAT_ALIASES: dict[str, str] = {
    "png": "png",
    "jpg": "jpg",
    "jpeg": "jpg",
    "svg": "svg",
    "pdf": "pdf",
}

_JPEG_QUALITY = 95

# DPI bounds mirror the GUI's number_input limits. The ceiling guards
# against runaway raster allocations (figsize × dpi pixels in RAM).
MIN_DPI = 30
MAX_DPI = 1200


def normalise_format(fmt: str) -> str:
    """Return the canonical extension for ``fmt`` (case-insensitive,
    optional leading dot). Raises ``ValueError`` for unknown formats."""
    key = fmt.strip().lower().lstrip(".")
    if key not in _FORMAT_ALIASES:
        raise ValueError(
            f"unsupported export format {fmt!r} — choose from: " f"{', '.join(SUPPORTED_FORMATS)}"
        )
    return _FORMAT_ALIASES[key]


def save_figure(
    fig: plt.Figure,
    out_path: str | Path,
    *,
    formats: Sequence[str] | None = None,
    dpi: int = 200,
) -> list[Path]:
    """Save ``fig`` to disk in one or more formats.

    If ``formats`` is ``None`` the figure is written exactly to
    ``out_path`` with the format inferred from its suffix — any suffix
    matplotlib supports is passed through unchanged (back-compat with
    the old direct ``fig.savefig`` behaviour, e.g. ``.eps`` / ``.tif``).
    Otherwise one file per entry in ``formats`` (which must come from
    ``SUPPORTED_FORMATS``) is written next to ``out_path`` (same
    directory, same stem, format-specific suffix) and ``out_path``'s own
    suffix is ignored. A bare string is accepted as a single format.

    Returns the list of written paths, in ``formats`` order.
    """
    if not MIN_DPI <= dpi <= MAX_DPI:
        raise ValueError(f"dpi must be between {MIN_DPI} and {MAX_DPI}, got {dpi}")
    out_path = Path(out_path)
    if isinstance(formats, str):
        formats = (formats,)
    if formats is None:
        suffix = (out_path.suffix or ".png").lstrip(".").lower()
        # Known formats go through our canonical path (jpg gets white
        # background + quality); anything else defers to matplotlib.
        fmt_list = [_FORMAT_ALIASES.get(suffix, suffix)]
        targets = [out_path]
    else:
        if not formats:
            raise ValueError("formats must contain at least one entry")
        fmt_list = [normalise_format(f) for f in formats]
        # Preserve order but drop duplicates (e.g. "jpg" + "jpeg").
        seen: dict[str, None] = {}
        for f in fmt_list:
            seen.setdefault(f)
        fmt_list = list(seen)
        targets = [out_path.with_suffix(f".{f}") for f in fmt_list]

    out_path.parent.mkdir(parents=True, exist_ok=True)

    written: list[Path] = []
    for fmt, target in zip(fmt_list, targets, strict=True):
        kwargs: dict = {"dpi": dpi, "bbox_inches": "tight"}
        if fmt == "jpg":
            # JPEG can't store alpha — flatten onto white explicitly.
            kwargs["facecolor"] = "white"
            kwargs["pil_kwargs"] = {"quality": _JPEG_QUALITY}
        fig.savefig(target, format=fmt, **kwargs)
        written.append(target)
    return written
