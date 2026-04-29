"""Shared constants: column names, default parameters, color palettes."""

from __future__ import annotations

from typing import Final

# ── Bedgraph schemas ──────────────────────────────────────────────────────
BEDGRAPH_COLS: Final[tuple[str, ...]] = ("chrom", "start", "end", "value")

# ── Modkit defaults ───────────────────────────────────────────────────────
DEFAULT_MODKIT_THREADS: Final[int] = 8
DEFAULT_MODKIT_QUEUE_SIZE: Final[int] = 1000

# ── Plot palette ─────────────────────────────────────────────────────────
PALETTE_CONTROL: Final[str] = "#1f77b4"
PALETTE_TARGET: Final[str] = "#ff7f0e"
PALETTE_MML: Final[str] = "#3498db"
PALETTE_HMC: Final[str] = "#27ae60"
PALETTE_C: Final[str] = "#bdc3c7"

# ── Tab labels ───────────────────────────────────────────────────────────
TAB_FILE_PREP: Final[str] = "1. File Preparation"
TAB_GRAPH_PREP: Final[str] = "2. Graph Preparation"

# ── Entropy modes ────────────────────────────────────────────────────────
# How 5hmC calls are interpreted when computing methylation entropy.
ENTROPY_MODE_TRUE_MC: Final[str] = "true_mc"
ENTROPY_MODE_BISULFITE: Final[str] = "bisulfite"
ENTROPY_MODE_TERNARY: Final[str] = "ternary"

ENTROPY_MODES: Final[tuple[str, ...]] = (
    ENTROPY_MODE_TRUE_MC,
    ENTROPY_MODE_BISULFITE,
    ENTROPY_MODE_TERNARY,
)

ENTROPY_MODE_LABELS: Final[dict[str, str]] = {
    ENTROPY_MODE_TRUE_MC: "True-mC (5hmC → unmodified C)",
    ENTROPY_MODE_BISULFITE: "Bisulfite-equivalent (5hmC counted as 5mC)",
    ENTROPY_MODE_TERNARY: "Ternary (C / 5mC / 5hmC, 3 states)",
}

ENTROPY_MODE_HELP: Final[dict[str, str]] = {
    ENTROPY_MODE_TRUE_MC: (
        "Filter modkit calls to 5mC only; CpGs with high 5hmC probability "
        "are treated as unmodified C. This isolates true 5mC signal and is "
        "the recommended default for most analyses."
    ),
    ENTROPY_MODE_BISULFITE: (
        "Sum 5mC + 5hmC probabilities into a single 'methylated' signal. "
        "Mirrors what bisulfite sequencing reports (BS cannot distinguish "
        "5mC from 5hmC). Use this when comparing nanopore to BS data."
    ),
    ENTROPY_MODE_TERNARY: (
        "Three states per CpG (0=unmod, 1=5mC, 2=5hmC). Entropy is "
        "normalised by k·log2(3) so values stay in [0,1]. Requires HIGH "
        "coverage: 3^k = 27× for 3 CpGs/bin, 81× for 4 CpGs/bin."
    ),
}
