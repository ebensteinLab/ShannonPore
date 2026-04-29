"""Visual identity for nanoentropy.

A small CSS injection + a header component to replace Streamlit's
default chrome with something that looks like a piece of laboratory
software, not a generic dashboard prototype.

Design notes
------------
- Single primary colour (#0f4c75 — deep prussian) used sparingly for
  emphasis only. Most chrome is neutral grey-on-cream.
- IBM Plex Mono everywhere data is shown (`code`, `pre`, metrics).
- Headers use a low-contrast subtitle line under each H1/H2.
- Metric tiles get a hairline border and aligned baselines.
- The default Streamlit "Made with Streamlit" footer and main-menu
  three-dots are hidden.

This module exposes one entry point: ``inject()``.
"""

from __future__ import annotations

import streamlit as st

# Pinned to a single google-fonts request so the look is the same on
# every machine without bundling fonts.
_FONT_HREF = (
    "https://fonts.googleapis.com/css2?"
    "family=IBM+Plex+Mono:wght@400;500;600&"
    "family=IBM+Plex+Sans:wght@400;500;600;700&display=swap"
)

_CSS = """
<style>
:root {
    --ne-bg:        #fbfbfa;
    --ne-bg-alt:    #f0f0eb;
    --ne-fg:        #1a1a1a;
    --ne-muted:     #7a7a78;
    --ne-line:      #d8d8d2;
    --ne-accent:    #0f4c75;
    --ne-accent-2:  #b86b3a;
    --ne-success:   #2c6e49;
    --ne-warning:   #b8860b;
    --ne-error:     #a04040;
}

html, body, [class*="css"]  {
    font-family: 'IBM Plex Sans', -apple-system, BlinkMacSystemFont, 'Helvetica Neue', sans-serif !important;
}
code, pre, kbd, samp, .stCode, [data-testid="stCodeBlock"] {
    font-family: 'IBM Plex Mono', 'JetBrains Mono', 'SF Mono', Consolas, monospace !important;
    font-size: 0.86rem !important;
}

/* Hide Streamlit chrome that screams "prototype" */
#MainMenu, footer, [data-testid="stToolbar"] { display: none !important; }
header[data-testid="stHeader"] { background: transparent; }

/* Block container — tighten the default 6rem padding */
.main .block-container {
    padding-top: 1.4rem !important;
    padding-bottom: 4rem !important;
    max-width: 1180px;
}

/* Hairline rule under every H1/H2 */
h1, h2 {
    letter-spacing: -0.015em !important;
    color: var(--ne-fg) !important;
    font-weight: 600 !important;
    margin-bottom: 0.2rem !important;
}
h1 { font-size: 1.65rem !important; padding-bottom: 0.4rem; border-bottom: 1px solid var(--ne-line); }
h2 { font-size: 1.20rem !important; margin-top: 1.6rem !important; }
h3 { font-size: 1.00rem !important; color: var(--ne-muted) !important; text-transform: uppercase; letter-spacing: 0.05em; font-weight: 600 !important; }

/* Sidebar branding */
section[data-testid="stSidebar"] {
    background: var(--ne-bg-alt);
    border-right: 1px solid var(--ne-line);
}
section[data-testid="stSidebar"] .block-container {
    padding-top: 1rem;
}

/* Tabs — flat, no rounded pill nonsense */
[data-baseweb="tab-list"] {
    gap: 0 !important;
    border-bottom: 1px solid var(--ne-line);
    margin-bottom: 1rem;
}
[data-baseweb="tab"] {
    padding: 0.55rem 1.1rem !important;
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.82rem !important;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    color: var(--ne-muted);
    border-bottom: 2px solid transparent !important;
}
[data-baseweb="tab"][aria-selected="true"] {
    color: var(--ne-accent) !important;
    border-bottom-color: var(--ne-accent) !important;
    background: transparent !important;
}

/* Buttons */
.stButton > button {
    border: 1px solid var(--ne-line) !important;
    border-radius: 2px !important;
    background: white !important;
    color: var(--ne-fg) !important;
    font-weight: 500 !important;
    padding: 0.4rem 1rem !important;
    transition: border-color 0.12s, background 0.12s !important;
}
.stButton > button:hover {
    border-color: var(--ne-accent) !important;
    color: var(--ne-accent) !important;
}
.stButton > button[kind="primary"] {
    background: var(--ne-accent) !important;
    color: #fbfbfa !important;
    border-color: var(--ne-accent) !important;
}
.stButton > button[kind="primary"]:hover { background: #0a3a5a !important; }

/* Inputs — kill the bouncy border-radius */
input, textarea, select,
[data-baseweb="input"], [data-baseweb="textarea"], [data-baseweb="select"] > div {
    border-radius: 2px !important;
}

/* Metric tiles */
[data-testid="stMetric"] {
    background: white;
    border: 1px solid var(--ne-line);
    border-radius: 2px;
    padding: 0.85rem 1rem 0.7rem;
}
[data-testid="stMetricLabel"] {
    font-family: 'IBM Plex Mono', monospace !important;
    font-size: 0.72rem !important;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    color: var(--ne-muted) !important;
}
[data-testid="stMetricValue"] {
    font-weight: 600 !important;
    font-feature-settings: "tnum" 1;
}

/* Alert boxes — flatter */
[data-testid="stAlert"] {
    border-radius: 2px !important;
    border-left: 3px solid;
    padding: 0.6rem 0.9rem !important;
}
[data-testid="stAlertContentInfo"]    { border-color: var(--ne-accent); }
[data-testid="stAlertContentSuccess"] { border-color: var(--ne-success); }
[data-testid="stAlertContentWarning"] { border-color: var(--ne-warning); }
[data-testid="stAlertContentError"]   { border-color: var(--ne-error); }

/* Code blocks */
[data-testid="stCodeBlock"] {
    background: var(--ne-bg-alt) !important;
    border: 1px solid var(--ne-line);
    border-radius: 2px;
}

/* Dataframe */
[data-testid="stDataFrame"] {
    border: 1px solid var(--ne-line);
    border-radius: 2px;
}

/* Status / progress */
[data-testid="stStatus"] {
    border: 1px solid var(--ne-line);
    border-radius: 2px;
    background: white;
}

/* Caption muted */
.caption, .stCaption, [data-testid="stCaptionContainer"] {
    color: var(--ne-muted) !important;
    font-size: 0.82rem !important;
}

/* Custom branding band at the very top */
.ne-band {
    display: flex;
    align-items: baseline;
    gap: 0.8rem;
    padding: 0.2rem 0 0.6rem;
    border-bottom: 1px solid var(--ne-line);
    margin-bottom: 0.8rem;
}
.ne-band .ne-mark {
    font-family: 'IBM Plex Mono', monospace;
    font-weight: 600;
    font-size: 0.78rem;
    letter-spacing: 0.18em;
    text-transform: uppercase;
    color: var(--ne-accent);
}
.ne-band .ne-version {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.74rem;
    color: var(--ne-muted);
}
.ne-band .ne-spacer { flex: 1; }
.ne-band .ne-tag {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.72rem;
    color: var(--ne-muted);
    letter-spacing: 0.04em;
}

/* Section dividers */
hr {
    border: none !important;
    border-top: 1px solid var(--ne-line) !important;
    margin: 1.4rem 0 !important;
}
</style>
"""


def inject() -> None:
    """Render the global stylesheet. Idempotent — calling twice is fine."""
    st.markdown(f'<link rel="stylesheet" href="{_FONT_HREF}">', unsafe_allow_html=True)
    st.markdown(_CSS, unsafe_allow_html=True)


def header_band(version: str, *, kind: str = "GUI") -> None:
    """Draw a thin product-identity band above the page title.

    Replaces the generic Streamlit `st.title` chrome with a small
    'NANOENTROPY · v4.0.0 · GUI' marker bar.
    """
    st.markdown(
        f"""
        <div class="ne-band">
          <span class="ne-mark">nanoentropy</span>
          <span class="ne-version">v{version}</span>
          <span class="ne-spacer"></span>
          <span class="ne-tag">{kind} · methylation entropy</span>
        </div>
        """,
        unsafe_allow_html=True,
    )


def section(title: str, subtitle: str | None = None) -> None:
    """Numbered section header used in tabs."""
    st.markdown(f"## {title}")
    if subtitle:
        st.markdown(
            f'<div class="caption">{subtitle}</div>',
            unsafe_allow_html=True,
        )
