"""shannonpore v4 — Streamlit entrypoint.

Slim entrypoint. All feature logic lives in `src/tabs/`, all state in
`src/state.py`, all plots in `src/plots/`, all pipelines in
`src/pipelines/`.

Run with:
    streamlit run app.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

# Ensure project root is importable as `src.*` when run via `streamlit run app.py`.
_THIS = Path(__file__).resolve().parent
if str(_THIS) not in sys.path:
    sys.path.insert(0, str(_THIS))

import streamlit as st  # noqa: E402

from src import __version__  # noqa: E402
from src.config import REFERENCE_DIR, RESULTS_DIR, ensure_dirs  # noqa: E402
from src.constants import TAB_FILE_PREP, TAB_GRAPH_PREP  # noqa: E402
from src.state import get_state, reset_state  # noqa: E402
from src.tabs import tab_file_prep, tab_graph_prep  # noqa: E402
from src.ui.style import header_band, inject  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)

# Browser tab favicon: prefer the bundled PNG; fall back to a glyph if
# it's not on disk (e.g. user deleted the asset).
_ICON_PATH = _THIS / "Shannonpore_icon.png"
_PAGE_ICON: object = str(_ICON_PATH) if _ICON_PATH.exists() else "◉"

st.set_page_config(
    page_title=f"shannonpore v{__version__}",
    page_icon=_PAGE_ICON,
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "Get Help": "https://github.com/uribertocchitau/shannonpore",
        "Report a bug": "https://github.com/uribertocchitau/shannonpore/issues",
        "About": (
            f"**shannonpore v{__version__}**  \n"
            "Nanopore methylation entropy analysis.  \n"
            "Ebenstein Lab · Tel Aviv University."
        ),
    },
)

inject()


def _sidebar() -> None:
    state = get_state()
    with st.sidebar:
        # ── Branding ────────────────────────────────────────────────
        if _ICON_PATH.exists():
            st.image(str(_ICON_PATH), width=120)
        st.markdown(
            f"""
            <div style="padding:0.4rem 0 0.6rem;">
              <div style="font-family:'IBM Plex Mono',monospace;
                          font-size:0.82rem;letter-spacing:0.18em;
                          text-transform:uppercase;color:#0f4c75;
                          font-weight:600;">shannonpore</div>
              <div style="font-family:'IBM Plex Mono',monospace;
                          font-size:0.72rem;color:#7a7a78;
                          margin-top:0.2rem;">
                v{__version__} · methylation entropy
              </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        st.divider()

        # ── Environment ────────────────────────────────────────────
        st.markdown("**Environment**")
        ref_ok = "✅" if REFERENCE_DIR.exists() else "❌"
        st.caption(f"{ref_ok} `ref_dir`")
        st.code(str(REFERENCE_DIR), language="text")
        st.caption("✅ `results`")
        st.code(str(RESULTS_DIR), language="text")

        if not REFERENCE_DIR.exists():
            st.warning(
                "REFERENCE_DIR does not exist. Set `SHANNONPORE_REF_DIR` "
                "to your reference files directory.",
                icon="⚠️",
            )

        st.divider()

        # ── Session controls ───────────────────────────────────────
        if st.button("Reset state", use_container_width=True):
            reset_state()
            st.rerun()

        if state.last_error:
            st.error(f"Last error: {state.last_error}", icon="❗")

        # ── Footer (in normal flow — NOT absolutely-positioned, so it
        #    can never overlap session widgets above it) ────────────
        st.divider()
        st.markdown(
            """
            <div style="font-family:'IBM Plex Mono',monospace;
                        font-size:0.7rem;color:#9a9a98;
                        letter-spacing:0.04em;line-height:1.7;
                        padding-top:0.4rem;padding-bottom:1rem;">
              ebenstein lab<br/>
              tel aviv university<br/>
              <a href="https://github.com/uribertocchitau/shannonpore"
                 style="color:#9a9a98;text-decoration:none;">
                github.com/uribertocchitau/shannonpore
              </a>
            </div>
            """,
            unsafe_allow_html=True,
        )


def main() -> None:
    ensure_dirs()
    _sidebar()
    header_band(__version__, kind="GUI")

    st.markdown("# Nanopore methylation entropy")
    st.markdown(
        '<div class="caption">Per-CpG entropy (ME), '
        "mean methylation level (MML), and coverage from modkit-extracted "
        "nanopore reads — three modes for handling 5hmC.</div>",
        unsafe_allow_html=True,
    )

    tabs = st.tabs([TAB_FILE_PREP, TAB_GRAPH_PREP])
    with tabs[0]:
        tab_file_prep.render()
    with tabs[1]:
        tab_graph_prep.render()


if __name__ == "__main__":
    main()
