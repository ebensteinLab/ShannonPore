"""nanoentropy v4 — Streamlit entrypoint.

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

st.set_page_config(
    page_title=f"nanoentropy v{__version__}",
    page_icon="◉",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        "Get Help": "https://github.com/uribertocchitau/nanoentropy",
        "Report a bug": "https://github.com/uribertocchitau/nanoentropy/issues",
        "About": (
            f"**nanoentropy v{__version__}**  \n"
            "Nanopore methylation entropy analysis.  \n"
            "Ebenstein Lab · Tel Aviv University."
        ),
    },
)

inject()


def _sidebar() -> None:
    state = get_state()
    with st.sidebar:
        st.markdown(
            f"""
            <div style="padding:0.4rem 0 1.2rem;">
              <div style="font-family:'IBM Plex Mono',monospace;
                          font-size:0.78rem;letter-spacing:0.18em;
                          text-transform:uppercase;color:#0f4c75;
                          font-weight:600;">nanoentropy</div>
              <div style="font-family:'IBM Plex Mono',monospace;
                          font-size:0.72rem;color:#7a7a78;
                          margin-top:0.15rem;">v{__version__} · methylation entropy</div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        st.markdown("### environment")
        ref_ok = "✓" if REFERENCE_DIR.exists() else "✗"
        st.markdown(
            f"""
            <div style="font-family:'IBM Plex Mono',monospace;font-size:0.78rem;
                        line-height:1.5;color:#1a1a1a;">
              <div><span style="color:#7a7a78;">ref_dir</span> {ref_ok}</div>
              <div style="color:#7a7a78;font-size:0.7rem;
                          margin-bottom:0.4rem;word-break:break-all;">
                {REFERENCE_DIR}
              </div>
              <div><span style="color:#7a7a78;">results</span> ✓</div>
              <div style="color:#7a7a78;font-size:0.7rem;word-break:break-all;">
                {RESULTS_DIR}
              </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        if not REFERENCE_DIR.exists():
            st.warning(
                "REFERENCE_DIR does not exist. Set `NANOENTROPY_REF_DIR` "
                "to your reference files directory."
            )

        st.markdown("### session")
        if st.button("reset state", use_container_width=True):
            reset_state()
            st.rerun()

        st.markdown(
            """
            <div style="position:absolute;bottom:1rem;left:1rem;right:1rem;
                        font-family:'IBM Plex Mono',monospace;font-size:0.66rem;
                        color:#9a9a98;letter-spacing:0.04em;line-height:1.6;">
              <div>ebenstein lab</div>
              <div>tel aviv university</div>
              <div style="margin-top:0.3rem;">
                <a href="https://github.com/uribertocchitau/nanoentropy"
                   style="color:#9a9a98;text-decoration:none;">
                  github.com/uribertocchitau/nanoentropy</a>
              </div>
            </div>
            """,
            unsafe_allow_html=True,
        )

        if state.last_error:
            st.error(f"last error: {state.last_error}")


def main() -> None:
    ensure_dirs()
    _sidebar()
    header_band(__version__, kind="GUI")

    st.markdown("# Nanopore methylation entropy")
    st.markdown(
        '<div class="caption">Per-CpG entropy (ME), '
        'mean methylation level (MML), and coverage from modkit-extracted '
        'nanopore reads — three modes for handling 5hmC.</div>',
        unsafe_allow_html=True,
    )

    tabs = st.tabs([TAB_FILE_PREP, TAB_GRAPH_PREP])
    with tabs[0]:
        tab_file_prep.render()
    with tabs[1]:
        tab_graph_prep.render()


if __name__ == "__main__":
    main()
