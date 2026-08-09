"""Unified error handling decorator for Streamlit callbacks.

Replaces the 44 bare `except Exception:` clauses in v3. The `@show_error`
decorator logs the full traceback (via stdlib logging) and surfaces a
short, user-friendly message in the Streamlit UI without crashing the
session.
"""

from __future__ import annotations

import functools
import logging
import traceback
from collections.abc import Callable
from typing import Any, TypeVar

logger = logging.getLogger("shannonpore")

F = TypeVar("F", bound=Callable[..., Any])


def show_error(*, user_message: str | None = None, reraise: bool = False) -> Callable[[F], F]:
    """Decorator: log full traceback, render st.error with friendly message.

    Args:
        user_message: short message shown to user. If None, the exception's
            class name is used.
        reraise: if True, re-raise after rendering. Useful for tests.
    """

    def decorate(fn: F) -> F:
        @functools.wraps(fn)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            try:
                return fn(*args, **kwargs)
            except Exception as exc:
                tb = traceback.format_exc()
                logger.error("Unhandled exception in %s: %s\n%s", fn.__name__, exc, tb)
                msg = user_message or f"{type(exc).__name__}: {exc}"
                try:
                    import streamlit as st

                    st.error(msg)
                    with st.expander("Traceback (for debugging)", expanded=False):
                        st.code(tb, language="text")
                except ImportError:
                    pass
                if reraise:
                    raise
                return None

        return wrapper  # type: ignore[return-value]

    return decorate
