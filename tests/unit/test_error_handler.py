"""Unit tests for the @show_error decorator."""

from __future__ import annotations

import pytest

from shannonpore.ui.error_handler import show_error


@pytest.mark.unit
def test_swallows_exception_and_returns_none() -> None:
    @show_error()
    def boom() -> int:
        raise ValueError("nope")

    assert boom() is None


@pytest.mark.unit
def test_reraise_propagates_exception() -> None:
    @show_error(reraise=True)
    def boom() -> int:
        raise ValueError("nope")

    with pytest.raises(ValueError):
        boom()


@pytest.mark.unit
def test_passes_through_successful_call() -> None:
    @show_error()
    def add(a: int, b: int) -> int:
        return a + b

    assert add(2, 3) == 5
