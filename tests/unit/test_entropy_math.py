"""Unit tests for entropy math primitives."""

from __future__ import annotations

import math

import numpy as np
import pytest

from shannonpore.pipelines.roi_entropy_pipeline import entropy_vec, vect_to_num
from shannonpore.pipelines.whole_genome_duckdb_pipeline import shannon_entropy_from_counts


@pytest.mark.unit
def test_vect_to_num_zero_vector() -> None:
    assert vect_to_num(np.array([0, 0, 0])) == 0


@pytest.mark.unit
def test_vect_to_num_all_ones_three_bits() -> None:
    assert vect_to_num(np.array([1, 1, 1])) == 7


@pytest.mark.unit
def test_vect_to_num_msb_first() -> None:
    # [1, 0, 0] is binary 100 = 4
    assert vect_to_num(np.array([1, 0, 0])) == 4
    # [0, 1, 0] is binary 010 = 2
    assert vect_to_num(np.array([0, 1, 0])) == 2


@pytest.mark.unit
def test_entropy_vec_uniform_distribution() -> None:
    # 4 unique values uniformly distributed → H = log2(4) = 2; normalised by k=2 -> 1.0
    values = np.array([0, 1, 2, 3])
    assert math.isclose(entropy_vec(values, config_options_num=2), 1.0, abs_tol=1e-9)


@pytest.mark.unit
def test_entropy_vec_single_value_zero_entropy() -> None:
    values = np.array([7, 7, 7, 7])
    assert entropy_vec(values, config_options_num=2) == pytest.approx(0.0)


@pytest.mark.unit
def test_entropy_vec_empty_array_returns_zero() -> None:
    assert entropy_vec(np.array([]), config_options_num=2) == 0.0


@pytest.mark.unit
def test_entropy_vec_zero_options_returns_zero() -> None:
    assert entropy_vec(np.array([1, 2, 3]), config_options_num=0) == 0.0


@pytest.mark.unit
def test_shannon_entropy_from_counts_uniform_two_bins() -> None:
    # H([1,1]) = 1
    assert shannon_entropy_from_counts(np.array([5, 5])) == pytest.approx(1.0)


@pytest.mark.unit
def test_shannon_entropy_from_counts_zero_total() -> None:
    assert shannon_entropy_from_counts(np.array([0, 0, 0])) == 0.0


@pytest.mark.unit
def test_shannon_entropy_skips_zero_counts() -> None:
    # Adding zero-count bins should not change the result
    h1 = shannon_entropy_from_counts(np.array([1, 1]))
    h2 = shannon_entropy_from_counts(np.array([1, 1, 0, 0]))
    assert h1 == pytest.approx(h2)
