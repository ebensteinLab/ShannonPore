"""Unit tests for entropy-mode constants + ternary helpers."""

from __future__ import annotations

import pytest

from src.constants import (
    ENTROPY_MODE_BISULFITE,
    ENTROPY_MODE_HELP,
    ENTROPY_MODE_LABELS,
    ENTROPY_MODE_TERNARY,
    ENTROPY_MODE_TRUE_MC,
    ENTROPY_MODES,
)
from src.pipelines.ternary_entropy import required_coverage_for_k
from src.pipelines.whole_genome_duckdb_pipeline import stage_a_build_duckdb_table
from src.state import FilePrepState


@pytest.mark.unit
def test_entropy_modes_defined() -> None:
    assert ENTROPY_MODES == (
        ENTROPY_MODE_TRUE_MC,
        ENTROPY_MODE_BISULFITE,
        ENTROPY_MODE_TERNARY,
    )


@pytest.mark.unit
def test_each_mode_has_label_and_help() -> None:
    for mode in ENTROPY_MODES:
        assert mode in ENTROPY_MODE_LABELS
        assert mode in ENTROPY_MODE_HELP
        assert ENTROPY_MODE_LABELS[mode]
        assert ENTROPY_MODE_HELP[mode]


@pytest.mark.unit
def test_default_state_is_true_mc() -> None:
    assert FilePrepState().entropy_mode == ENTROPY_MODE_TRUE_MC


@pytest.mark.unit
def test_required_coverage_grows_as_3_to_the_k() -> None:
    assert required_coverage_for_k(3) == 27
    assert required_coverage_for_k(4) == 81
    assert required_coverage_for_k(5) == 243


@pytest.mark.unit
def test_stage_a_rejects_ternary_mode(tmp_path) -> None:
    # The 2-state stage must refuse ternary; ternary has its own pipeline.
    with pytest.raises(ValueError, match="ternary"):
        stage_a_build_duckdb_table(
            tsv_path=str(tmp_path / "x.tsv"),
            db_path=str(tmp_path / "x.duckdb"),
            table_name="t",
            tmp_dir=str(tmp_path / "tmp"),
            threads=1,
            force=False,
            entropy_mode="ternary",
        )


@pytest.mark.unit
def test_stage_a_rejects_unknown_mode(tmp_path) -> None:
    # Need a non-empty tab-delimited "TSV" so we get past the file checks
    p = tmp_path / "fake.tsv"
    p.write_text("read_id\tchrom\tref_position\tref_strand\tmod_qual\tmod_code\n")
    with pytest.raises(ValueError, match="Unknown entropy_mode"):
        stage_a_build_duckdb_table(
            tsv_path=str(p),
            db_path=str(tmp_path / "x.duckdb"),
            table_name="t",
            tmp_dir=str(tmp_path / "tmp"),
            threads=1,
            force=False,
            entropy_mode="not-a-mode",
        )
