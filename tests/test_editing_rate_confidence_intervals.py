from pathlib import Path

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from CRISPRSCope.editing_rate_ci import (
    EditingRateCIConfig,
    compute_editing_rate_confidence_intervals,
    write_editing_rate_ci_plots,
)


def _example_inputs():
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [10, 10, 10, 1, 10],
            "modPct.ampA": [0.0, 50.0, 100.0, 100.0, "NA"],
            "totCount.ampB": [12, 12, 12, 12, 12],
            "modPct.ampB": [100.0, 50.0, 0.0, 50.0, 100.0],
        },
        index=["cell1", "cell2", "cell3", "cell4", "cell5"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["LQ_HI", "HQ_HI", "HQ_HI", "HQ_LO"]},
        index=["cell1", "cell2", "cell3", "cell4"],
    )
    return editing_summary, quality_scores


def test_compute_uses_cell_weighted_rates_and_configured_hq_codes():
    editing_summary, quality_scores = _example_inputs()
    config = EditingRateCIConfig(enabled=True, bootstrap_iterations=500, seed=7)

    result = compute_editing_rate_confidence_intervals(
        editing_summary=editing_summary,
        quality_scores=quality_scores,
        high_quality_codes=["HQ_HI"],
        min_reads_per_amplicon_per_cell=5,
        config=config,
        n_processes=1,
    ).set_index("amplicon")

    assert result.loc["ampA", "all_n_cells"] == 3
    assert result.loc["ampA", "hq_n_cells"] == 2
    assert result.loc["ampA", "all_estimate_pct"] == 50.0
    assert result.loc["ampA", "hq_estimate_pct"] == 75.0
    assert result.loc["ampA", "hq_minus_all_pct"] == 25.0
    assert result.loc["ampA", "status"] == "ok"
    assert 0 <= result.loc["ampA", "all_ci_lower_pct"] <= result.loc["ampA", "all_ci_upper_pct"] <= 100


def test_parallel_and_serial_results_are_identical():
    editing_summary, quality_scores = _example_inputs()
    config = EditingRateCIConfig(enabled=True, bootstrap_iterations=300, seed=42, batch_size=25)
    kwargs = {
        "editing_summary": editing_summary,
        "quality_scores": quality_scores,
        "high_quality_codes": ["HQ_HI", "HQ_LO"],
        "min_reads_per_amplicon_per_cell": 1,
        "config": config,
    }

    serial = compute_editing_rate_confidence_intervals(n_processes=1, **kwargs)
    parallel = compute_editing_rate_confidence_intervals(n_processes=2, **kwargs)

    pdt.assert_frame_equal(serial, parallel)


def test_insufficient_hq_cells_retains_estimate_but_not_interval():
    editing_summary = pd.DataFrame(
        {"totCount.ampA": [10, 10], "modPct.ampA": [0.0, 100.0]},
        index=["cell1", "cell2"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI", "LQ_HI"]},
        index=editing_summary.index,
    )

    result = compute_editing_rate_confidence_intervals(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        1,
        EditingRateCIConfig(enabled=True, bootstrap_iterations=200),
        n_processes=1,
    ).iloc[0]

    assert result["hq_n_cells"] == 1
    assert result["hq_estimate_pct"] == 0.0
    assert np.isnan(result["hq_ci_lower_pct"])
    assert np.isnan(result["delta_ci_lower_pct"])
    assert result["status"] == "insufficient_hq_cells"


def test_missing_matching_count_column_is_rejected():
    editing_summary = pd.DataFrame({"modPct.ampA": [0.0, 100.0]}, index=["a", "b"])
    quality_scores = pd.DataFrame({"Color": ["HQ_HI", "HQ_HI"]}, index=["a", "b"])

    with pytest.raises(ValueError, match="Missing matching count column"):
        compute_editing_rate_confidence_intervals(
            editing_summary,
            quality_scores,
            ["HQ_HI"],
            1,
            EditingRateCIConfig(enabled=True, bootstrap_iterations=200),
        )


def test_plot_writer_creates_primary_and_delta_artifacts(tmp_path):
    editing_summary, quality_scores = _example_inputs()
    results = compute_editing_rate_confidence_intervals(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        1,
        EditingRateCIConfig(enabled=True, bootstrap_iterations=200),
        n_processes=1,
    )

    metadata = write_editing_rate_ci_plots(results, str(tmp_path / "run"))

    assert len(metadata) == 2
    for suffix in [".10_EditingRateConfidenceIntervals", ".11_EditingRateQualityDelta"]:
        assert Path(str(tmp_path / "run") + suffix + ".png").is_file()
        assert Path(str(tmp_path / "run") + suffix + ".pdf").is_file()
