from pathlib import Path

import numpy as np
import pandas as pd
import pandas.testing as pdt
import pytest

from CRISPRSCope.editing_rate_ci import (
    EditingRateCIConfig,
    _benjamini_hochberg,
    _significant_plot_rows,
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

    assert result.loc["ampA", "all_cells_n_cells"] == 3
    assert result.loc["ampA", "in_group_n_cells"] == 2
    assert result.loc["ampA", "all_cells_estimate_pct"] == 50.0
    assert result.loc["ampA", "in_group_estimate_pct"] == 75.0
    assert result.loc["ampA", "in_group_minus_all_cells_pct"] == 25.0
    assert result.loc["ampA", "out_group_estimate_pct"] == 0.0
    assert result.loc["ampA", "out_group_n_cells"] == 1
    assert result.loc["ampA", "in_group_minus_out_group_pct"] == 75.0
    assert not any("hq" in column or "non_hq" in column for column in result.columns)
    assert np.isnan(result.loc["ampA", "permutation_p_value"])
    assert result.loc["ampA", "status"] == "insufficient_out_group_cells"
    assert 0 <= result.loc["ampA", "all_cells_ci_lower_pct"] <= result.loc["ampA", "all_cells_ci_upper_pct"] <= 100


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

    assert result["in_group_n_cells"] == 1
    assert result["in_group_estimate_pct"] == 0.0
    assert np.isnan(result["in_group_ci_lower_pct"])
    assert np.isnan(result["in_group_minus_all_cells_ci_lower_pct"])
    assert result["status"] == "insufficient_in_group_cells;insufficient_out_group_cells"


def test_two_sided_permutation_detects_separation_and_is_sign_symmetric():
    index = [f"cell{i}" for i in range(20)]
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI"] * 10 + ["LQ_HI"] * 10},
        index=index,
    )
    config = EditingRateCIConfig(
        enabled=True,
        bootstrap_iterations=200,
        permutation_iterations=2_000,
        seed=13,
    )

    def compute(hq_value, non_hq_value):
        editing_summary = pd.DataFrame(
            {
                "totCount.ampA": [10] * 20,
                "modPct.ampA": [hq_value] * 10 + [non_hq_value] * 10,
            },
            index=index,
        )
        return compute_editing_rate_confidence_intervals(
            editing_summary,
            quality_scores,
            ["HQ_HI"],
            1,
            config,
            n_processes=1,
        ).iloc[0]

    positive = compute(100.0, 0.0)
    negative = compute(0.0, 100.0)

    assert positive["permutation_p_value"] < 0.01
    assert positive["bh_adjusted_p_value"] == positive["permutation_p_value"]
    assert negative["permutation_p_value"] == positive["permutation_p_value"]
    assert positive["valid_permutation_replicates"] == 2_000


def test_identical_groups_have_permutation_p_value_one():
    index = [f"cell{i}" for i in range(8)]
    editing_summary = pd.DataFrame(
        {"totCount.ampA": [10] * 8, "modPct.ampA": [50.0] * 8},
        index=index,
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI"] * 4 + ["LQ_HI"] * 4},
        index=index,
    )

    result = compute_editing_rate_confidence_intervals(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        1,
        EditingRateCIConfig(
            enabled=True,
            bootstrap_iterations=200,
            permutation_iterations=200,
        ),
        n_processes=1,
    ).iloc[0]

    assert result["permutation_p_value"] == 1.0
    assert result["bh_adjusted_p_value"] == 1.0


def test_benjamini_hochberg_adjusts_only_finite_p_values():
    adjusted = _benjamini_hochberg([0.01, 0.04, 0.03, np.nan])

    assert adjusted[:3].tolist() == pytest.approx([0.03, 0.04, 0.04])
    assert np.isnan(adjusted[3])


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
    results.loc[:, "bh_adjusted_p_value"] = [0.01, 0.20]
    results.loc[:, "coverage_adjusted_bh_p_value"] = [0.01, 0.20]

    significant = _significant_plot_rows(results)
    assert significant["amplicon"].tolist() == ["ampA"]

    metadata = write_editing_rate_ci_plots(results, str(tmp_path / "run"))

    assert len(metadata) == 3
    for suffix in [
        ".10_EditingRateConfidenceIntervals",
        ".11_EditingRateQualityDelta",
        ".14_EditingRateCoverageAdjustedEffects",
    ]:
        assert Path(str(tmp_path / "run") + suffix + ".png").is_file()
        assert Path(str(tmp_path / "run") + suffix + ".pdf").is_file()


def test_plot_writer_skips_artifacts_when_no_amplicon_is_significant(tmp_path, caplog):
    caplog.set_level("INFO")
    editing_summary, quality_scores = _example_inputs()
    results = compute_editing_rate_confidence_intervals(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        1,
        EditingRateCIConfig(
            enabled=True,
            bootstrap_iterations=200,
            permutation_iterations=200,
        ),
        n_processes=1,
    )

    stale_paths = []
    for suffix in [
        ".10_EditingRateConfidenceIntervals",
        ".11_EditingRateQualityDelta",
        ".14_EditingRateCoverageAdjustedEffects",
    ]:
        for extension in [".png", ".pdf"]:
            path = Path(str(tmp_path / "run") + suffix + extension)
            path.write_text("stale")
            stale_paths.append(path)

    metadata = write_editing_rate_ci_plots(results, str(tmp_path / "run"))

    assert len(metadata) == 1
    assert metadata[0]["plot_name"].endswith(".14_EditingRateCoverageAdjustedEffects")
    assert "No amplicons passed the coverage-adjusted" in caplog.text
    assert all(
        not path.exists()
        for path in stale_paths
        if not str(path).endswith(
            (".14_EditingRateCoverageAdjustedEffects.png", ".14_EditingRateCoverageAdjustedEffects.pdf")
        )
    )
    assert Path(str(tmp_path / "run") + ".14_EditingRateCoverageAdjustedEffects.png").is_file()
