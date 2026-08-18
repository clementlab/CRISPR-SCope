import numpy as np
import pandas as pd
import pandas.testing as pdt

from CRISPRSCope.editing_rate_ci import (
    EditingRateCIConfig,
    _bootstrap_means,
    _bootstrap_stratified_delta,
    _coverage_bin_ids,
    _permutation_hq_minus_all,
    compute_editing_rate_confidence_intervals,
)


class RecordingRNG:
    """Minimal deterministic RNG that records requested resample shapes."""

    def __init__(self):
        self.sizes = []

    def integers(self, low, high, size):
        assert low == 0
        self.sizes.append(size)
        return np.zeros(size, dtype=int)


class PermutationRecordingRNG:
    """Deterministic random-key generator that records permutation shapes."""

    def __init__(self):
        self.sizes = []

    def random(self, size):
        self.sizes.append(size)
        return np.tile(np.arange(size[1], dtype=float), (size[0], 1))


def test_bootstrap_mean_draws_observed_eligible_sample_size_each_time():
    rng = RecordingRNG()

    bootstrap = _bootstrap_means(
        np.array([10.0, 20.0, 30.0]),
        iterations=5,
        rng=rng,
        batch_size=2,
    )

    assert rng.sizes == [(2, 3), (2, 3), (1, 3)]
    assert bootstrap.tolist() == [10.0] * 5


def test_stratified_delta_fixes_hq_and_non_hq_sample_sizes():
    rng = RecordingRNG()

    bootstrap = _bootstrap_stratified_delta(
        hq_values=np.array([10.0, 20.0]),
        non_hq_values=np.array([0.0, 5.0, 10.0]),
        iterations=5,
        rng=rng,
        batch_size=2,
    )

    assert rng.sizes == [(2, 2), (2, 3), (2, 2), (2, 3), (1, 2), (1, 3)]
    assert bootstrap.tolist() == [6.0] * 5


def test_permutation_preserves_hq_size_without_replacement():
    rng = PermutationRecordingRNG()

    null_effects = _permutation_hq_minus_all(
        all_values=np.array([0.0, 10.0, 20.0, 30.0]),
        n_hq=2,
        iterations=5,
        rng=rng,
        batch_size=2,
    )

    assert rng.sizes == [(2, 4), (2, 4), (1, 4)]
    assert null_effects.tolist() == [-10.0] * 5


def test_coverage_bins_are_exact_through_ceiling_then_use_fixed_widths():
    bin_ids = _coverage_bin_ids(
        np.array([0, 5, 10, 11, 13, 15, 16, 20, 21]),
        exact_max_reads=10,
        bin_width_reads=5,
    )

    assert bin_ids.tolist() == [0, 5, 10, 11, 11, 11, 12, 12, 13]
    assert bin_ids[3] == bin_ids[4]
    assert bin_ids[1] != bin_ids[7]


def test_coverage_adjustment_removes_effect_explained_by_coverage():
    index = [f"group{i}" for i in range(100)] + [
        f"rest_low{i}" for i in range(100)
    ] + [f"rest_matched{i}" for i in range(100)]
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [20] * 100 + [1] * 100 + [20] * 100,
            "modPct.ampA": [0.0] * 100 + [100.0] * 100 + [0.0] * 100,
        },
        index=index,
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI"] * 100 + ["LQ_LO"] * 200},
        index=index,
    )

    result = compute_editing_rate_confidence_intervals(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        min_reads_per_amplicon_per_cell=1,
        config=EditingRateCIConfig(
            bootstrap_iterations=200,
            permutation_iterations=500,
            seed=19,
        ),
        n_processes=1,
    ).iloc[0]

    assert result["permutation_p_value"] < 0.01
    assert result["coverage_adjusted_group_minus_non_group_pct"] == 0.0
    assert result["coverage_adjusted_permutation_p_value"] == 1.0
    assert result["coverage_adjusted_group_retained_pct"] == 100.0
    assert result["coverage_adjusted_non_group_retained_pct"] == 50.0


def test_coverage_adjustment_detects_within_bin_group_effect_for_configured_code():
    index = [f"group{i}" for i in range(50)] + [f"rest{i}" for i in range(50)]
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [13] * 50 + [11] * 50,
            "modPct.ampA": [100.0] * 50 + [0.0] * 50,
        },
        index=index,
    )
    quality_scores = pd.DataFrame(
        {"Color": ["LQ_LO"] * 50 + ["HQ_HI"] * 50},
        index=index,
    )

    result = compute_editing_rate_confidence_intervals(
        editing_summary,
        quality_scores,
        ["LQ_LO"],
        min_reads_per_amplicon_per_cell=1,
        config=EditingRateCIConfig(
            bootstrap_iterations=200,
            permutation_iterations=500,
            seed=23,
        ),
        n_processes=1,
    ).iloc[0]

    assert result["coverage_adjusted_group_minus_non_group_pct"] == 100.0
    assert result["coverage_adjusted_permutation_p_value"] < 0.01
    assert result["coverage_adjusted_bh_p_value"] < 0.01
    assert result["coverage_adjusted_mixed_bin_count"] == 1
    assert result["coverage_adjusted_within_bin_coverage_difference_reads"] == 2.0


def test_coverage_adjustment_reports_absent_common_support():
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [20, 20, 1, 1],
            "modPct.ampA": [0.0, 50.0, 50.0, 100.0],
        },
        index=["g1", "g2", "r1", "r2"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI", "HQ_HI", "LQ_LO", "LQ_LO"]},
        index=editing_summary.index,
    )

    result = compute_editing_rate_confidence_intervals(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        1,
        EditingRateCIConfig(bootstrap_iterations=200, permutation_iterations=200),
        n_processes=1,
    ).iloc[0]

    assert result["coverage_adjusted_mixed_bin_count"] == 0
    assert result["coverage_adjusted_group_n_cells"] == 0
    assert np.isnan(result["coverage_adjusted_group_minus_non_group_pct"])
    assert np.isnan(result["coverage_adjusted_permutation_p_value"])
    assert result["coverage_adjusted_status"].startswith("no_mixed_coverage_bins")


def test_adding_excluded_rows_does_not_change_any_interval():
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [10, 10, 10, 10],
            "modPct.ampA": [0.0, 50.0, 100.0, 50.0],
        },
        index=["cell1", "cell2", "cell3", "cell4"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["LQ_HI", "HQ_HI", "HQ_HI", "LQ_HI"]},
        index=editing_summary.index,
    )
    excluded = pd.DataFrame(
        {
            "totCount.ampA": [0, 10, np.nan],
            "modPct.ampA": [100.0, np.nan, 0.0],
        },
        index=["low_coverage", "missing_rate", "missing_count"],
    )
    excluded_quality = pd.DataFrame(
        {"Color": ["HQ_HI", "HQ_HI", "HQ_HI"]},
        index=excluded.index,
    )
    config = EditingRateCIConfig(
        enabled=True,
        bootstrap_iterations=500,
        permutation_iterations=500,
        seed=42,
        batch_size=25,
    )

    original = compute_editing_rate_confidence_intervals(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        min_reads_per_amplicon_per_cell=5,
        config=config,
        n_processes=1,
    )
    extended = compute_editing_rate_confidence_intervals(
        pd.concat([editing_summary, excluded]),
        pd.concat([quality_scores, excluded_quality]),
        ["HQ_HI"],
        min_reads_per_amplicon_per_cell=5,
        config=config,
        n_processes=1,
    )

    pdt.assert_frame_equal(original, extended)
    assert original.loc[0, "valid_all_bootstrap_replicates"] == 500
    assert original.loc[0, "valid_hq_bootstrap_replicates"] == 500
    assert original.loc[0, "valid_delta_bootstrap_replicates"] == 500
    assert original.loc[0, "valid_permutation_replicates"] == 500


def test_delta_is_zero_when_all_eligible_cells_are_high_quality():
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [10, 10, 10],
            "modPct.ampA": [0.0, 50.0, 100.0],
        },
        index=["cell1", "cell2", "cell3"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI", "HQ_HI", "HQ_HI"]},
        index=editing_summary.index,
    )

    result = compute_editing_rate_confidence_intervals(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        min_reads_per_amplicon_per_cell=1,
        config=EditingRateCIConfig(enabled=True, bootstrap_iterations=200),
        n_processes=1,
    ).iloc[0]

    assert result["hq_minus_all_pct"] == 0.0
    assert result["delta_ci_lower_pct"] == 0.0
    assert result["delta_ci_upper_pct"] == 0.0
    assert result["valid_delta_bootstrap_replicates"] == 200
    assert result["valid_permutation_replicates"] == 0
    assert np.isnan(result["permutation_p_value"])
    assert np.isnan(result["bh_adjusted_p_value"])
    assert result["status"] == "insufficient_non_hq_cells"
