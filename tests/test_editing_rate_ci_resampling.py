import numpy as np
import pandas as pd
import pandas.testing as pdt

from CRISPRSCope.editing_rate_ci import (
    EditingRateCIConfig,
    _bootstrap_means,
    _bootstrap_stratified_delta,
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
