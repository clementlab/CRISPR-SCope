import numpy as np
import pandas as pd
import pandas.testing as pdt

from CRISPRSCope.editing_rate_ci import (
    EditingRateCIConfig,
    _observed_centered_permutation_differences,
    _observed_centered_swarm_coordinates,
    compute_editing_rate_resampling_analyses,
    write_editing_rate_observed_centered_permutation_swarm_plot,
    write_editing_rate_unconditional_permutation_plot,
)


def _input_data(values, selected_count, *, count=10):
    index = [f"cell{index}" for index in range(len(values))]
    return (
        pd.DataFrame(
            {"totCount.ampA": [count] * len(values), "modPct.ampA": values},
            index=index,
        ),
        pd.DataFrame(
            {"Color": ["HQ_HI"] * selected_count + ["LQ_LO"] * (len(values) - selected_count)},
            index=index,
        ),
    )


def _compute(editing_summary, quality_scores, *, iterations=200, seed=17, n_processes=1):
    return compute_editing_rate_resampling_analyses(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        min_reads_per_amplicon_per_cell=1,
        config=EditingRateCIConfig(
            bootstrap_iterations=iterations,
            permutation_iterations=iterations,
            seed=seed,
        ),
        n_processes=n_processes,
    )


def test_unconditional_distribution_matches_reported_test_and_selected_group_size():
    editing_summary, quality_scores = _input_data([100.0] * 20 + [0.0] * 80, 20)

    results, summaries, simulations = _compute(editing_summary, quality_scores)

    result = results.iloc[0]
    summary = summaries.iloc[0]
    assert summary["all_n_cells"] == 100
    assert summary["selected_group_n_cells"] == 20
    assert summary["selected_group_estimate_pct"] == 100.0
    assert summary["selected_group_minus_all_pct"] == 80.0
    assert summary["permutation_p_value"] == result["permutation_p_value"]
    assert summary["bh_adjusted_p_value"] == result["bh_adjusted_p_value"]
    assert summary["permutation_p_value"] < 0.01
    assert simulations.shape[0] == 200
    assert simulations["selected_group_n_cells"].eq(20).all()
    assert simulations["eligible_all_n_cells"].eq(100).all()
    permutation_means = simulations["permuted_selected_estimate_pct"].to_numpy()
    extreme_count = np.count_nonzero(
        np.abs(permutation_means - summary["all_estimate_pct"])
        >= abs(summary["selected_group_minus_all_pct"])
    )
    assert summary["extreme_permutation_count"] == extreme_count
    assert summary["permutation_p_value"] == (extreme_count + 1) / 201


def test_unconditional_distribution_samples_without_replacement():
    editing_summary, quality_scores = _input_data([0.0, 1.0, 2.0, 4.0], 2)

    _, _, simulations = _compute(editing_summary, quality_scores, iterations=500)

    possible_without_replacement_means = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0}
    assert set(simulations["permuted_selected_estimate_pct"]).issubset(
        possible_without_replacement_means
    )


def test_unconditional_distribution_handles_invariant_and_insufficient_groups():
    invariant_summary, invariant_scores = _input_data([50.0] * 8, 4)
    results, summaries, simulations = _compute(invariant_summary, invariant_scores)
    assert results.iloc[0]["permutation_p_value"] == 1.0
    assert summaries.iloc[0]["status"] == "invariant_permutation_distribution"
    assert simulations["permuted_selected_estimate_pct"].eq(50.0).all()

    one_cell_summary, one_cell_scores = _input_data([100.0, 0.0, 0.0], 1)
    results, summaries, simulations = _compute(one_cell_summary, one_cell_scores)
    assert np.isnan(results.iloc[0]["permutation_p_value"])
    assert summaries.iloc[0]["status"] == "insufficient_selected_group_cells"
    assert summaries.iloc[0]["valid_permutations"] == 0
    assert simulations.empty

    missing_summary, missing_scores = _input_data([0.0, 100.0, 0.0], 2)
    missing_summary.loc[["cell0", "cell1"], "modPct.ampA"] = np.nan
    _, summaries, simulations = _compute(missing_summary, missing_scores)
    assert "insufficient_selected_group_cells" in summaries.iloc[0]["status"]
    assert simulations.empty


def test_unconditional_distribution_is_deterministic_in_parallel_and_plots_all(tmp_path):
    index = [f"cell{index}" for index in range(20)]
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [10] * 20,
            "modPct.ampA": [100.0] * 10 + [0.0] * 10,
            "totCount.ampB": [10] * 20,
            "modPct.ampB": [50.0] * 20,
        },
        index=index,
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI"] * 10 + ["LQ_LO"] * 10}, index=index
    )

    serial = _compute(editing_summary, quality_scores, seed=31, n_processes=1)
    parallel = _compute(editing_summary, quality_scores, seed=31, n_processes=2)
    for serial_frame, parallel_frame in zip(serial, parallel):
        pdt.assert_frame_equal(serial_frame, parallel_frame)

    metadata = write_editing_rate_unconditional_permutation_plot(
        serial[1], serial[2], str(tmp_path / "run")
    )
    assert len(metadata) == 1
    assert (tmp_path / "run.12_EditingRateUnconditionalPermutation.png").is_file()
    assert (tmp_path / "run.12_EditingRateUnconditionalPermutation.pdf").is_file()
    assert set(serial[1]["amplicon"]) == {"ampA", "ampB"}

    assert write_editing_rate_unconditional_permutation_plot(
        serial[1].iloc[0:0], serial[2].iloc[0:0], str(tmp_path / "run")
    ) == []
    assert not (tmp_path / "run.12_EditingRateUnconditionalPermutation.png").exists()
    assert not (tmp_path / "run.12_EditingRateUnconditionalPermutation.pdf").exists()


def test_observed_centered_permutation_differences_filter_and_center_draws():
    summaries = pd.DataFrame(
        {
            "amplicon": ["ampA", "ampB"],
            "selected_group_estimate_pct": [40.0, 60.0],
            "valid_permutations": [2, 0],
            "seed": [42, 42],
        }
    )
    simulations = pd.DataFrame(
        {
            "amplicon": ["ampA", "ampA", "ampB"],
            "permuted_selected_estimate_pct": [35.0, 45.0, 55.0],
        }
    )

    centered = _observed_centered_permutation_differences(summaries, simulations)

    assert centered["amplicon"].tolist() == ["ampA", "ampA"]
    assert centered["observed_in_group_estimate_pct"].tolist() == [40.0, 40.0]
    assert centered["simulated_minus_observed_pct"].tolist() == [-5.0, 5.0]


def test_observed_centered_swarm_is_deterministic_and_writes_artifacts(tmp_path):
    editing_summary, quality_scores = _input_data([100.0] * 10 + [0.0] * 10, 10)
    _, summaries, simulations = _compute(editing_summary, quality_scores, seed=31)
    centered = _observed_centered_permutation_differences(summaries, simulations)

    first_coordinates, first_order = _observed_centered_swarm_coordinates(centered)
    second_coordinates, second_order = _observed_centered_swarm_coordinates(centered)
    pdt.assert_frame_equal(first_coordinates, second_coordinates)
    pdt.assert_frame_equal(first_order, second_order)
    assert first_coordinates["simulated_minus_observed_pct"].eq(
        simulations["permuted_selected_estimate_pct"]
        - summaries.iloc[0]["selected_group_estimate_pct"]
    ).all()

    metadata = write_editing_rate_observed_centered_permutation_swarm_plot(
        summaries, simulations, str(tmp_path / "run")
    )
    assert len(metadata) == 1
    assert (tmp_path / "run.13_EditingRateObservedCenteredPermutationSwarm.png").is_file()
    assert (tmp_path / "run.13_EditingRateObservedCenteredPermutationSwarm.pdf").is_file()


def test_observed_centered_swarm_orders_by_largest_mean_relative_effect():
    centered = pd.DataFrame(
        {
            "amplicon": ["ampA", "ampA", "ampB", "ampB", "ampC", "ampC"],
            "observed_in_group_estimate_pct": [50.0] * 6,
            "simulated_minus_observed_pct": [1.0, 3.0, 2.0, 2.0, -1.0, -3.0],
            "seed": [42] * 6,
        }
    )

    _, order = _observed_centered_swarm_coordinates(centered)

    assert order["amplicon"].tolist() == ["ampA", "ampB", "ampC"]
