from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pandas.testing as pdt

from CRISPRSCope.editing_rate_ci import (
    EditingRateDepthStabilityConfig,
    _create_depth_stability_figure,
    _depth_stability_amplicon_order,
    _nested_subsample_means,
    compute_editing_rate_depth_stability,
    write_editing_rate_depth_stability_plot,
)


class ChoiceRNG:
    def __init__(self, samples):
        self.samples = iter(samples)
        self.calls = []

    def choice(self, population_size, size, replace, shuffle):
        sample = np.asarray(next(self.samples), dtype=int)
        assert len(sample) == size
        assert len(np.unique(sample)) == size
        assert np.all((0 <= sample) & (sample < population_size))
        self.calls.append((population_size, size, replace, shuffle))
        return sample


def _example_inputs():
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [10, 10, 10, 10, 1, 10],
            "modPct.ampA": [0.0, 20.0, 40.0, 60.0, 100.0, np.nan],
            "totCount.ampB": [10, 10, 10, 10, 10, 10],
            "modPct.ampB": [100.0, 80.0, 60.0, 40.0, 20.0, 0.0],
        },
        index=["cell4", "cell2", "cell6", "cell1", "cell5", "cell3"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["LQ_HI", "HQ_HI", "HQ_HI", "LQ_HI", "HQ_HI", "HQ_HI"]},
        index=editing_summary.index,
    )
    return editing_summary, quality_scores


def test_nested_subsamples_are_without_replacement_and_share_a_trajectory():
    rng = ChoiceRNG(
        [
            [0, 1, 2],
            [3, 2, 1],
        ]
    )

    means = _nested_subsample_means(
        np.array([10.0, 20.0, 30.0, 40.0]),
        sample_sizes=[1, 3],
        iterations=2,
        rng=rng,
    )

    assert means[1].tolist() == [10.0, 40.0]
    assert means[3].tolist() == [20.0, 30.0]
    assert rng.calls == [
        (4, 3, False, True),
        (4, 3, False, True),
    ]


def test_nested_subsamples_support_full_cohort_as_maximum_depth():
    rng = ChoiceRNG(
        [
            [0, 1, 2, 3],
            [3, 2, 1, 0],
        ]
    )

    means = _nested_subsample_means(
        np.array([10.0, 20.0, 30.0, 40.0]),
        sample_sizes=[2, 4],
        iterations=2,
        rng=rng,
    )

    assert means[2].tolist() == [15.0, 35.0]
    assert means[4].tolist() == [25.0, 25.0]
    assert rng.calls == [
        (4, 4, False, True),
        (4, 4, False, True),
    ]


def test_percentages_use_cohort_sizes_and_append_full_reference():
    editing_summary, quality_scores = _example_inputs()
    result = compute_editing_rate_depth_stability(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        min_reads_per_amplicon_per_cell=5,
        config=EditingRateDepthStabilityConfig(
            enabled=True,
            iterations=100,
            percentages=(25.0, 50.0, 90.0),
            seed=7,
        ),
        n_processes=1,
    )
    amp_a = result.loc[result["amplicon"] == "ampA"]

    all_rows = amp_a.loc[amp_a["cohort"] == "all"].set_index("sample_percent")
    hq_rows = amp_a.loc[amp_a["cohort"] == "hq"].set_index("sample_percent")
    assert all_rows["eligible_n_cells"].unique().tolist() == [4]
    assert all_rows["sample_n_cells"].tolist() == [1, 2, 4, 4]
    assert hq_rows["eligible_n_cells"].unique().tolist() == [2]
    assert hq_rows["sample_n_cells"].tolist() == [1, 1, 2, 2]
    assert all_rows.loc[100.0, "is_full_reference"]
    assert all_rows.loc[100.0, "median_deviation_pp"] == 0.0
    assert all_rows.loc[100.0, "valid_subsamples"] == 1

    duplicate_depth_rows = hq_rows.loc[[25.0, 50.0]]
    summary_columns = [
        "subsample_median_pct",
        "subsample_interval_lower_pct",
        "subsample_interval_upper_pct",
    ]
    assert duplicate_depth_rows[summary_columns].nunique().eq(1).all()


def test_relative_metrics_use_inclusive_hq_threshold_and_own_cohort_denominators():
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [10, 10, 10, 10],
            "modPct.ampA": [0.0, 0.0, 2.0, 2.0],
        },
        index=["cell1", "cell2", "cell3", "cell4"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI", "LQ_HI", "HQ_HI", "LQ_HI"]},
        index=editing_summary.index,
    )
    results = compute_editing_rate_depth_stability(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        1,
        EditingRateDepthStabilityConfig(
            iterations=100,
            percentages=(50.0,),
            relative_min_hq_edit_pct=1.0,
        ),
    )

    assert results["hq_full_estimate_pct"].eq(1.0).all()
    assert results["relative_plot_eligible"].all()
    nonreference = results.loc[results["sample_percent"] == 50.0]
    source_to_relative = {
        "median_deviation_pp": "median_relative_deviation_pct",
        "deviation_interval_lower_pp": "relative_deviation_interval_lower_pct",
        "deviation_interval_upper_pp": "relative_deviation_interval_upper_pct",
        "median_abs_deviation_from_full_pp": "median_absolute_relative_deviation_pct",
        "p95_abs_deviation_from_full_pp": "p95_absolute_relative_deviation_pct",
    }
    for source_column, relative_column in source_to_relative.items():
        expected = (
            100.0
            * nonreference[source_column]
            / nonreference["full_estimate_pct"]
        )
        np.testing.assert_allclose(nonreference[relative_column], expected)
    references = results.loc[results["is_full_reference"]]
    relative_columns = [
        "median_relative_deviation_pct",
        "relative_deviation_interval_lower_pct",
        "relative_deviation_interval_upper_pct",
        "median_absolute_relative_deviation_pct",
        "p95_absolute_relative_deviation_pct",
    ]
    assert references[relative_columns].eq(0.0).all().all()


def test_relative_metrics_are_na_below_hq_threshold_and_safe_for_zero_hq_rate():
    editing_summary = pd.DataFrame(
        {
            "totCount.below": [10, 10, 10, 10],
            "modPct.below": [0.0, 0.0, 1.0, 1.0],
            "totCount.zero": [10, 10, 10, 10],
            "modPct.zero": [0.0, 10.0, 0.0, 10.0],
        },
        index=["cell1", "cell2", "cell3", "cell4"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI", "LQ_HI", "HQ_HI", "LQ_HI"]},
        index=editing_summary.index,
    )
    results = compute_editing_rate_depth_stability(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        1,
        EditingRateDepthStabilityConfig(
            iterations=100,
            percentages=(50.0,),
            relative_min_hq_edit_pct=1.0,
        ),
    )

    assert not results["relative_plot_eligible"].any()
    relative_columns = [
        column for column in results.columns if "relative_deviation" in column
    ]
    assert results[relative_columns].isna().all().all()


def test_hq_ordering_takes_precedence_and_missing_hq_is_last():
    reference_rows = [
        ("ampA", "all", 90.0),
        ("ampA", "hq", 5.0),
        ("ampB", "all", 10.0),
        ("ampB", "hq", 10.0),
        ("ampC", "all", 20.0),
        ("ampC", "hq", 5.0),
        ("ampD", "all", 80.0),
        ("ampE", "all", 80.0),
    ]
    results = pd.DataFrame(
        [
            {
                "amplicon": amplicon,
                "cohort": cohort,
                "full_estimate_pct": estimate,
                "is_full_reference": True,
                "status": "full_reference",
            }
            for amplicon, cohort, estimate in reference_rows
        ]
    )

    assert _depth_stability_amplicon_order(results) == [
        "ampB",
        "ampA",
        "ampC",
        "ampD",
        "ampE",
    ]


def test_every_active_facet_has_x_ticks_and_each_row_has_one_axis_label():
    rows = []
    for amplicon_index in range(5):
        for cohort in ("all", "hq"):
            for sample_percent, deviation in ((10.0, 1.0), (100.0, 0.0)):
                rows.append(
                    {
                        "amplicon": f"amp{amplicon_index}",
                        "cohort": cohort,
                        "sample_percent": sample_percent,
                        "median_deviation_pp": deviation,
                        "deviation_interval_lower_pp": deviation - 0.5,
                        "deviation_interval_upper_pp": deviation + 0.5,
                    }
                )
    plotted = pd.DataFrame(rows)
    order = [f"amp{index}" for index in range(5)]
    figure = _create_depth_stability_figure(
        plotted,
        order,
        median_column="median_deviation_pp",
        lower_column="deviation_interval_lower_pp",
        upper_column="deviation_interval_upper_pp",
        title="Test stability",
        y_label="Deviation (percentage points)",
    )
    figure.canvas.draw()
    active_axes = [axis for axis in figure.axes if axis.get_visible()]

    assert len(active_axes) == 5
    assert all(
        any(label.get_visible() and label.get_text() for label in axis.get_xticklabels())
        for axis in active_axes
    )
    assert sum(bool(axis.get_xlabel()) for axis in active_axes) == 2
    assert sum(bool(axis.get_ylabel()) for axis in active_axes) == 2
    assert active_axes[4].get_xlabel() == "Eligible cohort retained (%)"
    plt.close(figure)


def test_results_ignore_excluded_rows_and_input_order():
    editing_summary, quality_scores = _example_inputs()
    config = EditingRateDepthStabilityConfig(
        enabled=True,
        iterations=100,
        percentages=(25.0, 75.0),
        seed=42,
    )
    original = compute_editing_rate_depth_stability(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        5,
        config,
        n_processes=1,
    )

    reordered_summary = editing_summary.sample(frac=1, random_state=9)
    reordered_quality = quality_scores.reindex(reordered_summary.index)
    reordered = compute_editing_rate_depth_stability(
        reordered_summary,
        reordered_quality,
        ["HQ_HI"],
        5,
        config,
        n_processes=1,
    )
    pdt.assert_frame_equal(original, reordered)

    excluded_summary = pd.DataFrame(
        {
            "totCount.ampA": [0],
            "modPct.ampA": [100.0],
            "totCount.ampB": [0],
            "modPct.ampB": [100.0],
        },
        index=["excluded"],
    )
    excluded_quality = pd.DataFrame({"Color": ["HQ_HI"]}, index=["excluded"])
    extended = compute_editing_rate_depth_stability(
        pd.concat([editing_summary, excluded_summary]),
        pd.concat([quality_scores, excluded_quality]),
        ["HQ_HI"],
        5,
        config,
        n_processes=1,
    )
    pdt.assert_frame_equal(original, extended)


def test_serial_parallel_and_repeated_results_are_identical():
    editing_summary, quality_scores = _example_inputs()
    config = EditingRateDepthStabilityConfig(
        enabled=True,
        iterations=100,
        percentages=(25.0, 75.0),
        seed=123,
    )
    kwargs = {
        "editing_summary": editing_summary,
        "quality_scores": quality_scores,
        "high_quality_codes": ["HQ_HI"],
        "min_reads_per_amplicon_per_cell": 1,
        "config": config,
    }

    first = compute_editing_rate_depth_stability(n_processes=1, **kwargs)
    repeated = compute_editing_rate_depth_stability(n_processes=1, **kwargs)
    parallel = compute_editing_rate_depth_stability(n_processes=2, **kwargs)

    pdt.assert_frame_equal(first, repeated)
    pdt.assert_frame_equal(first, parallel)


def test_insufficient_cohorts_are_retained_but_not_plotted(tmp_path):
    editing_summary = pd.DataFrame(
        {"totCount.ampA": [10, 10], "modPct.ampA": [0.0, 100.0]},
        index=["cell1", "cell2"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI", "LQ_HI"]},
        index=editing_summary.index,
    )
    results = compute_editing_rate_depth_stability(
        editing_summary,
        quality_scores,
        ["HQ_HI"],
        1,
        EditingRateDepthStabilityConfig(
            enabled=True,
            iterations=100,
            percentages=(50.0,),
        ),
    )

    hq_rows = results.loc[results["cohort"] == "hq"]
    assert hq_rows["status"].eq("insufficient_cells").all()
    assert hq_rows["subsample_median_pct"].isna().all()
    assert not results["relative_plot_eligible"].any()
    assert results["median_relative_deviation_pct"].isna().all()

    metadata = write_editing_rate_depth_stability_plot(results, str(tmp_path / "run"))
    assert len(metadata) == 1
    assert Path(str(tmp_path / "run") + ".12_EditingRateDepthStability.png").is_file()
    assert Path(str(tmp_path / "run") + ".12_EditingRateDepthStability.pdf").is_file()
    assert not Path(
        str(tmp_path / "run") + ".13_EditingRateRelativeDepthStability.png"
    ).exists()
    assert not Path(
        str(tmp_path / "run") + ".13_EditingRateRelativeDepthStability.pdf"
    ).exists()
