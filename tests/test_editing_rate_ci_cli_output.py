import pandas as pd

from CRISPRSCope import cli
from CRISPRSCope.editing_rate_ci import EditingRateCIConfig, EditingRateDepthStabilityConfig


def test_cli_writer_creates_table_plots_and_report_links(tmp_path):
    output_root = str(tmp_path / "run")
    index = [f"cell{i}" for i in range(20)]
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [10] * 20,
            "modPct.ampA": [100.0] * 10 + [0.0] * 10,
        },
        index=index,
    )
    quality_scores = pd.DataFrame(
        {"Color": ["HQ_HI"] * 10 + ["LQ_HI"] * 10},
        index=index,
    )
    editing_summary.to_csv(output_root + ".editingSummary.txt", sep="\t")
    quality_scores.to_csv(output_root + ".amplicon_score.txt", sep="\t")

    plot_objects = cli.write_editing_rate_ci_output(
        output_root=output_root,
        cell_quality_to_analyze=["HQ_HI"],
        min_reads_per_amplicon_per_cell=1,
        config=EditingRateCIConfig(
            enabled=True,
            bootstrap_iterations=200,
            permutation_iterations=2_000,
            seed=42,
        ),
        n_processes=1,
    )

    table_path = tmp_path / "run.editingRateConfidenceIntervals.txt"
    assert table_path.is_file()
    table = pd.read_csv(table_path, sep="\t")
    assert table.loc[0, "permutation_p_value"] < 0.01
    assert table.loc[0, "bh_adjusted_p_value"] <= 0.05
    assert len(plot_objects) == 2
    assert all(plot.datas == [("Editing-rate confidence intervals", str(table_path))] for plot in plot_objects)

    report_path = tmp_path / "run.html"
    cli.make_report(
        report_file=str(report_path),
        report_name="Dataset Summary Report",
        results_folder="",
        crispresso_run_names=[],
        crispresso_sub_html_files={},
        summary_plot_objects=plot_objects,
    )

    report_text = report_path.read_text()
    assert "Amplicon editing-rate confidence intervals" in report_text
    assert "run.editingRateConfidenceIntervals.txt" in report_text


def test_cli_writer_keeps_table_but_omits_report_plots_without_significance(tmp_path):
    output_root = str(tmp_path / "run")
    index = [f"cell{i}" for i in range(8)]
    pd.DataFrame(
        {"totCount.ampA": [10] * 8, "modPct.ampA": [50.0] * 8},
        index=index,
    ).to_csv(output_root + ".editingSummary.txt", sep="\t")
    pd.DataFrame(
        {"Color": ["HQ_HI"] * 4 + ["LQ_HI"] * 4},
        index=index,
    ).to_csv(output_root + ".amplicon_score.txt", sep="\t")

    plot_objects = cli.write_editing_rate_ci_output(
        output_root=output_root,
        cell_quality_to_analyze=["HQ_HI"],
        min_reads_per_amplicon_per_cell=1,
        config=EditingRateCIConfig(
            enabled=True,
            bootstrap_iterations=200,
            permutation_iterations=200,
        ),
        n_processes=1,
    )

    table_path = tmp_path / "run.editingRateConfidenceIntervals.txt"
    assert table_path.is_file()
    assert pd.read_csv(table_path, sep="\t")["amplicon"].tolist() == ["ampA"]
    assert plot_objects == []
    assert not list(tmp_path.glob("run.10_EditingRateConfidenceIntervals.*"))
    assert not list(tmp_path.glob("run.11_EditingRateQualityDelta.*"))


def test_depth_stability_writer_creates_table_plot_and_report_link(tmp_path):
    output_root = str(tmp_path / "run")
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [10, 10, 10, 10],
            "modPct.ampA": [0.0, 25.0, 75.0, 100.0],
        },
        index=["cell1", "cell2", "cell3", "cell4"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["LQ_HI", "HQ_HI", "HQ_HI", "LQ_HI"]},
        index=editing_summary.index,
    )
    editing_summary.to_csv(output_root + ".editingSummary.txt", sep="\t")
    quality_scores.to_csv(output_root + ".amplicon_score.txt", sep="\t")

    plot_objects = cli.write_editing_rate_depth_stability_output(
        output_root=output_root,
        cell_quality_to_analyze=["HQ_HI"],
        min_reads_per_amplicon_per_cell=1,
        config=EditingRateDepthStabilityConfig(
            enabled=True,
            iterations=100,
            percentages=(25.0, 50.0, 90.0),
            seed=42,
        ),
        n_processes=1,
    )

    table_path = tmp_path / "run.editingRateDepthStability.txt"
    assert table_path.is_file()
    assert (tmp_path / "run.12_EditingRateDepthStability.png").is_file()
    assert (tmp_path / "run.12_EditingRateDepthStability.pdf").is_file()
    assert (tmp_path / "run.13_EditingRateRelativeDepthStability.png").is_file()
    assert (tmp_path / "run.13_EditingRateRelativeDepthStability.pdf").is_file()
    assert len(plot_objects) == 2
    assert all(
        plot.datas == [("Editing-rate depth stability", str(table_path))]
        for plot in plot_objects
    )

    table = pd.read_csv(table_path, sep="\t")
    assert {
        "hq_full_estimate_pct",
        "relative_min_hq_edit_pct",
        "relative_plot_eligible",
        "median_relative_deviation_pct",
        "relative_deviation_interval_lower_pct",
        "relative_deviation_interval_upper_pct",
        "median_absolute_relative_deviation_pct",
        "p95_absolute_relative_deviation_pct",
    }.issubset(table.columns)

    report_path = tmp_path / "run.html"
    cli.make_report(
        report_file=str(report_path),
        report_name="Dataset Summary Report",
        results_folder="",
        crispresso_run_names=[],
        crispresso_sub_html_files={},
        summary_plot_objects=plot_objects,
    )

    report_text = report_path.read_text()
    assert "Editing-rate cell-depth stability" in report_text
    assert "Relative editing-rate cell-depth stability" in report_text
    assert "run.editingRateDepthStability.txt" in report_text
