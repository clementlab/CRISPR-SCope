import pandas as pd

from CRISPRSCope import cli
from CRISPRSCope.editing_rate_ci import EditingRateCIConfig


def test_cli_writer_creates_table_plots_and_report_links(tmp_path):
    output_root = str(tmp_path / "run")
    editing_summary = pd.DataFrame(
        {
            "totCount.ampA": [10, 10, 10],
            "modPct.ampA": [0.0, 50.0, 100.0],
        },
        index=["cell1", "cell2", "cell3"],
    )
    quality_scores = pd.DataFrame(
        {"Color": ["LQ_HI", "HQ_HI", "HQ_HI"]},
        index=editing_summary.index,
    )
    editing_summary.to_csv(output_root + ".editingSummary.txt", sep="\t")
    quality_scores.to_csv(output_root + ".amplicon_score.txt", sep="\t")

    plot_objects = cli.write_editing_rate_ci_output(
        output_root=output_root,
        cell_quality_to_analyze=["HQ_HI"],
        min_reads_per_amplicon_per_cell=1,
        config=EditingRateCIConfig(enabled=True, bootstrap_iterations=200, seed=42),
        n_processes=1,
    )

    table_path = tmp_path / "run.editingRateConfidenceIntervals.txt"
    assert table_path.is_file()
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
