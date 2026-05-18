import matplotlib.pyplot as plt
import pandas as pd

from CRISPRSCope.cli import (
	_numeric_mod_pct_columns,
	_numeric_tot_count_columns,
	amp_per_cell_filtered,
	cell_per_amp_filtered,
	generate_edit_histogram,
	generate_upset_plot,
	log_log_plot,
	mod_per_amp_filtered,
)


def _parsed_information_with_missing_amplicon():
	return pd.DataFrame(
		{
			"totCount.amp_ok": [12, 0],
			"modPct.amp_ok": [25.0, 0],
			"totCount.amp_failed": ["NA", "NA"],
			"modPct.amp_failed": ["NA", "NA"],
			"Color": ["HQ_HI", "HQ_LO"],
		},
		index=["cellA", "cellB"],
	)


def test_numeric_plot_helpers_coerce_missing_amplicon_values():
	parsed_information = _parsed_information_with_missing_amplicon()

	totals = _numeric_tot_count_columns(parsed_information)
	mods = _numeric_mod_pct_columns(parsed_information)

	assert totals.loc["cellA", "totCount.amp_ok"] == 12
	assert pd.isna(totals.loc["cellA", "totCount.amp_failed"])
	assert totals.sum(axis=1).to_dict() == {"cellA": 12, "cellB": 0}
	assert (totals >= 1).sum(axis=1).to_dict() == {"cellA": 1, "cellB": 0}
	assert (mods > 0).sum(axis=1).to_dict() == {"cellA": 1, "cellB": 0}


def test_in_memory_plots_handle_missing_amplicon_values(tmp_path, monkeypatch):
	monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: None)
	parsed_information = _parsed_information_with_missing_amplicon()
	output_root = str(tmp_path / "run")
	cell_quality_to_analyze = ["HQ_HI", "HQ_LO"]

	log_log_plot(parsed_information.copy(), output_root, cell_quality_to_analyze, filtered=False)
	cell_per_amp_filtered(parsed_information.copy(), output_root, cell_quality_to_analyze)
	amp_per_cell_filtered(parsed_information.copy(), output_root, cell_quality_to_analyze)
	mod_per_amp_filtered(parsed_information.copy(), output_root, cell_quality_to_analyze)


def test_file_based_edit_plots_treat_missing_mod_pct_as_unedited(tmp_path, monkeypatch):
	monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: None)
	output_root = str(tmp_path / "run")
	parsed_information = _parsed_information_with_missing_amplicon()
	parsed_information.index.name = "cell"
	parsed_information.drop(columns=["Color"]).to_csv(output_root + ".filteredEditingSummary.txt", sep="\t")

	amplicon_score = pd.DataFrame(
		{"Color": ["HQ_HI", "HQ_LO"]},
		index=["cellA", "cellB"],
	)
	amplicon_score.to_csv(output_root + ".amplicon_score.txt", sep="\t")

	generate_edit_histogram(output_root, ["HQ_HI", "HQ_LO"])
	generate_upset_plot(output_root, ["HQ_HI", "HQ_LO"])
