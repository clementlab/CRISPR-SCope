"""Editing-rate table, plot, and report-object output adapters."""

import logging

import pandas as pd

from CRISPRSCope.output_artifacts import OutputContext

from .plots_and_report import PlotObject


def write_editing_rate_ci_output(
	output_root,
	cell_quality_to_analyze,
	min_reads_per_amplicon_per_cell,
	config,
	n_processes,
):
	"""Write editing-rate outputs and return plots plus significant amplicons."""
	from CRISPRSCope.editing_rate_ci import (
		_significant_plot_rows,
		compute_editing_rate_resampling_analyses,
		write_editing_rate_ci_plots,
		write_editing_rate_observed_centered_permutation_swarm_plot,
		write_editing_rate_unconditional_permutation_plot,
	)

	outputs = OutputContext(output_root)
	editing_summary_path = outputs.path("editing_summary")
	quality_scores_path = outputs.path("amplicon_score")
	output_path = outputs.path("editing_rate_ci")
	permutation_output_path = outputs.path("editing_rate_unconditional_permutation")
	permutation_simulations_output_path = outputs.path("editing_rate_unconditional_simulations")
	editing_summary = pd.read_csv(editing_summary_path, sep="\t", index_col=0)
	quality_scores = pd.read_csv(quality_scores_path, sep="\t", index_col=0)
	results, permutation_results, permutation_simulations = compute_editing_rate_resampling_analyses(
		editing_summary=editing_summary,
		quality_scores=quality_scores,
		high_quality_codes=cell_quality_to_analyze,
		min_reads_per_amplicon_per_cell=min_reads_per_amplicon_per_cell,
		config=config,
		n_processes=n_processes,
	)
	results.to_csv(output_path, sep="\t", index=False, na_rep="NA", float_format="%.6f")
	logging.info("Wrote editing-rate confidence interval table to %s", output_path)
	permutation_results.to_csv(
		permutation_output_path,
		sep="\t",
		index=False,
		na_rep="NA",
		float_format="%.6f",
	)
	permutation_simulations.to_csv(
		permutation_simulations_output_path,
		sep="\t",
		index=False,
		na_rep="NA",
		float_format="%.6f",
	)
	logging.info("Wrote unconditional permutation summary to %s", permutation_output_path)
	logging.info(
		"Wrote unconditional permutation simulations to %s",
		permutation_simulations_output_path,
	)

	plot_metadata = write_editing_rate_ci_plots(results, output_root)
	plot_metadata.extend(
		write_editing_rate_unconditional_permutation_plot(
			permutation_results,
			permutation_simulations,
			output_root,
		)
	)
	plot_metadata.extend(
		write_editing_rate_observed_centered_permutation_swarm_plot(
			permutation_results,
			permutation_simulations,
			output_root,
		)
	)
	plot_objects = []
	for metadata in plot_metadata:
		artifact_key = metadata["artifact_key"]
		declared = outputs.plot_metadata(artifact_key)
		plot_objects.append(
			PlotObject(
				plot_name=declared["plot_name"],
				plot_title=metadata["plot_title"],
				plot_label=metadata["plot_label"],
				plot_datas=declared["plot_datas"],
			)
		)
	significant_amplicons = _significant_plot_rows(results)["amplicon"].astype(str).tolist()
	return plot_objects, significant_amplicons




def write_editing_rate_depth_stability_output(
	output_root,
	cell_quality_to_analyze,
	min_reads_per_amplicon_per_cell,
	config,
	n_processes,
	significant_amplicons=None,
):
	"""Compute, write, and plot first-pass editing-rate depth stability."""
	from CRISPRSCope.editing_rate_ci import (
		compute_editing_rate_depth_stability,
		write_editing_rate_depth_stability_plot,
	)

	outputs = OutputContext(output_root)
	editing_summary_path = outputs.path("editing_summary")
	quality_scores_path = outputs.path("amplicon_score")
	output_path = outputs.path("editing_rate_depth_stability")
	editing_summary = pd.read_csv(editing_summary_path, sep="\t", index_col=0)
	quality_scores = pd.read_csv(quality_scores_path, sep="\t", index_col=0)
	results = compute_editing_rate_depth_stability(
		editing_summary=editing_summary,
		quality_scores=quality_scores,
		high_quality_codes=cell_quality_to_analyze,
		min_reads_per_amplicon_per_cell=min_reads_per_amplicon_per_cell,
		config=config,
		n_processes=n_processes,
	)
	results.to_csv(output_path, sep="\t", index=False, na_rep="NA", float_format="%.6f")
	logging.info("Wrote editing-rate depth stability table to %s", output_path)

	plot_metadata = write_editing_rate_depth_stability_plot(
		results,
		output_root,
		significant_amplicons=significant_amplicons,
	)
	return [
		PlotObject(
			plot_name=outputs.plot_metadata(metadata["artifact_key"])["plot_name"],
			plot_title=metadata["plot_title"],
			plot_label=metadata["plot_label"],
			plot_datas=outputs.plot_metadata(metadata["artifact_key"])["plot_datas"],
		)
		for metadata in plot_metadata
	]
