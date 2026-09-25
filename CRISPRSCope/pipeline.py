"""Shared imports for CRISPRSCope pipeline-stage modules."""
from __future__ import annotations

import errno
import gzip
import hashlib
import json
import logging
import multiprocessing as mp
import os
import random
import re
import shlex
import shutil
import signal
import subprocess as sb
import sys
import threading
import time
import traceback
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from functools import partial
from typing import Optional

import dnaio
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from CRISPResso2 import CRISPRessoShared
from adjustText import adjust_text
from scipy.stats import multinomial
from upsetplot import UpSet

from CRISPRSCope import __version__
from CRISPRSCope.cache import CacheManager, OutputRootLock
from CRISPRSCope.io_utils import open_text_maybe_gzip
from CRISPRSCope.output_artifacts import OutputContext, OutputManifest

from .amplicon_assignment import split_reads_by_amplicon
from .crispresso import (
    filter_amplicon_reads, parse_crispresso_outputs, prune_removed_amplicon_caches,
    run_crispresso_commands,
    write_filtered_editing_summary_from_filtered_crispresso,
)
from .fastq_processing import parse_and_align_reads
from .editing_rate_outputs import (
    write_editing_rate_ci_output,
    write_editing_rate_depth_stability_output,
)
from .paths import validate_output_root
from .plots_and_report import (
    PlotObject, amp_per_cell_filtered, cell_per_amp_filtered, declared_plot_object,
    generate_amplicon_coverage_plot, generate_cell_coverage_plot,
    generate_edit_histogram, generate_read_depth_boxplots, generate_upset_plot,
    log_log_plot, make_report, mod_per_amp_filtered, plot_amp_score,
)
from .settings import (
    _parse_amplicon_score_config, _parse_editing_rate_ci_config,
    _parse_editing_rate_depth_stability_config, _parse_cache_config,
    _parse_settings_file, _resolve_settings_path, parse_settings,
)

def _require_selected_barcodes(parsed_information, cell_quality_to_analyze):
	"""Fail early when configured quality groups contain no barcodes."""
	selected_barcode_count = int(
		parsed_information['Color'].isin(cell_quality_to_analyze).sum()
	)
	if selected_barcode_count == 0:
		raise ValueError(
			"No barcodes match the configured analysis groups "
			f"{cell_quality_to_analyze}. Lower the supported-breadth or barcode-rank "
			"thresholds, or select additional cell-quality groups."
		)
	return selected_barcode_count



def write_h5ad_output(output_root, settings_file, h5ad_output=None, h5ad_export_config=None, n_processes=None):
	"""
	Build and save an h5ad file from pipeline outputs for a single run.
	"""
	if h5ad_output is None:
		h5ad_output = output_root + ".h5ad"

	try:
		from CRISPRSCope.h5ad.api import build_h5ad_from_output_root
	except ModuleNotFoundError as exc:
		if exc.name != "CRISPRSCope":
			raise RuntimeError(
				"h5ad export requires optional dependencies: anndata and pyarrow."
			) from exc
		try:
			from h5ad.api import build_h5ad_from_output_root
		except ImportError as inner_exc:
			raise RuntimeError(
				"h5ad export requires optional dependencies: anndata and pyarrow."
			) from inner_exc
	except ImportError as exc:
		raise RuntimeError(
			"h5ad export requires optional dependencies: anndata and pyarrow."
		) from exc

	build_h5ad_from_output_root(
		output_root=output_root,
		output_path=h5ad_output,
		config=h5ad_export_config,
		settings_path=settings_file,
		n_processes=n_processes,
	)
	return h5ad_output


def _run_pipeline_with_manifest_finalization(manifest_observer=None):
	"""Run and finalize the manifest while the caller still owns the run lock."""
	active_manifest = None

	def observe_manifest(manifest):
		nonlocal active_manifest
		active_manifest = manifest
		if manifest_observer is not None:
			manifest_observer(manifest)

	try:
		result = _run_pipeline_unlocked(observe_manifest)
	except BaseException as error:
		if active_manifest is not None:
			active_manifest.fail(active_manifest.active_stage or "initialization", error)
			try:
				active_manifest.write()
			except BaseException:
				logging.exception("Failed to write output manifest after pipeline failure")
		raise
	else:
		if active_manifest is not None:
			active_manifest.complete()
			active_manifest.write()
		return result


def run_pipeline(manifest_observer=None):
	"""Resolve and lock the output root before entering the pipeline."""
	if len(sys.argv) > 1 and sys.argv[1] in {"--version", "-V"}:
		print(__version__)
		return
	if len(sys.argv) < 2:
		return _run_pipeline_with_manifest_finalization(manifest_observer)
	settings_file = os.path.abspath(sys.argv[1])
	settings = _parse_settings_file(settings_file)
	output_root = settings_file
	if "output_root" in settings:
		output_root = _resolve_settings_path(settings["output_root"], os.path.dirname(settings_file))
	output_root = validate_output_root(output_root)
	with OutputRootLock(output_root):
		return _run_pipeline_with_manifest_finalization(manifest_observer)


def _run_pipeline_unlocked(manifest_observer=None):
	"""
	Top-level pipeline entry point for the CRISPRSCope processing workflow.

	High-level summary
	------------------
	Orchestrates a full end-to-end run that takes sequencing input (FASTQ or an
	aligned BAM), assigns reads to amplicons, runs CRISPResso2 per-amplicon,
	aggregates per-cell / per-amplicon allele statistics, produces summary
	plots and an HTML report, and (optionally) cleans up intermediate files.

	Steps performed
	---------------
	1. Parse configuration / command-line arguments.
	2. Parse raw reads and/or run alignment to produce a name-sorted BAM and
	   a `reads_per_cell` mapping of barcode → read counts.
	3. Build primer/amplicon lookup tables and split reads into per-amplicon
	   FASTQ files (writes an ampliconInfo file to speed re-runs).
	4. Launch CRISPResso2 runs (one per amplicon), possibly in parallel.
	5. Parse CRISPResso outputs and map allele/edit calls back to cells.
	6. Aggregate results, generate plots (coverage, edit histograms, UpSet,
	   etc.), and assemble a final HTML report; optionally zip outputs
	"""
	if len(sys.argv) > 1 and sys.argv[1] in {"--version", "-V"}:
		print(__version__)
		return

	start_settings = time.time()
	(r1, r2, constant1, constant2, allow_barcode_mismatches,barcode_file, amplicon_file,
		primer_lookup_len,adapter_DNA, amp_file_dir, alt_alleles_file, bowtie2_index,
		crispresso_dir, output_root, n_processes, keep_intermediate_files,
		ignore_substitutions, assign_reads_to_all_possible_amplicons, suppress_sub_crispresso_plots,
		min_total_reads_per_barcode, min_reads_per_amplicon_per_cell, cell_quality_to_analyze,
		write_h5ad, h5ad_output, h5ad_export_config, debug_rescued_reads_bam,
		debug_rejected_rescue_reads_bam, debug_require_strict_amplicon_alignment,
		partial_rescue_min_mean_read_quality, write_output_manifest, settings_file
		) = parse_settings(sys.argv)
	amplicon_score_config = _parse_amplicon_score_config(settings_file)
	editing_rate_ci_config = _parse_editing_rate_ci_config(settings_file)
	editing_rate_depth_stability_config = _parse_editing_rate_depth_stability_config(settings_file)
	cache_config = _parse_cache_config(settings_file)
	end_settings = time.time() - start_settings
	#print(f"Parse Settings: {end_settings}")

	output_root = validate_output_root(output_root)
	outputs = OutputContext(output_root, h5ad_output=h5ad_output)
	manifest = OutputManifest(outputs) if write_output_manifest else None
	if manifest is not None:
		manifest.configure_cache(cache_config.mode.value)
	if manifest_observer is not None:
		manifest_observer(manifest)
	cache_manager = CacheManager(
		output_root,
		cache_config,
		producer_version=__version__,
		event_callback=manifest.record_cache_event if manifest is not None else None,
	)

	def record_plot(key, plot_object, reason="plot requirements were not met"):
		if manifest is None:
			return
		if plot_object is None:
			manifest.mark_skipped(key, reason)
		else:
			manifest.mark_written(key)

	def remove_optional_artifacts(keys, reason="stale output from an earlier run"):
		removed = set(outputs.remove_optional(keys))
		if manifest is None:
			return set()
		for key in keys:
			if any(path in removed for path in outputs.paths(key)):
				manifest.mark_removed_stale(key, reason)
		return {key for key in keys if any(path in removed for path in outputs.paths(key))}

	def mark_written(*keys):
		if manifest is not None:
			for key in keys:
				manifest.mark_written(key)

	def mark_skipped_unless_removed(keys, reason, removed_keys):
		if manifest is not None:
			for key in keys:
				if key not in removed_keys:
					manifest.mark_skipped(key, reason)

	if manifest is not None:
		manifest.set_stage("parse_and_align_reads")
	start_parse_and_align = time.time()
	aligned_bam, reads_per_cell = parse_and_align_reads(r1,r2,constant1,constant2,output_root,barcode_file,allow_barcode_mismatches,adapter_DNA,bowtie2_index,n_processes,keep_intermediate_files, cache_manager=cache_manager)
	end_parse_and_align = time.time() - start_parse_and_align
	logging.info(f"Parse and Align Reads: {end_parse_and_align}")


	if manifest is not None:
		manifest.set_stage("split_reads_by_amplicon")
	start_split_reads = time.time()
	amplicon_names, amplicon_information, amplicon_info_file = split_reads_by_amplicon(aligned_bam, output_root, amplicon_file, alt_alleles_file, primer_lookup_len, amp_file_dir, bowtie2_index, adapter_DNA, n_processes, keep_intermediate_files, reads_per_cell, min_total_reads_per_barcode, assign_reads_to_all_possible_amplicons, debug_rescued_reads_bam, debug_require_strict_amplicon_alignment, debug_rejected_rescue_reads_bam, partial_rescue_min_mean_read_quality, cache_manager=cache_manager)
	prune_removed_amplicon_caches(
		cache_manager, amplicon_names, output_root, crispresso_dir,
	)
	end_split_reads = time.time() - start_split_reads
	logging.info(f"Split Reads by Amplicon: {end_split_reads}")
	mark_written(
		"valid_amplicons",
		"aligned_read_counts",
		"unaligned_read_counts",
		"amplicon_classification",
	)

#    print(f"Line266\n{amplicon_names=}\n{amplicon_information=}\n{amplicon_info_file=}")

	if manifest is not None:
		manifest.set_stage("run_crispresso")
	start_crispresso = time.time()
	crispresso_information = run_crispresso_commands(amplicon_names,amplicon_information,output_root,crispresso_dir,suppress_sub_crispresso_plots,n_processes, alleles = False, cache_manager=cache_manager)
	end_crispresso = time.time() - start_crispresso
	logging.info(f"Run CRISPResso: {end_crispresso}")

	if manifest is not None:
		manifest.set_stage("parse_crispresso_outputs")
	start_parse_crispresso = time.time()
	parsed_information = parse_crispresso_outputs(amplicon_names,amplicon_information,amplicon_info_file,
											  crispresso_information,output_root, min_total_reads_per_barcode, min_reads_per_amplicon_per_cell,
											  n_processes=n_processes,
											  ignore_substitutions=ignore_substitutions,
											  amplicon_score_config=amplicon_score_config,
											  cache_manager=cache_manager)
	end_parse_crispresso = time.time() - start_parse_crispresso
	logging.info(f"Parse CRISPResso Outputs: {end_parse_crispresso}")
	mark_written(
		"editing_summary",
		"editing_summary_pseudobulk",
		"amplicon_score",
		"filtered_editing_summary_pseudobulk",
	)

	if manifest is not None:
		manifest.set_stage("filter_amplicon_reads")
	start_filter_amplicon = time.time()
	_require_selected_barcodes(parsed_information, cell_quality_to_analyze)
	filter_amplicon_reads(output_root, parsed_information, amplicon_names, cell_quality_to_analyze, n_processes, cache_manager=cache_manager)
	end_filter_amplicon = time.time() - start_filter_amplicon
	logging.info(f"Filter Amplicon Reads: {end_filter_amplicon}")

	if manifest is not None:
		manifest.set_stage("run_filtered_crispresso")
	start_run_crispresso2 = time.time()
	crispresso_filtered_information = run_crispresso_commands(amplicon_names,amplicon_information,output_root,crispresso_dir,suppress_sub_crispresso_plots,n_processes, alleles = True, cache_manager=cache_manager)
	end_run_crispresso2 = time.time() - start_run_crispresso2
	logging.info(f"Run CRISPResso 2: {end_run_crispresso2}")

	if manifest is not None:
		manifest.set_stage("write_filtered_editing_summary")
	start_filtered_summary = time.time()
	filtered_parsed_information = write_filtered_editing_summary_from_filtered_crispresso(
		amplicon_names=amplicon_names,
		crispresso_information=crispresso_information,
		crispresso_filtered_information=crispresso_filtered_information,
		output_root=output_root,
		ignore_substitutions=ignore_substitutions,
	)
	end_filtered_summary = time.time() - start_filtered_summary
	logging.info(f"Write Filtered Editing Summary: {end_filtered_summary}")
	mark_written("filtered_editing_summary")

	if manifest is not None:
		manifest.set_stage("generate_summary_plots")
	filtered_summary_plot_objects = []

	filtered_read_count_plot_obj = generate_read_depth_boxplots(output_root, cell_quality_to_analyze)
	record_plot("cell_coverage_boxplot", filtered_read_count_plot_obj)
	if filtered_read_count_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_read_count_plot_obj)

	# Cell Coverage Bar Chart
	filtered_cell_cov_plot_obj = generate_cell_coverage_plot(output_root, cell_quality_to_analyze)
	record_plot("cell_coverage_plot", filtered_cell_cov_plot_obj)
	if filtered_cell_cov_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_cell_cov_plot_obj)

	# Avg Read Count Per Amplicon Boxplot
	filtered_amp_cov_plot_obj = generate_amplicon_coverage_plot(output_root, cell_quality_to_analyze)
	record_plot("amplicon_coverage_plot", filtered_amp_cov_plot_obj)
	if filtered_amp_cov_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_amp_cov_plot_obj)

	# Upset plot displaying edit site combinations
	filtered_upset_plot_obj = generate_upset_plot(output_root, cell_quality_to_analyze)
	record_plot("edit_combinations_plot", filtered_upset_plot_obj)
	if filtered_upset_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_upset_plot_obj)

	# Histogram displaying the frequency of edited site number
	filtered_hist_plot_obj = generate_edit_histogram(output_root, cell_quality_to_analyze)
	record_plot("edit_histogram_plot", filtered_hist_plot_obj)
	if filtered_hist_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_hist_plot_obj)

	# Generate a log read count (y) vs log barcode rank (x) colored by cell quality category
	filtered_log_log_plot_obj = log_log_plot(parsed_information, output_root, cell_quality_to_analyze, filtered = False)
	record_plot("log_log_plot", filtered_log_log_plot_obj)
	if filtered_log_log_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_log_log_plot_obj)

	#
	filtered_cell_per_amp_obj = cell_per_amp_filtered(filtered_parsed_information, output_root, cell_quality_to_analyze)
	record_plot("cell_count_per_amplicon_plot", filtered_cell_per_amp_obj)
	if filtered_cell_per_amp_obj is not None:
		filtered_summary_plot_objects.append(filtered_cell_per_amp_obj)

	#
	filtered_amp_per_cell_obj = amp_per_cell_filtered(filtered_parsed_information, output_root, cell_quality_to_analyze)
	record_plot("amplicon_covered_per_cell_plot", filtered_amp_per_cell_obj)
	if filtered_amp_per_cell_obj is not None:
		filtered_summary_plot_objects.append(filtered_amp_per_cell_obj)

	#
	filtered_mod_pct_plot_obj = mod_per_amp_filtered(filtered_parsed_information, output_root, cell_quality_to_analyze)
	record_plot("modification_percentage_plot", filtered_mod_pct_plot_obj)
	if filtered_mod_pct_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_mod_pct_plot_obj)

	#
	amp_score_plot_obj = plot_amp_score(output_root, config=amplicon_score_config)
	record_plot("amplicon_score_plot", amp_score_plot_obj)
	if amp_score_plot_obj is not None:
		filtered_summary_plot_objects.append(amp_score_plot_obj)

	from CRISPRSCope.editing_rate_ci import remove_editing_rate_plot_artifacts
	remove_optional_artifacts(
		(
			"editing_rate_confidence_intervals_plot",
			"editing_rate_coverage_adjusted_effects_plot",
			"editing_rate_unconditional_permutation_plot",
			"editing_rate_observed_centered_permutation_swarm_plot",
			"editing_rate_depth_stability_plot",
		)
	)
	remove_editing_rate_plot_artifacts(output_root)

	editing_rate_significant_amplicons = None
	if editing_rate_ci_config.enabled:
		start_editing_rate_ci = time.time()
		(
			editing_rate_ci_plot_objects,
			editing_rate_significant_amplicons,
		) = write_editing_rate_ci_output(
			output_root=output_root,
			cell_quality_to_analyze=cell_quality_to_analyze,
			min_reads_per_amplicon_per_cell=min_reads_per_amplicon_per_cell,
			config=editing_rate_ci_config,
			n_processes=n_processes,
		)
		filtered_summary_plot_objects.extend(editing_rate_ci_plot_objects)
		mark_written(
			"editing_rate_ci",
			"editing_rate_unconditional_permutation",
			"editing_rate_unconditional_simulations",
		)
		if manifest is not None:
			created_roots = {plot.name for plot in editing_rate_ci_plot_objects}
			for key in (
				"editing_rate_confidence_intervals_plot",
				"editing_rate_coverage_adjusted_effects_plot",
				"editing_rate_unconditional_permutation_plot",
				"editing_rate_observed_centered_permutation_swarm_plot",
			):
				if outputs.plot_root(key) in created_roots:
					manifest.mark_written(key)
				else:
					manifest.mark_skipped(key, "plot requirements were not met")
		logging.info(
			"Generated editing-rate confidence intervals in %.2f seconds",
			time.time() - start_editing_rate_ci,
		)
	else:
		removed_keys = remove_optional_artifacts(
			(
				"editing_rate_ci",
				"editing_rate_unconditional_permutation",
				"editing_rate_unconditional_simulations",
			),
			"write_editing_rate_ci is disabled",
		)
		mark_skipped_unless_removed(
			(
				"editing_rate_ci",
				"editing_rate_unconditional_permutation",
				"editing_rate_unconditional_simulations",
				"editing_rate_confidence_intervals_plot",
				"editing_rate_coverage_adjusted_effects_plot",
				"editing_rate_unconditional_permutation_plot",
				"editing_rate_observed_centered_permutation_swarm_plot",
			),
			"write_editing_rate_ci is disabled",
			removed_keys,
		)

	if editing_rate_depth_stability_config.enabled:
		start_editing_rate_depth_stability = time.time()
		depth_stability_plot_objects = write_editing_rate_depth_stability_output(
			output_root=output_root,
			cell_quality_to_analyze=cell_quality_to_analyze,
			min_reads_per_amplicon_per_cell=min_reads_per_amplicon_per_cell,
			config=editing_rate_depth_stability_config,
			n_processes=n_processes,
			significant_amplicons=editing_rate_significant_amplicons,
		)
		filtered_summary_plot_objects.extend(depth_stability_plot_objects)
		mark_written("editing_rate_depth_stability")
		if manifest is not None:
			if depth_stability_plot_objects:
				manifest.mark_written("editing_rate_depth_stability_plot")
			else:
				manifest.mark_skipped(
					"editing_rate_depth_stability_plot", "plot requirements were not met"
				)
		logging.info(
			"Generated editing-rate depth stability analysis in %.2f seconds",
			time.time() - start_editing_rate_depth_stability,
		)
	else:
		removed_keys = remove_optional_artifacts(
			("editing_rate_depth_stability",),
			"write_editing_rate_depth_stability is disabled",
		)
		mark_skipped_unless_removed(
			("editing_rate_depth_stability", "editing_rate_depth_stability_plot"),
			"write_editing_rate_depth_stability is disabled",
			removed_keys,
		)


	# # filtered_read_count_plot_obj = generate_read_depth_boxplots(output_root, cell_quality_to_analyze)
	# filtered_summary_plot_objects.append(filtered_read_count_plot_obj)
   #
	# # Cell Coverage Bar Chart
	# filtered_cell_cov_plot_obj = generate_cell_coverage_plot(output_root, cell_quality_to_analyze)
	# filtered_summary_plot_objects.append(filtered_cell_cov_plot_obj)
#
	# # Avg Read Count Per Amplicon Boxplot
	# filtered_amp_cov_plot_obj = generate_amplicon_coverage_plot(output_root, cell_quality_to_analyze)
	# filtered_summary_plot_objects.append(filtered_amp_cov_plot_obj)
	#
	# # Upset plot displaying edit site combinations
	# filtered_upset_plot_obj = generate_upset_plot(output_root, cell_quality_to_analyze)
	# filtered_summary_plot_objects.append(filtered_upset_plot_obj)
   #
	# # Histogram displaying the frequency of edited site number
	# filtered_hist_plot_obj = generate_edit_histogram(output_root, cell_quality_to_analyze)
	# filtered_summary_plot_objects.append(filtered_hist_plot_obj)
   #
	# # Generate a log read count (y) vs log barcode rank (x) colored by cell quality category
	# filtered_log_log_plot_obj = log_log_plot(parsed_information, output_root, cell_quality_to_analyze, filtered = False)
	# filtered_summary_plot_objects.append(filtered_log_log_plot_obj)
	  #
	# #
	# filtered_cell_per_amp_obj = cell_per_amp_filtered(parsed_information, output_root, cell_quality_to_analyze)
	# filtered_summary_plot_objects.append(filtered_cell_per_amp_obj)
		#
	# #
	# filtered_amp_per_cell_obj = amp_per_cell_filtered(parsed_information, output_root, cell_quality_to_analyze)
	# filtered_summary_plot_objects.append(filtered_amp_per_cell_obj)
#
	# #
	# filtered_mod_pct_plot_obj = mod_per_amp_filtered(parsed_information, output_root, cell_quality_to_analyze)
	# filtered_summary_plot_objects.append(filtered_mod_pct_plot_obj)
#
	# #
	# amp_score_plot_obj = plot_amp_score(output_root)
	# filtered_summary_plot_objects.append(amp_score_plot_obj)


	if suppress_sub_crispresso_plots:
		crispresso_run_names = []
		crispresso_sub_html_files = {}
		filtered_crispresso_sub_html_files = {}
	else:
		crispresso_run_names = amplicon_names #this is a list of target names to display on the report
		crispresso_sub_html_files = {}
		filtered_crispresso_sub_html_files = {}
		for amplicon_name in amplicon_names:
			relative_crispresso_dir = os.path.relpath(crispresso_dir,os.path.dirname(output_root))
			relative_filtered_crispresso_dir = relative_crispresso_dir + ".filtered"
			crispresso_sub_html_files[amplicon_name] = relative_crispresso_dir + "/" "CRISPResso_on_" + amplicon_name + ".html"
			filtered_crispresso_sub_html_files[amplicon_name] = relative_filtered_crispresso_dir + "/" "CRISPResso_on_" + amplicon_name + ".html"


	make_report(report_file=outputs.path("report"),
				report_name = "Dataset Summary Report",
				results_folder='',
				crispresso_run_names=crispresso_run_names,
				crispresso_sub_html_files=filtered_crispresso_sub_html_files,
				summary_plot_objects=filtered_summary_plot_objects)
	mark_written("report")

	if write_h5ad:
		if manifest is not None:
			manifest.set_stage("write_h5ad")
		start_h5ad = time.time()
		h5ad_file = write_h5ad_output(
			output_root=output_root,
			settings_file=settings_file,
			h5ad_output=h5ad_output,
			h5ad_export_config=h5ad_export_config,
			n_processes=n_processes,
		)
		end_h5ad = time.time() - start_h5ad
		mark_written("h5ad")
		logging.info("Generated h5ad output at %s in %.2f seconds", h5ad_file, end_h5ad)
	elif manifest is not None:
		manifest.mark_skipped("h5ad", "write_h5ad is disabled")

	if manifest is not None:
		manifest.set_stage("finalize")
	logging.info('Finished')
