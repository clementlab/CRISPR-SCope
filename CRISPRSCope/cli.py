import logging
from collections import defaultdict
from datetime import datetime
import multiprocessing as mp
import gzip
import os
import re
import sys
import signal
import threading
import subprocess as sb
import traceback
import re
import errno
from typing import Optional
from scipy.stats import multinomial
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import seaborn as sns
import time
from functools import partial
from CRISPResso2 import CRISPRessoShared
import zipfile
from upsetplot import UpSet
from adjustText import adjust_text
import random
import dnaio
from CRISPRSCope import __version__

# Constants and default settings
MIN_TOTAL_READS_PER_BARCODE_DEFAULT = 10
MIN_READS_PER_AMPLICON_PER_CELL_DEFAULT = 0
H5AD_ZYGOSITY_DEFAULTS = {
	"wt_max_mod_pct": 20.0,
	"het_max_mod_pct": 80.0,
	"hom_min_mod_pct": 80.0,
	"compound_het_min_allele2_pct": 20.0,
}
CELL_QUALITY_CODES = {
	"HQ_HI": {"label": "HighScore_HighDepth", "display": "High score / High depth", "color": "#1f77b4"},
	"HQ_LO": {"label": "HighScore_LowDepth",  "display": "High score / Low depth",  "color": "#2ca02c"},
	"LQ_HI": {"label": "LowScore_HighDepth",  "display": "Low score / High depth",  "color": "#ff7f0e"},
	"LQ_LO": {"label": "LowScore_LowDepth",   "display": "Low score / Low depth",   "color": "#d62728"},
}


# Pipeline stage numbering
STAGE_PARSE = 1
STAGE_ALIGN = 2
STAGE_SPLIT = 3
STAGE_FILTER = 4


_SAFE_FN_RE = re.compile(r"[^A-Za-z0-9._-]+")

def _sanitize_token(token: Optional[str]) -> str:
	"""
	Convert an arbitrary string into a filesystem-safe token.

	Replaces whitespace and non-alphanumeric characters with underscores,
	collapses repeated underscores, and removes leading/trailing underscores.

	Parameters
	----------
	token : str or None
		Input string to sanitize (e.g., amplicon name, sample ID,
		output prefix component).

	Returns
	-------
	str
		Sanitized string suitable for use in filenames.
		Returns an empty string if `token` is falsy.

	Notes
	-----
	- Allowed characters after sanitization:
   - Does not modify filesystem state.
	"""
	if not token:
		return ""
	t = str(token).strip()
	# Replace whitespace and unsafe characters with underscore : collaspe multiple underscores
	t = _SAFE_FN_RE.sub("_", t)
	t = re.sub(r"_+", "_", t)
	# Avoid leading/trailing underscores
	return t.strip("_") 

def build_stage_filename(stage: int,
						 tag: str,
						 amplicon: Optional[str] = None,
						 read: Optional[str] = None,
						 ext: str = "fastq.gz",
						 sep: str = ".",
						 output_root: str = None) -> str:
	"""
	Construct a standardized pipeline filename for a processing stage.

	Supports two modes based on `output_root`:

	Directory mode:
		If `output_root` ends with os.sep or is an existing directory:
			<output_root>/<stage>_<tag>.<amplicon>.<read>.<ext>

	Prefix mode:
		Otherwise:
			<parent_dir>/<basename(output_root)>_<stage>_<tag>.<amplicon>.<read>.<ext>

	Parameters
	----------
	stage : int
		Pipeline stage number (non-negative).
		Zero-padded to two digits.

	tag : str
		Short descriptor for the file's purpose
		(e.g., 'parsed_worker0', 'reads_qc_cells').

	amplicon : str, optional
		Amplicon identifier.

	read : str, optional
		Read label (e.g., 'r1', 'r2').

	ext : str, optional
		File extension (without leading dot).

	sep : str, optional
		Separator used between filename components.

	output_root : str
		Base output root or directory path.

	Returns
	-------
	str
		Fully constructed absolute output path.

	Raises
	------
	ValueError
		If stage or tag are invalid.
	ValueError
		If output_root is None.

	Notes
	-----
	- All components are sanitized using `_sanitize_token`.
	- Does not create directories.
	"""
	if output_root is None:
		raise ValueError("output_root is required for build_stage_filename()")

	if not isinstance(stage, int) or stage < 0:
		raise ValueError("stage must be a non-negative integer")
	if not tag or not isinstance(tag, str):
		raise ValueError("tag must be a non-empty string")

	# normalize extension and root
	if ext and ext.startswith("."):
		ext = ext.lstrip(".")
	raw_root = str(output_root)
	abs_root = os.path.abspath(raw_root)

	# determine mode: directory mode if user ended with sep OR path exists as a directory
	explicit_dir = raw_root.endswith(os.sep)
	exists_dir = os.path.isdir(abs_root)
	is_dir_mode = explicit_dir or exists_dir

	# build descriptive base (e.g. "01_parsed_worker0.AMP1.r1")
	parts = [f"{stage:02d}_{_sanitize_token(tag)}"]
	amp = _sanitize_token(amplicon) if amplicon else ""
	if amp:
		parts.append(amp)
	rd = _sanitize_token(read) if read else ""
	if rd:
		parts.append(rd)
	base = sep.join(parts)

	if is_dir_mode:
		out_dir = abs_root.rstrip(os.sep)
		filename = f"{base}.{ext}"
		return os.path.join(out_dir, filename)
	else:
		parent = os.path.dirname(abs_root) or os.getcwd()
		root_base = _sanitize_token(os.path.basename(abs_root))
		filename = f"{root_base}_{base}.{ext}"
		return os.path.join(parent, filename)



def safe_write_path(path: str) -> None:
	"""
	Ensure a file path is writable.

	Creates parent directories if necessary and verifies write
	permission by creating and deleting a temporary test file.

	Parameters
	----------
	path : str
		Target file path.

	Returns
	-------
	None

	Raises
	------
	IOError
		If the path cannot be written to.

	Notes
	-----
	- Creates parent directories if they do not exist.
	- Does not leave temporary files behind.
	"""
	parent = os.path.dirname(path)
	if parent and not os.path.isdir(parent):
		os.makedirs(parent, exist_ok=True)

	# quick write test (create temp file next to target then remove it)
	test_path = path + ".crispresso_write_test"
	try:
		with open(test_path, "w") as fh:
			fh.write("test")
		os.remove(test_path)
	except Exception as e:
		raise IOError(f"Cannot write to path {path!r}: {e}")

def safe_remove(path: str, silent: bool = False) -> bool:
	"""
	Remove a file if it exists.

	Parameters
	----------
	path : str
		File path to remove.

	silent : bool, optional
		If True, suppress exceptions and return False on failure.

	Returns
	-------
	bool
		True if file was removed.
		False if file did not exist or removal failed (silent=True).

	Notes
	-----
	- Logs deletion events.
	- Does not remove directories.
	"""
	if not path:
		return False
	try:
		if os.path.exists(path):
			os.remove(path)
			logging.info("Deleted intermediate file: %s", path)
			return True
		else:
			logging.debug("safe_remove(): file not found, skipping: %s", path)
			return False
	except Exception as e:
		logging.exception("safe_remove(): failed to delete %s: %s", path, e)
		if not silent:
			raise
		return False


def validate_output_root(output_root: str) -> str:
	"""
	Validate and normalize an output root path.

	Determines whether `output_root` represents:
		- A directory mode (ends with os.sep or existing directory), or
		- A prefix mode (file prefix whose parent directory must exist).

	Also verifies write permission by creating and removing a temporary file.

	Parameters
	----------
	output_root : str
		User-provided output root. May be:
			- A directory path (ending with os.sep or existing directory)
			- A file prefix whose parent directory must exist

	Returns
	-------
	str
		Absolute path to the validated output root.

	Raises
	------
	ValueError
		If `output_root` is empty or None.
	FileNotFoundError
		If required directory does not exist.
	PermissionError
		If write permission check fails.

	Notes
	-----
	- Performs a write test by creating a temporary file.
	- Does not create missing directories.
	"""
	if not output_root:
		raise ValueError("output_root must be provided")

	raw = str(output_root)
	abs_root = os.path.abspath(raw)

	explicit_dir = raw.endswith(os.sep)
	exists_dir = os.path.isdir(abs_root)
	is_dir_mode = explicit_dir or exists_dir

	if is_dir_mode:
		target_dir = abs_root.rstrip(os.sep)
		if not os.path.isdir(target_dir):
			raise FileNotFoundError(f"Requested output directory '{target_dir}' does not exist")
		test_path = os.path.join(target_dir, f".cr_write_test_{os.getpid()}")
	else:
		parent = os.path.dirname(abs_root) or os.getcwd()
		if not os.path.isdir(parent):
			raise FileNotFoundError(f"Parent directory '{parent}' for output_root does not exist")
		test_path = os.path.join(parent, f".cr_write_test_{os.getpid()}")

	try:
		with open(test_path, "w") as fh:
			fh.write("ok")
		os.remove(test_path)
	except Exception as e:
		raise PermissionError(f"No write permission for '{output_root}': {e}")

	return abs_root




def main():
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
		write_h5ad, h5ad_output, h5ad_export_config, settings_file
		) = parse_settings(sys.argv)
	end_settings = time.time() - start_settings
	#print(f"Parse Settings: {end_settings}")

	output_root = validate_output_root(output_root)
	
	start_parse_and_align = time.time() 
	aligned_bam, reads_per_cell = parse_and_align_reads(r1,r2,constant1,constant2,output_root,barcode_file,allow_barcode_mismatches,adapter_DNA,bowtie2_index,n_processes,keep_intermediate_files)
	end_parse_and_align = time.time() - start_parse_and_align
	logging.info(f"Parse and Align Reads: {end_parse_and_align}")


	start_split_reads = time.time()
	amplicon_names, amplicon_information, amplicon_info_file = split_reads_by_amplicon(aligned_bam, output_root, amplicon_file, alt_alleles_file, primer_lookup_len, amp_file_dir, bowtie2_index, adapter_DNA, n_processes, keep_intermediate_files, reads_per_cell, min_total_reads_per_barcode, assign_reads_to_all_possible_amplicons)
	end_split_reads = time.time() - start_split_reads
	logging.info(f"Split Reads by Amplicon: {end_split_reads}")
	
#    print(f"Line266\n{amplicon_names=}\n{amplicon_information=}\n{amplicon_info_file=}")

	start_crispresso = time.time()
	crispresso_information = run_crispresso_commands(amplicon_names,amplicon_information,output_root,crispresso_dir,suppress_sub_crispresso_plots,n_processes, alleles = False)
	end_crispresso = time.time() - start_crispresso
	logging.info(f"Run CRISPResso: {end_crispresso}")

	start_parse_crispresso = time.time()
	parsed_information = parse_crispresso_outputs(amplicon_names,amplicon_information,amplicon_info_file,
												  crispresso_information,output_root, min_total_reads_per_barcode, min_reads_per_amplicon_per_cell,
												  n_processes=n_processes,
												  ignore_substitutions=ignore_substitutions)
	end_parse_crispresso = time.time() - start_parse_crispresso
	logging.info(f"Parse CRISPResso Outputs: {end_parse_crispresso}")

	start_filter_amplicon = time.time()
	filter_amplicon_reads(output_root, parsed_information, amplicon_names, cell_quality_to_analyze, n_processes)
	end_filter_amplicon = time.time() - start_filter_amplicon
	logging.info(f"Filter Amplicon Reads: {end_filter_amplicon}")
	
	start_run_crispresso2 = time.time()
	crispresso_filtered_information = run_crispresso_commands(amplicon_names,amplicon_information,output_root,crispresso_dir,suppress_sub_crispresso_plots,n_processes, alleles = True)
	end_run_crispresso2 = time.time() - start_run_crispresso2
	logging.info(f"Run CRISPResso 2: {end_run_crispresso2}")

	filtered_summary_plot_objects = []

	filtered_read_count_plot_obj = generate_read_depth_boxplots(output_root, cell_quality_to_analyze)
	if filtered_read_count_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_read_count_plot_obj)

	# Cell Coverage Bar Chart
	filtered_cell_cov_plot_obj = generate_cell_coverage_plot(output_root, cell_quality_to_analyze)
	if filtered_cell_cov_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_cell_cov_plot_obj)

	# Avg Read Count Per Amplicon Boxplot
	filtered_amp_cov_plot_obj = generate_amplicon_coverage_plot(output_root, cell_quality_to_analyze)
	if filtered_amp_cov_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_amp_cov_plot_obj)

	# Upset plot displaying edit site combinations
	filtered_upset_plot_obj = generate_upset_plot(output_root, cell_quality_to_analyze)
	if filtered_upset_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_upset_plot_obj)

	# Histogram displaying the frequency of edited site number
	filtered_hist_plot_obj = generate_edit_histogram(output_root, cell_quality_to_analyze)
	if filtered_hist_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_hist_plot_obj)

	# Generate a log read count (y) vs log barcode rank (x) colored by cell quality category
	filtered_log_log_plot_obj = log_log_plot(parsed_information, output_root, cell_quality_to_analyze, filtered = False)
	if filtered_log_log_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_log_log_plot_obj)

	#
	filtered_cell_per_amp_obj = cell_per_amp_filtered(parsed_information, output_root, cell_quality_to_analyze)
	if filtered_cell_per_amp_obj is not None:
		filtered_summary_plot_objects.append(filtered_cell_per_amp_obj)

	#
	filtered_amp_per_cell_obj = amp_per_cell_filtered(parsed_information, output_root, cell_quality_to_analyze)
	if filtered_amp_per_cell_obj is not None:
		filtered_summary_plot_objects.append(filtered_amp_per_cell_obj)

	#
	filtered_mod_pct_plot_obj = mod_per_amp_filtered(parsed_information, output_root, cell_quality_to_analyze)
	if filtered_mod_pct_plot_obj is not None:
		filtered_summary_plot_objects.append(filtered_mod_pct_plot_obj)

	#
	amp_score_plot_obj = plot_amp_score(output_root)
	if amp_score_plot_obj is not None:
		filtered_summary_plot_objects.append(amp_score_plot_obj)


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
  
	
	make_report(report_file=output_root+".html",
				report_name = "Dataset Summary Report",
				results_folder='',
				crispresso_run_names=crispresso_run_names,
				crispresso_sub_html_files=filtered_crispresso_sub_html_files,
				summary_plot_objects=filtered_summary_plot_objects)

	if write_h5ad:
		start_h5ad = time.time()
		h5ad_file = write_h5ad_output(
			output_root=output_root,
			settings_file=settings_file,
			h5ad_output=h5ad_output,
			h5ad_export_config=h5ad_export_config,
			n_processes=n_processes,
		)
		end_h5ad = time.time() - start_h5ad
		logging.info("Generated h5ad output at %s in %.2f seconds", h5ad_file, end_h5ad)

	logging.info('Finished')
		
	
def generate_amplicon_coverage_plot(output_root, cell_quality_to_analyze):
	"""
		Generate a box-and-swarm plot of average read coverage per amplicon.

	This plot displays the distribution of mean read counts across
	amplicons, restricted to barcodes whose quality category is in
	`cell_quality_to_analyze`.

	Coverage is computed as:
		(sum of totCount.<amplicon> across selected cells)
		divided by
		number of selected cells.

	The plot highlights the three highest and three lowest coverage
	amplicons.

	Parameters
	----------
	output_root : str
		Base output prefix used to locate input summary files and
		write plot outputs.
		Required input files:
			- <output_root>.editingSummary.txt
			- <output_root>.amplicon_score.txt

	cell_quality_to_analyze : list[str]
		List of quality category short codes to include.
		Expected values:
			{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Reads:
		- <output_root>.editingSummary.txt
		- <output_root>.amplicon_score.txt

	Writes:
		- <output_root>.09_AmpliconCoverage.pdf
		- <output_root>.09_AmpliconCoverage.png

	Notes
	-----
	- Only barcodes present in the amplicon score file and matching
	  the requested quality categories are included.
	- Does not modify input data.
	- Uses matplotlib and adjustText for annotation.
	"""

	plt.clf()
	plt.cla()
	# Read in editingSummary file
	editing = pd.read_csv(output_root + ".editingSummary.txt", sep = "\t", index_col = 0)
	# Read in amplicon score file
	amplicon = pd.read_csv(output_root + ".amplicon_score.txt", sep = "\t", index_col = 0)
	
	#valid_barcodes = amplicon[amplicon['Color'] == "High Score / High Reads"].index.tolist()
	valid_barcodes = amplicon[amplicon['Color'].isin(cell_quality_to_analyze)].index.tolist()


	editing = editing[editing.index.isin(valid_barcodes)]
	editingSummary_count = editing.filter(like = "totCount")
	
	amplicon_column_sum = editingSummary_count.sum(axis = 0) / len(editingSummary_count)
	amplicon_column_sum_sorted = amplicon_column_sum.sort_values(ascending = False)
	amplicon_column_sum_sorted.index = amplicon_column_sum_sorted.index.str.replace('totCount.', '')
	
	# Extract top and bottom 3 amplicons into a dictionary
	indices = list(range(0, 3)) + list(range(len(amplicon_column_sum_sorted) - 3, len(amplicon_column_sum_sorted)))
	points_to_label = amplicon_column_sum_sorted[indices].to_dict()
	
	# Create the figure and axis
	fig, ax = plt.subplots(figsize = (12,12))

	# Overlay the boxplot
	ax.boxplot(amplicon_column_sum_sorted.values, vert=True, patch_artist=True, widths=0.1, boxprops=dict(facecolor='white'))

	# Create the swarm plot
	y = amplicon_column_sum_sorted.values
	x = [1 + random.uniform(-0.04, 0.04) for _ in range(len(y))]  # Add jitter for the swarm plot
	ax.scatter(x, y, alpha=0.6, zorder = 2)

	texts = []
	for i, target in enumerate(amplicon_column_sum_sorted.index):
		if target in points_to_label.keys():
			texts.append(ax.text(x[i], points_to_label[target], target, fontsize = 12))

	# Adjust text to avoid overlap
	adjust_text(texts, arrowprops=dict(arrowstyle="->", color = 'red', lw=1))

	# Add title and labels
	ax.set_title("Average Read Count Per Amplicon", fontsize = 24)
	ax.set_ylabel("Read Count", fontsize = 22)
	ax.set_xlabel("Amplicons", fontsize = 22)
	ax.set_xticks([])
	
	amplicon_cov_plot_root = output_root + ".09_AmpliconCoverage"
	plt.savefig(amplicon_cov_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(amplicon_cov_plot_root+".png", pad_inches = 1, bbox_inches = "tight") 
	
	logging.info("Finished generating the average amplicon coverage plot.")
	summary_plot_obj = PlotObject(
			plot_name = amplicon_cov_plot_root,
			plot_title = 'Average Amplicon Read Coverage',
			plot_label = 'The average read counts covering each amplicon.',
			plot_datas = [
				("Amplicon Read Counts (totCounts)", output_root + ".editingSummaryPseudobulk.txt"),
				("Filtered cell barcodes", output_root + ".amplicon_score.txt"),
				]
			)
	return summary_plot_obj 

def generate_read_depth_boxplots(output_root, cell_quality_to_analyze):
	"""
		Generate a boxplot comparing per-barcode read depth between
	selected quality categories and all remaining barcodes.

	This function separates barcodes into two groups:
		- "High" group: barcodes whose quality category is in
		  `cell_quality_to_analyze`.
		- "Low" group: all other barcodes.

	For each group, total read counts per amplicon (columns matching
	'totCount.*') are reshaped into long format and plotted as a
	box-and-whisker comparison.

	Parameters
	----------
	output_root : str
		Base output prefix used to locate input summary files and
		write plot outputs.
		Required input files:
			- <output_root>.editingSummary.txt
			- <output_root>.amplicon_score.txt

	cell_quality_to_analyze : list[str]
		List of quality category short codes defining the "High" group.
		Expected values:
			{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Reads:
		- <output_root>.editingSummary.txt
		- <output_root>.amplicon_score.txt

	Writes:
		- <output_root>.09_CellCoverageBoxplot.pdf
		- <output_root>.09_CellCoverageBoxplot.png

	Notes
	-----
	- Uses seaborn.boxplot for visualization.
	- All amplicons are pooled when computing per-barcode read depth.
	- Barcodes not present in the amplicon score file are excluded.
	- Does not modify input data.
	"""
	plt.clf()
	plt.cla()
	
	# Read in editingSummary file
	editing = pd.read_csv(output_root + ".editingSummary.txt", sep = "\t", index_col = 0)
	# Read in amplicon score file
	amplicon = pd.read_csv(output_root + ".amplicon_score.txt", sep = "\t", index_col = 0)
	
	valid_colors = set(CELL_QUALITY_CODES)
	high_colors = set(cell_quality_to_analyze)
	low_colors = valid_colors - high_colors

	classified_barcodes = amplicon[amplicon['Color'].isin(valid_colors)]
	editing = editing[editing.index.isin(classified_barcodes.index)]

	high_barcodes = classified_barcodes[classified_barcodes['Color'].isin(high_colors)].index
	low_barcodes = classified_barcodes[classified_barcodes['Color'].isin(low_colors)].index

	high_totals = editing[editing.index.isin(high_barcodes)].filter(like="totCount").sum(axis=1)
	low_totals = editing[editing.index.isin(low_barcodes)].filter(like="totCount").sum(axis=1)

	plot_df = pd.concat([
		pd.DataFrame({"Read Count": high_totals, "Barcode Quality": "High"}),
		pd.DataFrame({"Read Count": low_totals, "Barcode Quality": "Low"}),
	], ignore_index=True)

	# Generate boxplot
	plt.figure(figsize=(12, 12))
	sns.boxplot(x="Barcode Quality", y="Read Count", data=plot_df)
	plt.ylabel("Read Count", fontsize=22)
	plt.xlabel("Barcode Quality", fontsize=22)
	plt.title("Read Counts by Barcode Quality", fontsize=24)
	plt.yticks(fontsize=16)
	plt.xticks(fontsize=16)
	plt.tight_layout()

	cell_cov_plot_root = output_root + ".09_CellCoverageBoxplot"
	plt.savefig(cell_cov_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(cell_cov_plot_root+".png", pad_inches = 1, bbox_inches = "tight") 
	
	logging.info("Finished generating the read counts per barcode barplot plot.")
	summary_plot_obj = PlotObject(
			plot_name = cell_cov_plot_root,
			plot_title = 'Read Counts per Barcode Boxplot',
			plot_label = 'A boxplot displaying the read counts per barcode in high quality and low quality cells.',
			plot_datas = [
				("Barcode Read Counts (totCounts)", output_root + ".editingSummaryPseudobulk.txt"),
				("Filtered cell barcodes", output_root + ".amplicon_score.txt"),
				]
			)
	return summary_plot_obj   
	
def generate_cell_coverage_plot(output_root, cell_quality_to_analyze):
	"""
	Generate a bar plot of average total read count per barcode
	for selected and non-selected quality categories.

	Barcodes are separated into two groups:
		- "High Quality" group: barcodes whose 'Color' value is in
		  `cell_quality_to_analyze`.
		- "Low Quality" group: all remaining barcodes.

	For each group, the total read count per barcode is computed as
	the row-wise sum of all columns matching 'totCount.*'. The mean
	total read count across barcodes in each group is then plotted.

	Parameters
	----------
	output_root : str
		Base output prefix used to locate input summary files and
		write plot outputs.
		Required input files:
			- <output_root>.editingSummary.txt
			- <output_root>.amplicon_score.txt

	cell_quality_to_analyze : list[str]
		List of quality category short codes defining the "High Quality" group.
		Expected values:
			{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Reads:
		- <output_root>.editingSummary.txt
		- <output_root>.amplicon_score.txt

	Writes:
		- <output_root>.08_CellCoverage.pdf
		- <output_root>.08_CellCoverage.png

	Notes
	-----
	- Uses seaborn.barplot for visualization.
	- The "Low Quality" group includes all barcodes not explicitly
	  listed in `cell_quality_to_analyze`.
	- Barcodes not present in the amplicon score file are excluded.
	- Does not modify input data.
	"""
	plt.clf()
	plt.cla()
	
	# Read in editingSummary file
	editing = pd.read_csv(output_root + ".editingSummary.txt", sep = "\t", index_col = 0)
	# Read in amplicon score file
	amplicon = pd.read_csv(output_root + ".amplicon_score.txt", sep = "\t", index_col = 0)
	
	valid_colors = set(CELL_QUALITY_CODES)
	high_colors = set(cell_quality_to_analyze)
	low_colors = valid_colors - high_colors

	classified_barcodes = amplicon[amplicon['Color'].isin(valid_colors)]
	editing = editing[editing.index.isin(classified_barcodes.index)]

	high_barcodes = classified_barcodes[classified_barcodes['Color'].isin(high_colors)].index
	low_barcodes = classified_barcodes[classified_barcodes['Color'].isin(low_colors)].index

	high_qual_edit = editing[editing.index.isin(high_barcodes)].filter(like="totCount")
	low_qual_edit = editing[editing.index.isin(low_barcodes)].filter(like="totCount")
	
	# Calculate read count average in high quality barcodes
	high_qual_sum = high_qual_edit.sum(axis = 1)
	high_qual_avg = high_qual_sum.mean()
	
	# Calculate read count average in non-high quality barcodes
	low_qual_sum = low_qual_edit.sum(axis = 1)
	low_qual_avg = low_qual_sum.mean()
	
	# Create a DF for plotting
	results = pd.DataFrame({"Average Read Count": [high_qual_avg, low_qual_avg], "Data Source": ["High Quality Barcodes", "Low Quality Barcodes"]})
	
	# Generate barplot
	plt.figure(figsize = (12,12))
	sns.barplot(x = "Data Source", y = "Average Read Count", data = results)
	plt.ylabel("Read Count", fontsize = 22)
	plt.xlabel("Data Source", fontsize = 22)
	plt.title("Average Read Count Per Barcode", fontsize = 24)
	plt.yticks(fontsize = 16)
	plt.xticks(fontsize = 16)
	plt.tight_layout()
	
	
	cell_cov_plot_root = output_root + ".08_CellCoverage"
	plt.savefig(cell_cov_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(cell_cov_plot_root+".png", pad_inches = 1, bbox_inches = "tight") 
	
	logging.info("Finished generating the read counts per barcode barplot plot.")
	summary_plot_obj = PlotObject(
			plot_name = cell_cov_plot_root,
			plot_title = 'Read Counts per Barcode',
			plot_label = 'A barplot displaying the difference in average read count per barcode in high quality and low quality cells.',
			plot_datas = [
				("Barcode Read Counts (totCounts)", output_root + ".editingSummaryPseudobulk.txt"),
				("Filtered cell barcodes", output_root + ".amplicon_score.txt"),
				]
			)
	return summary_plot_obj 
	
def generate_edit_histogram(output_root, cell_quality_to_analyze):
	"""
	Generate a histogram showing the number of edited amplicons per barcode.

	For each selected barcode (based on `cell_quality_to_analyze`),
	an amplicon is considered "edited" if its corresponding
	'modPct.<amplicon>' value is greater than zero.

	The number of edited amplicons is counted per barcode, and a
	histogram is plotted showing the distribution of edited-site
	counts across barcodes.

	Parameters
	----------
	output_root : str
		Base output prefix used to locate input summary files and
		write plot outputs.
		Required input files:
			- <output_root>.editingSummary.txt
			- <output_root>.amplicon_score.txt

	cell_quality_to_analyze : list[str]
		List of quality category short codes to include.
		Expected values:
			{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Reads:
		- <output_root>.editingSummary.txt
		- <output_root>.amplicon_score.txt

	Writes:
		- <output_root>.07_EditHistogram.pdf
		- <output_root>.07_EditHistogram.png

	Notes
	-----
	- Only barcodes whose 'Color' is in `cell_quality_to_analyze`
	  are included.
	- An amplicon is counted as edited if modPct > 0 (no additional
	  frequency threshold is applied here).
	- Histogram bins are centered on integer counts of edited sites.
	- Does not modify input data.
	"""
	plt.clf()
	plt.cla()
	# Read in editingSummary file
	editing = pd.read_csv(output_root + ".editingSummary.txt", sep = "\t", index_col = 0) 
	# Read in amplicon score file
	amplicon = pd.read_csv(output_root + ".amplicon_score.txt", sep = "\t", index_col = 0)
	
	barcodes = amplicon[amplicon['Color'].isin(cell_quality_to_analyze)].index.tolist()
	editing = editing[editing.index.isin(barcodes)]
		
	editing = editing.filter(like = "modPct")
	
	editing = editing > 0
	
	row_sum = editing.sum(axis = 1)
	
	# Calculate the bin edges
	bin_edges = [x - 0.5 for x in range(0, max(row_sum) + 2)]
   
	plt.figure(figsize = (12,12)) 
	plt.hist(row_sum, bins = bin_edges, edgecolor = "black")
	plt.xticks(range(0, max(row_sum) + 1), fontsize  = 16)
	plt.yticks(fontsize = 16)
	plt.title("Number of Edited Sites in a Barcode", fontsize = 24)
	plt.xlabel("Number of Edited Sites", fontsize = 22)
	plt.ylabel("Number of Barcodes", fontsize = 22)
	plt.tight_layout()
	
	hist_plot_root = output_root + ".07_EditHistogram"
	plt.savefig(hist_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(hist_plot_root+".png", pad_inches = 1, bbox_inches = "tight") 
	
	logging.info("Finished generating the edit count histogram plot")
	summary_plot_obj = PlotObject(
			plot_name = hist_plot_root,
			plot_title = 'Editing Count Histogram',
			plot_label = 'A histogram that displays the number of edited sites in each barcode.',
			plot_datas = [
				("Modification Percentages (modPct)", output_root + ".editingSummary.txt"),
				("Filtered cell barcodes", output_root + ".amplicon_score.txt")
				]
			)
	return summary_plot_obj 
	
def generate_upset_plot(output_root, cell_quality_to_analyze):
	"""
	Generate an UpSet plot of the most frequent edit-site combinations.

	For barcodes whose 'Color' value is in `cell_quality_to_analyze`,
	each amplicon is considered "edited" if its corresponding
	'modPct.<amplicon>' value is greater than zero.

	A boolean edit matrix is constructed (barcode × amplicon),
	and the five most frequent edit combinations (based on exact
	boolean tuples across amplicons) are identified. The dataset
	is restricted to these top five combinations, and an UpSet
	plot is generated to visualize the intersections.

	Parameters
	----------
	output_root : str
		Base output prefix used to locate input summary files and
		write plot outputs.
		Required input files:
			- <output_root>.editingSummary.txt
			- <output_root>.amplicon_score.txt

	cell_quality_to_analyze : list[str]
		List of quality category short codes to include.
		Expected values:
			{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Reads:
		- <output_root>.editingSummary.txt
		- <output_root>.amplicon_score.txt

	Writes:
		- <output_root>.06_EditCombinations.pdf
		- <output_root>.06_EditCombinations.png

	Notes
	-----
	- Only barcodes matching `cell_quality_to_analyze` are included.
	- An amplicon is considered edited if modPct > 0.
	- Only the five most frequent exact edit combinations are shown.
	- Amplicons with no edits across the selected combinations are removed.
	- Uses the `upsetplot.UpSet` visualization library.
	- Does not modify input data.
	"""
	plt.clf()
	plt.cla() 
	# Read in editingSummary file
	editing = pd.read_csv(output_root + ".editingSummary.txt", sep = "\t", index_col = 0) 
	# Read in amplicon score file
	amplicon = pd.read_csv(output_root + ".amplicon_score.txt", sep = "\t", index_col = 0)
	
	# Filter barcodes if required 
	#if filtered:
	#    barcodes = amplicon[amplicon['Color'].isin(["High Score / High Reads", "High Score / Low Reads"])].index.tolist()
	#    editing = editing[editing.index.isin(barcodes)]
	barcodes = amplicon[amplicon['Color'].isin(cell_quality_to_analyze)].index.tolist()
	editing = editing[editing.index.isin(barcodes)]
	
	#high_qual = amplicon[amplicon['Color'].isin(cell_quality_to_analyze)]
 
	# Reduce to modified columns
	editing = editing.filter(like = "modPct") 
	editing.columns = editing.columns.str.replace('modPct.', '')
	editing.fillna(0,inplace=True)
	
	# Convert to a boolean matrix
	editing = editing.applymap(lambda x: False if x == 0 else True)
	
	# Find top 5 combinations of edit values
	counts = editing.apply(tuple, axis = 1).value_counts()
	top5 = counts.nlargest(5)
	
	# Filter editing matrix to top 5 combinations
	editing = editing[editing.apply(tuple, axis = 1).isin(top5.index)]
	
	# Remove columns where all modification values == False
	editing = editing.loc[:, ~(editing == False).all()] 
	
	first = True
	for col in editing.columns:
		if first:
			editing = editing.set_index(editing[col] > 0)
			first = False
		else:
			editing = editing.set_index(editing[col] > 0, append = True)
	
	upset = UpSet(editing, orientation = "horizontal", sort_by = "cardinality", show_counts = True)
	
	fig, ax = plt.subplots(figsize = (12, 12))
	upset.plot(fig = fig)
	
	for spine in ax.spines.values():
		spine.set_visible(False)
		ax.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False, labelsize = 16)
	plt.title("Editing Sites and Intersections", fontsize=24)

	
	upset_plot_root = output_root + ".06_EditCombinations"
	plt.savefig(upset_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(upset_plot_root+".png", pad_inches = 1, bbox_inches = "tight") 
	
	logging.info("Finished generating the edit combination upset plot")
	summary_plot_obj = PlotObject(
			plot_name = upset_plot_root,
			plot_title = 'Editing Sites and Intersections',
			plot_label = 'An upset plot that displays the most common edits and edit combinations.',
			plot_datas = [
				("Modification Percentages (modPct)", output_root + ".editingSummary.txt"),
				("Filtered cell barcodes", output_root + ".amplicon_score.txt"),
				]
			)
	return summary_plot_obj


def _filter_to_hq_reads_at_single_amplicon(args):
	"""
	Filter paired per-amplicon FASTQ files (R1/R2) to only reads whose barcode
	is in `barcodes`.

	Parameters
	----------
	args : tuple
		(amp, amplicon_dir, barcodes)

	Returns
	-------
	dict
		On success:
			{
				'Amplicon': amp,
				'Status': 'Success',
				'R1': out1_name,
				'R2': out2_name
			}

		On failure:
			{
				'Amplicon': amp,
				'Status': 'Failed',
				'Reason': 'MissingInput' | 'BarcodeMismatch'
			}
	"""    

	amp, amplicon_dir, barcodes = args
	#out1_name = os.path.join(amplicon_dir, f"filtered.{amp}.r1.fq.gz")
	#out2_name = os.path.join(amplicon_dir, f"filtered.{amp}.r2.fq.gz")
	out1_name = build_stage_filename(
		stage = STAGE_FILTER,
		tag = "reads_qc_cells",
		amplicon = amp,
		read = "r1",
		ext = "fq.gz",
		output_root = amplicon_dir
	)
	out2_name = build_stage_filename(
		stage = STAGE_FILTER,
		tag = "reads_qc_cells",
		amplicon = amp,
		read = "r2",
		ext = "fq.gz",
		output_root = amplicon_dir
	)
	

	# Skip if already exists
	if os.path.isfile(out1_name) and os.path.isfile(out2_name):
		return {
			"Amplicon" : amp,
			"Status" : "Success",
			"R1" : out1_name,
			"R2" : out2_name,
			"Reason" : "AlreadyExists"
		}
		
	#in_r1 = os.path.join(amplicon_dir, f"{amp}.r1.fq.gz")
	#in_r2 = os.path.join(amplicon_dir, f"{amp}.r2.fq.gz")
	
	in_r1 = build_stage_filename(
		stage = STAGE_SPLIT,
		tag = "reads_all_cells",
		amplicon = amp,
		read = "r1",
		ext = "fq.gz",
		output_root = amplicon_dir
	)
	in_r2 = build_stage_filename(
		stage = STAGE_SPLIT,
		tag = "reads_all_cells",
		amplicon = amp,
		read = "r2",
		ext = "fq.gz",
		output_root = amplicon_dir
	)

   
	# Checking if reads exist
	if not os.path.isfile(in_r1) or not os.path.isfile(in_r2):
		missing = [p for p in (in_r1, in_r2) if not os.path.isfile(p)]
		return {
			'Amplicon': amp,
			'Status': 'Failed',
			'Reason': 'MissingInput',
			'MissingFiles': missing
		}
		
	# Main processing
	with dnaio.open(in_r1, fileformat="fastq") as reader1, \
		 dnaio.open(in_r2, fileformat="fastq") as reader2, \
		 dnaio.open(out1_name, mode="w", fileformat="fastq") as writer1, \
		 dnaio.open(out2_name, mode="w", fileformat="fastq") as writer2:

		for rec1, rec2 in zip(reader1, reader2):
			barcode1 = rec1.name.split(":")[-1]
			barcode2 = rec2.name.split(":")[-1]
			
			if barcode1 != barcode2:
				logging.error(
					"Amplicon %s: barcode mismatch (%s != %s)",
					amp, barcode1, barcode2
				)
				return {
					"Amplicon" : amp,
					"Status" : "Failed",
					"Reason" : "BarcodeMismatch",
					"Barcode1" : barcode1,
					"Barcode2" : barcode2
				}
		
			if barcode1 in barcodes:
				writer1.write(rec1)
				writer2.write(rec2)
				
	return {
		"Amplicon" : amp,
		"Status" : "Success",
		"R1" : out1_name,
		"R2" : out2_name
	}
				

def _filter_to_hq_alleles_at_single_amplicon(args):
	"""
	Filter a single per-amplicon allele FASTQ to only reads with barcodes in `barcodes`.

	This function reads one input FASTQ (may be gzipped), writes one filtered
	output FASTQ (gzipped), and is safe to run inside a multiprocessing pool.
	It returns a structured result dict rather than raising exceptions so that
	a parent caller can aggregate successes/failures without the pool crashing.

	Parameters
	----------
	args : tuple
		A 3-tuple ``(amp, amplicon_dir, barcodes)``:
		- amp : str
			Amplicon name (used to build filenames).
		- amplicon_dir : str
			Directory where per-amplicon FASTQ files live and where filtered
			output will be written.
		- barcodes : set
			Set of valid barcodes to keep (strings).

	Returns
	-------
	dict
		Dictionary with at least the following keys:
		- 'Amplicon': amp
		- 'Status': 'Success' or 'Failed'
		- 'In': path to the input FASTQ (if found)
		- 'Out': path to the written output FASTQ (on success)
		- 'Reason': short reason if failed (e.g., 'MissingInput', 'Exception')
		- 'Exception': optional long exception string (traceback) on unexpected errors
	"""
	# Writes out alleles for each high quality barcode at an amplicon
	amp_dict, amplicon_dir, barcodes = args
	amp = amp_dict['Amplicon']
   
	amplicon_out_file = build_stage_filename(
		stage = STAGE_FILTER,
		tag = "alleles_qc_cells",
		amplicon = amp,
		ext = "fq.gz",
		output_root = amplicon_dir
	) 
	
	amplicon_in_file = build_stage_filename(
		stage = STAGE_SPLIT,
		tag = "alleles_all_cells",
		amplicon = amp,
		ext = "fq",
		output_root = amplicon_dir
	)
	
	
	# If outfile already exists, return success immediately
	if os.path.isfile(amplicon_out_file):
		return{
			"Amplicon" : amp,
			"Status" : "Success",
			"In" : amplicon_in_file,
			"Out" : amplicon_out_file,
			"Reason" : "AlreadyExists"
		}
	# If input file does not exists, return a warning
	if not os.path.isfile(amplicon_in_file):
		logging.warning("Amplicon %s: allele input file missing: %s", amp, amplicon_in_file)
		return {
			"Amplicon" : amp,
			"Status" : "Failed",
			"In" : amplicon_in_file,
			"Out" : amplicon_out_file,
			"Reason" : "MissingInput"
		}
	
	try:
		# Read input allele fq and return those with the barcodes of interest
		with dnaio.open(amplicon_in_file, fileformat = 'fastq') as reader, dnaio.open(amplicon_out_file, mode = 'w', fileformat = 'fastq') as writer:
			written = 0
			for rec in reader:
				try:
					barcode = rec.name.split(":")[-2]
					#print(f"Got Here : 657 : Barcode : {barcode}")
					if barcode in barcodes:
						writer.write(rec)
						written += 1
				except Exception:
					logging.debug(
						"Amplicon %s: skipping record with malformed read name: %r",
						amp, getattr(rec, "name", None)
					)
					continue
						
			#logging.info("Amplicon %s: wrote %d filtered alleles to %s", amp, written, amplicon_out_file)
		return {
			"Amplicon" : amp,
			'Status' : 'Success',
			'In' : amplicon_in_file,
			'Out' : amplicon_out_file,
			'Written' : written
		}
	except FileNotFoundError as e:
		logging.warning("Amplicon %s: FileNotFoundError process %s: %s", amp, amplicon_in_file, e)
		return {
			"Amplicon" : amp, 
			'Status' : 'Failed',
			'In' : amplicon_in_file,
			'Out' : amplicon_out_file,
			'Reason' : "FileNotFoundError",
			'Exception' : str(e)
		}
	
	except Exception:
		logging.exception("Amplicon %s: unexpected error while filtering alleles.", amp)
	return {
		"Amplicon": amp,
		"Status": "Failed",
		"In": amplicon_in_file,
		"Out": amplicon_out_file,
		"Reason": "Exception",
		"Exception": traceback.format_exc()
	}
   
	
   
def filter_amplicon_reads(output_root, parsed_information, amplicon_names,
						  cell_quality_to_analyze, n_processes):
	"""
	Filter per-amplicon FASTQs and allele files to only include reads from
	barcodes classified as high-quality, running the work in parallel.

	The function performs two phases:
	  1. For each amplicon, filter per-amplicon R1/R2 FASTQs to keep only
		 reads from barcodes in the requested quality set.
	  2. For each amplicon, filter per-amplicon allele FASTQ similarly.

	Parameters
	----------
	output_root : str
		Root path for pipeline outputs (used to construct amplicon directory).
	parsed_information : pandas.DataFrame
		DataFrame indexed by barcode with a 'Color' column used to select
		barcodes to keep (high-quality ones).
	amplicon_names : list[str]
		Names of amplicons available under the amplicon directory.
	cell_quality_to_analyze : list-like
		Values of the 'Color' column that indicate high-quality cells to keep.
	n_processes : int
		Number of worker processes to use for parallel filtering.

	Returns
	-------
	dict
		Summary with the following keys:
		  - 'read_filter_successes' : list of success dicts from read filtering
		  - 'read_filter_failures' : list of failure dicts from read filtering
		  - 'allele_filter_successes': list of success dicts from allele filtering
		  - 'allele_filter_failures': list of failure dicts from allele filtering

	Notes
	-----
	- Worker functions should return structured dictionaries with keys:
	  'Amplicon', 'Status', and other metadata (In, Out, Reason, Exception).
	- This function does not change the parsed_metrics/reads_per_cell mapping; it
	  only writes filtered per-amplicon FASTQs for downstream CRISPResso steps.
	"""
	amplicon_dir = output_root + ".seq_by_amplicon/"
	barcodes = parsed_information.loc[parsed_information['Color'].isin(cell_quality_to_analyze)].index.tolist()
	barcodes = set(barcodes)

	logging.info("Filtering for reads from %d high quality barcodes...", len(barcodes))
	
	worker_args = [(amp, amplicon_dir, barcodes) for amp in amplicon_names]
	
	# Filtering to reads from high quality barcodes at each amplicon
	with mp.Pool(n_processes) as pool:
		amp_results = pool.map(_filter_to_hq_reads_at_single_amplicon, worker_args)
		

	amp_success = [r for r in amp_results if r.get('Status') == 'Success']
	amp_failure = [r for r in amp_results if r.get('Status') != 'Success']
	
	logging.info("Read filtering for high quality barcodes complete: %s successes, %s faltures", len(amp_success), len(amp_failure))
	
	if amp_failure:
		logging.warning("High-quality read filtering failed: %s",[(f.get("Amplicon"), f.get("Reason")) for f in amp_failure[:10]])
	
	# === Phase 2: Filter from all reads of high quality barcodes at an amplicon to the alleles of the barcodes at the amplicon ===
	# Only run allele filtering for amplicons where the read-filter step succeeded.
	successful_amplicons = [r.get('Amplicon') for r in amp_success if r.get('Status') == 'Success']

	if not successful_amplicons: # No amplicons filtered read amplicon files were created
		failures_file = os.path.join(output_root + ".read_filter_failure.txt")
		try:
			with open(failures_file, "w") as fh:
				fh.write("Amplicon\tStatus\tReason\n")
				for f in amp_failure:
					fh.write(f"{f.get('Amplicon')}\t{f.get('Status')}\t{f.get('Reason')}\n")
			logging.error("No amplicons succeeded in read-filtering to high quality barcodes. Wrote failure details to %s", failures_file)
		except Exception:
			logging.exception("Failed to write failure summary to disk (continuing to raise error).")

		raise RuntimeError(
		   f"No amplicons completed read filtering to high quality barcodes successfully. "
		   f"See {failures_file} for further details."
		)
	   
	logging.info("Filtering to the alleles of high quality barcodes across %d amplicons. Skipping %d failed amplicons", len(amp_success), len(amp_failure))
	
	allele_worker_args = [(amp, amplicon_dir, barcodes) for amp in amp_success]
	
	with mp.Pool(n_processes) as pool:
		allele_results = pool.map(_filter_to_hq_alleles_at_single_amplicon, allele_worker_args)
	

def parse_settings(args):
	
	"""
	 Parse pipeline settings from a settings file and optional CLI flags.

	The first CLI argument must be a tab-delimited settings file. The function:
	  1. Parses key/value pairs from the file.
	  2. Validates required inputs and file existence.
	  3. Applies defaults where appropriate.
	  4. Creates required output directories if missing.
	  5. Verifies required external tools (bowtie2, CRISPResso2).

	Parameters
	----------
	args : list[str]
		Command-line arguments (typically sys.argv).

	Returns
	-------
	tuple
		(
			r1 : str,
			r2 : str,
			constant1 : str,
			constant2 : str,
			allow_barcode_mismatches : bool,
			barcode_file : str,
			amplicon_file : str,
			primer_lookup_len : int,
			adapter_DNA : str,
			amp_file_dir : str,
			alt_alleles_file : str,
			bowtie2_index : str,
			crispresso_dir : str,
			output_root : str,
			n_processes : int,
			keep_intermediate_files : bool,
			ignore_substitutions : bool,
			assign_reads_to_all_possible_amplicons : bool,
			suppress_sub_crispresso_plots : bool,
			min_total_reads_per_barcode : int,
			min_reads_per_amplicon_per_cell : int,
			cell_quality_to_analyze : list[str],
			write_h5ad : bool,
			h5ad_output : str,
			h5ad_export_config : dict,
			settings_file : str
		)

	Raises
	------
	Exception
		If required settings are missing or input files do not exist.
	ValueError
		If numeric settings are invalid.
	PermissionError
		If output directories cannot be created or written.
	
	Notes
	-----
	- The function performs side effects:
		* Configures global logging.
		* Creates output directories.
		* Verifies required software availability.
	- Returned values are positional and must be unpacked consistently.
	"""

	# Checking if a settings file is passed into the function
	if len(sys.argv) < 2:
		raise Exception('usage: '+sys.argv[0]+' {settings file}')

	# First argument passed on command line is the settings file #
	settings_file = sys.argv[1]
	settings_file = os.path.abspath(settings_file)
	settings = {}
	
	# Checking if a debug parameter has been provided to the function call #
	logging_level = logging.INFO
	if len(args) > 2 and 'debug' in args[2].lower():
		logging_level=logging.DEBUG
		settings['debug'] = True

	# Settings up logging formatting and parameters
	log_formatter = logging.Formatter("%(asctime)s:%(levelname)s: %(message)s")
	logging.basicConfig(
			level=logging_level,
			format="%(asctime)s: %(levelname)s: %(message)s",
			filename=settings_file+".log",
			filemode='w'
			)
	ch = logging.StreamHandler()
	ch.setFormatter(log_formatter)
	logging.getLogger().addHandler(ch)

	logging.info('Parsing settings file..')

	# Parse the settings file and write to a dictionary {key '\t' value}
	with open(settings_file,'r') as sin:
		for line in sin:
			line = line.strip()
			if line.startswith("#"):
				continue
			if line == "":
				continue
			(key,value) = line.split("\t")
			settings[key] = value

	# Checking for various required settings and raising exceptions for missing values
	if 'r1' not in settings:
		raise Exception('Settings file must contain an entry for r1')
	r1 = settings['r1'].strip()
	# Verifying that all of the read files exist
	for r1_file in r1.split(","):
		if not os.path.isfile(r1_file):
			raise Exception('Input r1 ' + r1_file + ' does not exist')

	# Same as above for r1
	if 'r2' not in settings:
		raise Exception('Settings file must contain an entry for r2')
	r2 = settings['r2'].strip()
	for r2_file in r2.split(","):
		if not os.path.isfile(r2_file):
			raise Exception('Input r2 ' + r2_file + ' does not exist')

	if r1 == r2:
		raise Exception('Input r1 and r2 are the same file.')

	if 'constant1' not in settings:
		raise Exception('Settings file must contain an entry for constant1')
	constant1 = settings['constant1'].strip()
	if 'constant2' not in settings:
		raise Exception('Settings file must contain an entry for constant2')
	constant2 = settings['constant2'].strip()

	allow_barcode_mismatches = False
	if 'allowBarcodeMismatches' in settings:
		allow_barcode_mismatches = True

	# Checking for inclusion and existence of a barcode file
	if 'barcodes' not in settings:
		raise Exception('Settings file must contain an entry for barcodes')
	barcode_file = settings['barcodes'].strip()
	if not os.path.isfile(barcode_file):
		raise Exception('Barcode file ' + barcode_file + ' does not exist')


	# Checking for existence and inclusion of an amplicon file
	if 'amplicons' not in settings:
		raise Exception('Settings file must contain an entry for amplicons')
	amplicon_file = os.path.abspath(settings['amplicons'])
	if not os.path.isfile(amplicon_file):
		raise Exception('Amplicon file does not exist at ' + amplicon_file)

	# Settings primer lookup length
	primer_lookup_len = 18
	if 'primerLookupLen' in settings:
		primer_lookup_len = int(settings['primer_lookup_len'])

	# Settings adapter DNA seq
	adapter_DNA = "TGTCTCTTATACACATCTCCGAGCCCACGAG"
	if 'plamsid_DNA' in settings:
		adapter_DNA = settings['adapter_DNA']

	n_processes = mp.cpu_count()
	if 'processes' in settings and settings['processes'] != 'max':
		n_processes=int(settings['processes'])
	if not n_processes or n_processes < 1:
		raise ValueError("n_processes must be >= 1")

	keep_intermediate_files = False
	if 'keep_intermediate_files' in settings and settings['keep_intermediate_files'].lower() == 'true':
		keep_intermediate_files = True

	ignore_substitutions = False
	if 'ignore_substitutions' in settings and settings['ignore_substitutions'].lower() == 'true':
		ignore_substitutions = True

	assign_reads_to_all_possible_amplicons = False # Default behavior is NOT permissive
	if 'assign_reads_to_all_possible_amplicons' in settings and settings['assign_reads_to_all_possible_amplicons'].lower() == 'true':
		assign_reads_to_all_possible_amplicons = True

	suppress_sub_crispresso_plots = False
	if 'suppress_sub_crispresso_plots' in settings and settings['suppress_sub_crispresso_plots'].lower() == 'true':
		suppress_sub_crispresso_plots = True

	write_h5ad = True
	if 'write_h5ad' in settings:
		write_h5ad = settings['write_h5ad'].strip().lower() == 'true'
	
	# --- normalized cutoffs (use canonical defaults and validate) ---
	try:
		# read user-provided values but store in canonical internal names
		min_total_reads_per_barcode = int(settings.get('min_total_reads_per_barcode', MIN_TOTAL_READS_PER_BARCODE_DEFAULT))
		if min_total_reads_per_barcode < 0:
			raise ValueError("min_total_reads_per_barcode must be >= 0")
	except Exception as e:
		raise ValueError(f"Invalid min_total_reads_per_barcode: {e}")

	try:
		min_reads_per_amplicon_per_cell = int(settings.get('min_reads_per_amplicon_per_cell', MIN_READS_PER_AMPLICON_PER_CELL_DEFAULT))
		if min_reads_per_amplicon_per_cell < 0:
			raise ValueError("min_reads_per_amplicon_per_cell must be >= 0")
	except Exception as e:
		raise ValueError(f"Invalid min_reads_per_amplicon_per_cell: {e}")
	
	# --- parse cell-quality selection (explicit booleans, depth-only terminology) ---

	# Mapping from settings keys -> internal short codes
	_cell_quality_flag_map = {
		"include_high_score_high_depth": "HQ_HI",
		"include_high_score_low_depth":  "HQ_LO",
		"include_low_score_high_depth":  "LQ_HI",
		"include_low_score_low_depth":   "LQ_LO",
	}

	cell_quality_to_analyze = set()

	for key, code in _cell_quality_flag_map.items():
		if key in settings:
			val = settings[key].strip().lower()
			if val == "true":
				cell_quality_to_analyze.add(code)
			elif val == "false":
				pass  # explicitly excluded
			else:
				raise ValueError(
					f"Invalid value for {key}: '{settings[key]}'. "
					"Expected True or False."
				)

	# Default behavior if user specifies none explicitly:
	# Include High_score_High_depth only
	if not cell_quality_to_analyze:
		cell_quality_to_analyze.add("HQ_HI")

	# Deterministic order for downstream logic
	cell_quality_to_analyze = sorted(cell_quality_to_analyze)

	
	# Checking for existence of pregenerated bowtie2 index files if provided
		# .bt2 / .bt21 are the index files generated by bowtie2-build
	bowtie2_index = ""
	if 'bowtie2_index' in settings:
		bowtie2_index = settings['bowtie2_index'].replace(".fa","")
	if 'genome' in settings:
		bowtie2_index = settings['genome'].replace(".fa","")
	if not os.path.isfile(bowtie2_index+".1.bt2") and not os.path.isfile(bowtie2_index+".1.bt2l"):
		raise Exception('bowtie2_index file does not exist at ' + bowtie2_index + ".bt2 or " + bowtie2_index + ".bt2l")

	# Checking for existence and inclusion of an alternate alleles file
	alt_alleles_file = ""
	if 'alt_alleles_file' in settings:
		alt_alleles_file = settings['alt_alleles_file']
		if not os.path.isfile(alt_alleles_file):
			raise Exception('Alt alleles file does not exist at ' + alt_alleles_file)

	# Checking for a user provided output_root
		# if not provided, an suffix is appended onto the settings file name
	output_root = settings_file
	if 'output_root' in settings:
		output_root = settings['output_root']

	h5ad_output = output_root + ".h5ad"
	if 'h5ad_output' in settings:
		h5ad_output = settings['h5ad_output']

	h5ad_zygosity = {}
	for key, default in H5AD_ZYGOSITY_DEFAULTS.items():
		settings_key = f"h5ad_{key}"
		raw_value = settings.get(settings_key, default)
		try:
			h5ad_zygosity[key] = float(raw_value)
		except Exception as e:
			raise ValueError(f"Invalid value for {settings_key}: {raw_value!r} ({e})")

	h5ad_export_config = {
		"analysis_parameters": {
			"zygosity": h5ad_zygosity,
		}
	}

	# Generating amplicon output directory if it does not exist
	amp_file_dir = output_root + ".seq_by_amplicon"
	if not os.path.isdir(amp_file_dir):
		os.mkdir(amp_file_dir)

	crispresso_dir = output_root + ".crispresso"
	if not os.path.isdir(crispresso_dir):
		os.mkdir(crispresso_dir)


	#check software
	#check bowtie2
	try:
		bowtie_result = sb.check_output('bowtie2 --version', stderr=sb.STDOUT,shell=True)
	except Exception:
		raise Exception('Error: bowtie2 is required')

	#check crispresso
	try:
		crispresso_result = sb.check_output('CRISPResso --version', stderr=sb.STDOUT,shell=True)
	except Exception:
		raise Exception('Error: CRISPResso2 is required')


	return (r1, r2, constant1, constant2, allow_barcode_mismatches,barcode_file, amplicon_file, primer_lookup_len, adapter_DNA, amp_file_dir, alt_alleles_file, bowtie2_index, crispresso_dir, output_root, n_processes, keep_intermediate_files, ignore_substitutions, assign_reads_to_all_possible_amplicons, suppress_sub_crispresso_plots, min_total_reads_per_barcode, min_reads_per_amplicon_per_cell, cell_quality_to_analyze, write_h5ad, h5ad_output, h5ad_export_config, settings_file)


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


def reverse_complement(seq):
	"""
	Compute the reverse complement of a DNA sequence.

	Parameters
	----------
	seq : str
		Input nucleotide sequence.

	Returns
	-------
	str
		Reverse complement sequence.

	Notes
	-----
	- Supports characters: A, C, G, T, N, '_', '-'.
	- Output is uppercase.
	"""
	nt_complement=dict({'A':'T','C':'G','G':'C','T':'A','N':'N','_':'_','-':'-'})
	return "".join([nt_complement[c] for c in seq.upper()[-1::-1]])

def add_primer_dict(primer_seq,name,primer_seqs,bad_primer_seqs,add_off_by_1=True):
	"""
	Add a primer sequence and optional one-base variants to lookup tables.

	This function inserts a primer sequence into `primer_seqs` and,
	optionally, generates:
		- All single-nucleotide substitutions
		- All one-base left/right shifts

	Sequences that collide with existing primers are moved to
	`bad_primer_seqs` to prevent ambiguous assignment.

	Parameters
	----------
	primer_seq : str
		Primer sequence.
	name : str
		Amplicon name associated with the primer.
	primer_seqs : dict
		Mapping primer_seq -> amplicon name.
	bad_primer_seqs : dict
		Mapping of ambiguous primer sequences.
	add_off_by_1 : bool, optional
		If True, generate single-base substitutions and shifts.

	Returns
	-------
	int
		Number of primer clashes detected.

	Notes
	-----
	- Modifies `primer_seqs` and `bad_primer_seqs` in place.
	- Collisions result in removal from `primer_seqs`.
	"""
	primer_clashes = 0

	primer_seqs[primer_seq] = name
	if add_off_by_1:
		#first add 1 mismatch
		for i in range(len(primer_seq)):
			for nuc in ['A','C','T','G','N']:
				primer_sub = primer_seq[:i]+nuc+primer_seq[i+1:]
				if primer_sub == primer_seq:
					continue
				if primer_sub in bad_primer_seqs:
					primer_clashes += 1
					continue
				if primer_sub in primer_seqs:
					bad_primer_seqs[primer_sub] = 1
					logging.error("clash between " + primer_sub + "(" + name + ") and " + primer_seqs[primer_sub])
					del primer_seqs[primer_sub]
					primer_clashes += 2
					continue
				primer_seqs[primer_sub] = name

		#next, add shift by 1
		for nuc in ['A','C','T','G','N']:
			primer_shift = nuc + primer_seq[:-1]
			if primer_shift in bad_primer_seqs:
				primer_clashes += 1
				continue
			if primer_shift in primer_seqs:
				bad_primer_seqs[primer_shift] = 1
				del primer_seqs[primer_shift]
				primer_clashes += 2
				continue
			primer_seqs[primer_shift] = name
		for nuc in ['A','C','T','G','N']:
			primer_shift = primer_seq[1:]+nuc
			if primer_shift in bad_primer_seqs:
				primer_clashes += 1
				continue
			if primer_shift in primer_seqs:
				bad_primer_seqs[primer_shift] = 1
				del primer_seqs[primer_shift]
				primer_clashes += 2
				continue
			primer_seqs[primer_shift] = name

	return primer_clashes

def get_valid_barcodes(barcode_file,allow_barcode_mismatches=True):
	"""
	Load valid barcodes and optionally generate single-mismatch variants.

	Each line in `barcode_file` is treated as a canonical barcode.
	If `allow_barcode_mismatches` is True, all single-nucleotide
	substitutions are added as aliases unless ambiguous.

	Parameters
	----------
	barcode_file : str
		Path to file containing one barcode per line.

	allow_barcode_mismatches : bool, optional
		If True, include single-base substitution variants.

	Returns
	-------
	dict[str, str]
		Mapping of observed barcode -> canonical barcode.

	Notes
	-----
	- Ambiguous variants are excluded.
	- Logs number of barcodes and mismatch variants added.
	""" 
	valid_barcodes = {}
	bad_barcodes = {} #barcodes that match with more than one barcode (with 1 mismatch)
	read_barcodes_count = 0
	mismatch_barcodes_count = 0
	barcode_clashes = 0
	with open(barcode_file,'r') as fin:
		for line in fin:
			read_barcodes_count += 1
			barcode = line.strip()
			valid_barcodes[barcode] = barcode
			if allow_barcode_mismatches:
				for i in range(len(barcode)):
					for nuc in ['A','C','T','G','N']:
						new_barcode = barcode[:i]+nuc+barcode[i+1:]
						if new_barcode == barcode:
							continue
						if new_barcode in bad_barcodes:
							barcode_clashes += 1
							continue
						if new_barcode in valid_barcodes:
							bad_barcodes[new_barcode] = 1
							valid_barcodes[new_barcode] = valid_barcodes.pop(new_barcode)
							barcode_clashes += 2
							del valid_barcodes[new_barcode]
							continue
						valid_barcodes[new_barcode] = barcode
						mismatch_barcodes_count += 1


	logging.info('Read ' + str(read_barcodes_count) + ' barcodes from ' + barcode_file)
	logging.info('Added ' + str(mismatch_barcodes_count) + ' barcodes with mismatches')

	return valid_barcodes

def get_primer_seqs(amplicon_file,primer_lookup_len,adapter_DNA,add_off_by_1=True):
	"""
	Build primer lookup table for amplicon assignment.

	Extracts primer sequences from both ends of each amplicon and,
	optionally, generates single-base mismatches and 1-bp shifts
	to increase matching robustness.

	Parameters
	----------
	amplicon_file : str
		Tab-delimited file with columns:
			name    amplicon_sequence    [guide...]

	primer_lookup_len : int
		Number of bases used from each amplicon end.

	adapter_DNA : str
		Adapter sequence used to detect plasmid contamination.

	add_off_by_1 : bool, optional
		If True, generate single-base substitutions and shifts.

	Returns
	-------
	dict[str, str]
		Mapping primer_sequence -> amplicon name.

	Raises
	------
	Exception
		If adapter DNA is detected in amplicon sequence.

	Notes
	-----
	- Primer clashes are logged.
	- Ambiguous primers are excluded.
	"""
	adapter_DNA_rc = reverse_complement(adapter_DNA)
	#add primer sequences to identify reads
	primer_seqs = {}
	bad_primer_seqs = {} #primer seqs that match with more than one amplicon
	primer_clashes = 0
	amps_read_count = 0
	with open(amplicon_file,'r') as fin:
		for line in fin:
			line = line.rstrip()
			line_els = line.split("\t")
			name = line_els[0]
			amplicon = line_els[1]
			#guide = line_els[2]

			amps_read_count += 1

			if adapter_DNA in amplicon or adapter_DNA_rc in amplicon:
				raise Exception("Plasmid DNA is in amplicon " + name + "(" + amplicon + ")")

			primer1 = amplicon[0:primer_lookup_len]
			primer1_rc = reverse_complement(primer1)
			primer2 = amplicon[(-1*primer_lookup_len):]
			primer2_rc = reverse_complement(primer2)
			primer_clashes += add_primer_dict(primer1,name,primer_seqs,bad_primer_seqs,add_off_by_1)
			primer_clashes += add_primer_dict(primer1_rc,name,primer_seqs,bad_primer_seqs,add_off_by_1)
			primer_clashes += add_primer_dict(primer2,name,primer_seqs,bad_primer_seqs,add_off_by_1)
			primer_clashes += add_primer_dict(primer2_rc,name,primer_seqs,bad_primer_seqs,add_off_by_1)

	logging.info('Read ' + str(amps_read_count) + ' amplicons from ' + amplicon_file)
	primer_count = len(primer_seqs.keys())
	mismatch_string = ""
	if add_off_by_1:
		mismatch_string = " with mismatches"
	logging.info('Added ' + str(primer_count) + ' primer seqs' + mismatch_string + ' for aligning reads to amplicons')
	logging.info("Got " + str(primer_clashes) + " primer clashes (off-by-one mismatches)")
	return primer_seqs


class Metrics():
	"""
	Container for counters and aggregation used during FASTQ parsing.

	This small, pickle-safe object is created per-worker (or per-file) during
	FASTQ parsing and then aggregated in the parent process.  It intentionally
	stores only primitive containers (ints and a defaultdict) so that it can be
	transferred safely across multiprocessing boundaries.

	Important contract (invariants)
	------------------------------
	- All counters are **read-based** (i.e., counts of sequencing reads),
	  not barcode- or cell-level aggregates (except where documented).
	- `reads_per_cell` is the sole barcode-level mapping and **must** hold:
		barcode (str) -> integer read count assigned to that barcode.
	  Example: {"AAAGGTT...": 143, "CCGTA...": 12}
	- Percentages reported in the parse log are computed with read-count
	  denominators (e.g. tot_reads), not barcode counts.
	- Metrics are *raw* at this stage: no downstream filtering (depth or
	  quality) should be applied here - filtering occurs later

	Attributes
	----------
	reads_per_cell : collections.defaultdict
		Mapping barcode -> integer read count for that cell (raw, unfiltered).
	tot_reads : int
		Total number of reads processed by this worker.
	has_constant1_count : int
		Number of reads containing constant region 1.
	has_constant2_count : int
		Number of reads containing constant region 2.
	barcodes_valid_count : int
		Number of reads with valid barcodes (after mapping barcode halves).
	barcodes_valid_error_correction_count : int
		Number of reads whose barcode required correction (mismatch mapping).
	long_enough_r1_count : int
		Number of R1 reads that are long-enough for downstream processing.
	no_adapter_read_count : int
		Number of reads that do NOT contain adapter sequences.
	printed_reads : int
		(auxiliary) counter for printed reads, not used for logic.

	Notes
	-----
	Keep this class minimal. The parent process aggregates worker Metrics by
	calling `gather_metrics()` and then uses `create_log_str()` to emit a
	human-readable summary immediately after parsing finishes.
	""" 
	def __init__(self):
		# A dictionary containing a map of (str) barcode -> (int) read count
		self.reads_per_cell = defaultdict(int)
		self.tot_reads = 0
		self.has_constant1_count = 0
		self.has_constant2_count = 0
		self.barcodes_valid_count = 0
		self.barcodes_valid_error_correction_count = 0
		self.long_enough_r1_count = 0
		self.no_adapter_read_count = 0
		self.printed_reads = 0

	def num_cells(self):
		return len(self.reads_per_cell)

	def gather_metrics(self, other):
		"""
		Add metrics from another Metrics object into this one.

		Parameters
		----------
		other : Metrics
		Another Metrics instance whose counters will be added into this one.
		"""
		self.tot_reads += other.tot_reads
		self.has_constant1_count += other.has_constant1_count
		self.has_constant2_count += other.has_constant2_count
		self.barcodes_valid_count += other.barcodes_valid_count
		self.barcodes_valid_error_correction_count += other.barcodes_valid_error_correction_count
		self.long_enough_r1_count += other.long_enough_r1_count
		self.no_adapter_read_count += other.no_adapter_read_count
		self.printed_reads += other.printed_reads
		for key in other.reads_per_cell:
			self.reads_per_cell[key] += other.reads_per_cell[key]
			
	def create_log_str(self):
		"""
		Create a human-readable summary of parse-stage metrics.

		Notes (explicit read-based semantics)
		------------------------------------
		- long_enough_r1_pct is defined as:
			long_enough_r1_count / tot_reads
		i.e. percent of *total reads* that had sufficiently-long R1s.
		- no_adapter_rc_pct is defined as:
			no_adapter_read_count / long_enough_r1_count
		i.e. percent of long-enough reads that did not contain adapters.
		- The only barcode-level aggregation reported is the number of distinct
		barcodes observed (self.num_cells()).

		Returns
		-------
		str
			Multi-line human-readable summary string suitable for logging.
		""" 
		def safe_pct(part: int, whole: int) -> str:
			if whole <= 0:
				return "0.00"
			return f"{100.0 * part / whole:.2f}"
		long_enough_r1_pct = safe_pct(self.long_enough_r1_count, self.tot_reads)
		no_adapter_rc_pct = safe_pct(self.no_adapter_read_count, self.long_enough_r1_count)

		info_str = "Read "+str(self.tot_reads)+" reads\n" + \
	"  Of these, "+str(self.has_constant1_count) + " have the constant region 1\n"+ \
	"    Of these, "+str(self.has_constant2_count) + " have the constant region 2\n"+ \
	"      Of these, "+str(self.barcodes_valid_count) + " have valid cell barcodes (" + str(self.barcodes_valid_error_correction_count) + " with error correction)\n"+ \
	"        Of these, "+str(self.long_enough_r1_count) + " ("+str(long_enough_r1_pct)+"%) have sufficiently-long R1's\n"+ \
	"          Of these, "+ str(self.no_adapter_read_count) + " ("+str(no_adapter_rc_pct)+"%) did not contain adapter sequences\n"\
	"             which are assigned to "+str(self.num_cells()) + " cells\n"  
		return info_str
		
	  
	
def parse_fq_file_pair(args):
	"""
	Parse a paired-end FASTQ file pair and extract reads that can be assigned to cells.

	This function is intended to be executed in a worker process (via
	multiprocessing.Pool). It reads paired FASTQ records in lockstep, scans for
	primer/constant sequences, extracts barcodes, applies basic length / adapter
	filters, writes per-worker parsed FASTQ outputs, and accumulates a Metrics
	object describing the parsing results.

	Parameters
	----------
	args : tuple
		Tuple containing:
		  - r1_path (str): path to R1 FASTQ
		  - r2_path (str): path to R2 FASTQ
		  - fq_index (int): ordinal index used to name per-job outputs
		  - valid_barcodes (dict): mapping of barcode halves -> canonical barcode
		  - constant1 (str): first constant primer sequence to search for
		  - constant2 (str): second constant primer sequence to search for
		  - constant1_len (int): length of constant1
		  - constant2_len (int): length of constant2
		  - adapter_DNA (str): adapter sequence to screen out
		  - adapter_DNA_rc (str): reverse-complement adapter sequence
		  - output_root (str): root path for per-job outputs

	Returns
	-------
	tuple
		(out1_name, out2_name, bam_name, metrics)
		  - out1_name (str): path to worker R1 parsed FASTQ
		  - out2_name (str): path to worker R2 parsed FASTQ
		  - bam_name (str): placeholder name for downstream alignment output
		  - metrics (Metrics): populated Metrics instance for this input pair

	Invariants / notes
	------------------
	- This function MUST return a Metrics instance (not a raw dict) so the
	  parent process can aggregate with `gather_metrics()`.
	- `metrics.reads_per_cell` is incremented as `metrics.reads_per_cell[barcode] += 1`
	  for every read that survived filtering and was assigned to `barcode`.
	- No global state is modified: all outputs are file-based and Metrics is
	  returned to the parent for aggregation.
	"""
	
	r1_path,r2_path,fq_index,valid_barcodes,constant1,constant2,constant1_len, constant2_len,adapter_DNA, adapter_DNA_rc,output_root = args
	
   
	metrics = Metrics()
	
	out1_name = build_stage_filename(
		STAGE_PARSE,
		f"parsed_worker{fq_index}",
		read="r1",
		ext="fastq",
		output_root = output_root
	)

	out2_name = build_stage_filename(
			STAGE_PARSE,
			f"parsed_worker{fq_index}",
			read="r2",
			ext="fastq",
			output_root = output_root
		)

	bam_name = build_stage_filename(
			STAGE_ALIGN,
			f"align_worker{fq_index}",
			ext="bam",
			output_root = output_root
		)

	
	out1 = dnaio.open(out1_name, mode = 'w', fileformat = 'fastq')
	out2 = dnaio.open(out2_name, mode = 'w', fileformat = 'fastq')
 
	with dnaio.open(r1_path, fileformat = 'fastq') as f1, dnaio.open(r2_path, fileformat = 'fastq') as f2:
		for rec1, rec2 in zip(f1, f2):
			metrics.tot_reads += 1
			
			info1 = rec1.name.strip()
			seq1 = str(rec1.sequence)
			plus1 = '+'
			qual1 = rec1.qualities
			
			info2 = rec2.name.strip()
			seq2 = str(rec2.sequence)
			plus2 = '+'
			qual2 = rec2.qualities
			
			const1_loc = seq1.find(constant1)
			const2_loc = seq1.find(constant2)
			if const1_loc < 0:
				continue
			metrics.has_constant1_count += 1
			if const2_loc < 0:
				continue
			metrics.has_constant2_count += 1

			barcode_part1 = seq1[0:9]
			barcode_part2 = seq1[const1_loc+constant1_len:const1_loc+constant1_len+9]

			if barcode_part1 not in valid_barcodes or barcode_part2 not in valid_barcodes:
				continue
			metrics.barcodes_valid_count += 1
			
			barcode = valid_barcodes[barcode_part1] + valid_barcodes[barcode_part2]
			if barcode != barcode_part1 + barcode_part2:
				metrics.barcodes_valid_error_correction_count += 1
				
			last_r1_bit = seq1[const2_loc+constant2_len:]
			last_r1_qual = qual1[const2_loc+constant2_len:]
			
			if len(last_r1_bit) < 20 or len(barcode) < 18:
				continue
			
			metrics.long_enough_r1_count += 1
			
			if adapter_DNA in seq1 or adapter_DNA_rc in seq1 or adapter_DNA in seq2 or adapter_DNA_rc in seq2:
				continue
			metrics.no_adapter_read_count += 1
			
			
			info1_first_bit = info1.split(' ')[0]
			info2_first_bit = info2.split(' ')[0]


			
			out1.write(dnaio.SequenceRecord(
					name = f"{info1_first_bit}:{barcode}",
					sequence= last_r1_bit,
					qualities = last_r1_qual
				))
			out2.write(dnaio.SequenceRecord(
				name = f"{info2_first_bit}:{barcode}",
				sequence= seq2,
				qualities = qual2
			))
			# assign read to barcode; metrics.reads_per_cell maps barcode -> read count
			metrics.reads_per_cell[barcode] += 1
			
	out1.close()
	out2.close()
	logging.info("Read %d reads from %s and %s", metrics.tot_reads, r1_path, r2_path)

	return (out1_name, out2_name, bam_name, metrics)

def run_alignment(args):
	"""
	Execute bowtie2 alignment for a paired FASTQ file pair.

	Parameters
	----------
	args : tuple
		(
			r1_path : str,
			r2_path : str,
			bowtie2_index : str,
			threads : int,
			out_name : str
		)

	Returns
	-------
	None

	Notes
	-----
	- Runs bowtie2 piped into samtools view to generate a BAM file.
	- Writes bowtie2 stderr to `<out_name>.bowtie2.log`.
	- Uses `run_command` internally.
	- Intended for multiprocessing execution.
	"""
	r1_path, r2_path, bowtie2_index, threads, out_name = args
	align_cmd = f'bowtie2 -x {bowtie2_index} -p {threads} -1 {r1_path} -2 {r2_path} 2>>{out_name}.bowtie2.log | samtools view -bS - > {out_name}'
	run_command(align_cmd)

def parse_and_align_reads(r1_fastqs,r2_fastqs,constant1,constant2,
						  output_root,barcode_file,allow_barcode_mismatches,
						  adapter_DNA,bowtie2_index,n_processes,keep_intermediate_files=False):
	"""
	Parse input FASTQs (possibly in parallel), align parsed reads, and produce a name-sorted BAM
	and a cell-count mapping.

	High-level behavior
	-------------------
	1. If parse outputs already exist (info, aligned BAM, and cell_count file),
	   return them (fast path).
	2. Otherwise:
	   - Build valid barcode lookup table.
	   - Spawn parsing workers (one per FASTQ pair); each worker writes per-job parsed FASTQs
		 and returns a Metrics object describing its file.
	   - Aggregate worker Metrics using `gather_metrics()` into a single merged_metrics.
	   - Emit a human-readable parsing summary (merged_metrics.create_log_str()) immediately.
	   - Run alignment on per-job parsed FASTQs (bowtie2 + samtools).
	   - Merge/sort BAMs, write cell-count file (`<output_root>.parseReads.cellCount.txt`),
		 and return `(aligned_bam, merged_metrics.reads_per_cell)`.

	Parameters
	----------
	r1_fastqs : str
		Comma-separated list of R1 FASTQ paths.
	r2_fastqs : str
		Comma-separated list of R2 FASTQ paths (must match r1_fastqs in length).
	constant1, constant2 : str
		Primer/constant sequences used to locate barcodes.
	output_root : str
		Root path/prefix for pipeline outputs.
	barcode_file : str
		Path to barcode lookup file; read by `get_valid_barcodes()`.
	allow_barcode_mismatches : bool
		Whether to accept barcodes with a small number of mismatches.
	adapter_DNA : str
		Adapter sequence used to screen contaminating reads.
	bowtie2_index : str
		Prefix of bowtie2 index for alignment.
	n_processes : int
		Number of worker processes to use for parsing/parallelizable steps.
	keep_intermediate_files : bool, optional
		If True, do not delete intermediate per-job files.

	Returns
	-------
	tuple
		(aligned_bam_path (str), reads_per_cell (dict))
		  - reads_per_cell is a mapping barcode -> integer read count (unfiltered).

	Invariants / notes
	------------------
	- The human-readable parsing summary (create_log_str) is printed and logged
	  immediately after parsing and before alignment. This ensures the summary
	  appears even if later stages (alignment/merging) fail.
	- `reads_per_cell` returned is raw (no filtering). Downstream filtering (depth,
	  per-amplicon thresholds) is performed by `split_reads_by_amplicon` and
	  functions that consume the returned `reads_per_cell`.
	"""
	info_file = output_root+".parseReads.info.txt"
	#aligned_bam = output_root+".alignReads.bam"
	cell_file = output_root + ".parseReads.cellCount.txt"
	aligned_bam = build_stage_filename(
		 STAGE_ALIGN,
		 "align_merged",
		 ext = "bam",
		 output_root = output_root
	)
	
	
	if os.path.isfile(info_file) and os.path.isfile(aligned_bam) and os.path.isfile(cell_file):
		reads_per_cell = {}
		with open(cell_file, 'r') as fin:
			for line in fin:
				reads_per_cell[line.split('\t')[0]] = int(line.split('\t')[1].strip())
		logging.info ("Finished parsing reads")
		return (aligned_bam, reads_per_cell)
	
	# valid_barcodes: dict of valid cell barcodes
	valid_barcodes = get_valid_barcodes(barcode_file,allow_barcode_mismatches)

	logging.info("Parsing reads..")
	
	adapter_DNA_rc = reverse_complement(adapter_DNA)
	constant1_len = len(constant1)
	constant2_len = len(constant2)

	r1_fastq_arr = r1_fastqs.split(",")
	r2_fastq_arr = r2_fastqs.split(",")
	if len(r1_fastq_arr) != len(r2_fastq_arr):
		raise Exception("Incorrect number of fastq paired files")

	# Create parsed fq files
	start_fq_parsing = time.time()
	with mp.Pool(n_processes) as pool:
		parsed_results = pool.map(
			parse_fq_file_pair, [(r1_fastq_arr[i], r2_fastq_arr[i],
								  i, valid_barcodes,
								  constant1, constant2,
								  constant1_len, constant2_len,
								  adapter_DNA, adapter_DNA_rc,
								  output_root) for i in range(len(r1_fastq_arr))]
		)
	# Merge parsed metrics
	merged_metrics = Metrics()
	for res in parsed_results:
		merged_metrics.gather_metrics(res[3])

	logging.info("%s", merged_metrics.create_log_str())

	end_fq_parsing = time.time() - start_fq_parsing
	logging.info("FQ parsing finished in %.3f seconds", end_fq_parsing)        
	# Calculate cores for bowtie2
	thread_count = int(os.environ.get("SLURM_CPUS_PER_TASK", n_processes))
	threads_per_job = max(1, thread_count // len(r1_fastq_arr))
	
	logging.info(f"Running bowtie2 with {threads_per_job} threads per job on {len(r1_fastq_arr)} jobs")
	start_bowtie_alignment = time.time()
	# Build alignment commands
	with mp.Pool(len(r1_fastq_arr)) as pool:
		pool.map(
			run_alignment, [(parsed_results[i][0], parsed_results[i][1],
							  bowtie2_index, threads_per_job,
							  parsed_results[i][2]) for i in range(len(parsed_results))]
		)
	end_bowtie_alignment = time.time() - start_bowtie_alignment
	logging.info("Bowtie2 alignment finished in %.3f seconds", end_bowtie_alignment)

	# Merge BAM files
	bam_threads = min(n_processes, 12)
	#bam_output =  output_root+".alignReads.bam"

	# Check for multiple parsed_results
	if len(parsed_results) > 1: 
		input_bams = " ".join(result[2] for result in parsed_results)
		#inter_bam = output_root+".merged.bam"
		inter_bam = build_stage_filename(STAGE_ALIGN, "intermediate_merged", ext="bam", output_root=output_root)
		start_bam_cat = time.time()
		run_command(f"samtools cat -o {inter_bam} {input_bams}")
		end_bam_cat = time.time() - start_bam_cat
		logging.info("BAM cat ended in %.3f seconds", end_bam_cat)
		logging.info("Inside parsed_results > 1")
	if len(parsed_results) == 1:
		inter_bam = parsed_results[0][2]
		logging.info("Inside parsed_results == 1")

	bam_threads = min(12, n_processes)
	sort_bam_cmd = f'samtools sort -n -@ {bam_threads} -o {aligned_bam} {inter_bam}'
		 
	start_bam_sort = time.time() 
	run_command(sort_bam_cmd)
	end_bam_sort = time.time() - start_bam_sort
	logging.info("BAM sort finished in %.3f seconds", end_bam_sort)


	# Optionally remove per worker parsed fq files
	if not keep_intermediate_files:
		for res in parsed_results:
			#print("Got to 1983")
			parsed_r1 = res[0]
			parsed_r2 = res[1]
			#print(f"{parsed_r1=}\n{parsed_r2=}")
			try:
				safe_remove(parsed_r1)
				safe_remove(parsed_r2)
			except Exception:
				logging.warning("Failed to remove parsed fastq for workers (continuing)")

	try:
		if inter_bam != aligned_bam:
			safe_remove(inter_bam)
	except Exception:
		logging.warning("Failed to remove intermediate BAM %s (continuing)", inter_bam)

	
	# Write out cell count file
	with open(cell_file, 'w') as fout:
		for cell in merged_metrics.reads_per_cell:
			fout.write(f"{cell}\t{merged_metrics.reads_per_cell[cell]}\n") 
	
	return (aligned_bam, merged_metrics.reads_per_cell)    



def alignment_end(pos, cigar):
	ref_len = 0
	for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
		length = int(length)
		if op in ('M', 'D', 'N', '=', 'X'):
			ref_len += length
	return pos + ref_len - 1



def split_reads_by_amplicon(aligned_bam, output_root,amplicon_file,alt_alleles_file,primer_lookup_len,amp_file_dir,bowtie2_index,adapter_DNA,n_processes,keep_intermediate_files, reads_per_cell, min_total_reads_per_barcode, assign_reads_to_all_possible_amplicons=False):
	"""
	Split reads from a name-sorted aligned BAM into per-amplicon FASTQ files.

	Behavior summary
	----------------
	- Build primer lookup tables and optionally align amplicons to the genome.
	- Iterate through the name-sorted BAM and assign each read pair to one or
	  more amplicons based on primer matches and alignment position.
	- Write per-amplicon R1/R2 FASTQ files and an amplicon info file used to
	  accelerate re-runs.
	- Uses `reads_per_cell` (barcode -> read count) as input for optional filtering
	  and to prioritize which barcodes to extract.

	Parameters
	----------
	aligned_bam : str
		Path to name-sorted aligned BAM.
	output_root : str
		Base path for outputs.
	amplicon_file : str
		Tab-delimited file: name\tamplicon_sequence\tguide
	alt_alleles_file : str or None
		Optional alternate allele sequences file.
	primer_lookup_len : int
		Number of bases to use from read ends for primer matching.
	amp_file_dir : str
		Directory to write per-amplicon FASTQ files.
	bowtie2_index : str
		Prefix to bowtie2 genome index.
	adapter_DNA : str
		Adapter sequence to screen against.
	n_processes : int
		Number of worker processes where parallelism is applicable.
	keep_intermediate_files : bool
		Keep per-amplicon intermediate files.
	reads_per_cell : dict
		Mapping barcode -> integer read count (raw, unfiltered).
	min_total_reads_per_barcode : int
		Minimum total reads to consider a barcode in downstream steps.
	assign_reads_to_all_possible_amplicons : bool
		If True, assign a read to every amplicon it plausibly matches;
		otherwise pick a single best amplicon.

	Returns
	-------
	tuple
		(amplicon_names, amplicon_information, amplicon_info_file)
		- amplicon_information: a structure describing written amplicon FASTQs
		- amplicon_info_file: path to the amplicon info JSON/text used to speed re-runs

	Notes
	-----
	- This function assumes `aligned_bam` is name-sorted (paired reads one after another).
	- Assignment heuristics are documented in the code near the primer match logic;
	  keep those comments in sync with this docstring if you change heuristics.
	"""
	info_file = output_root+".splitReads.ampliconInfo.txt"
	if os.path.isfile(info_file):
		with open(info_file,'r') as fin:
			head = fin.readline().strip()
			head_els = head.split("\t")
			amplicon_names = []
			amplicon_information = {}
			for line in fin:
				line_els = line.strip().split("\t")
				amp_info = dict(zip(head_els,line_els))
				amp_name = line_els[0]
				if 'reads_r1_file' not in amp_info or amp_info['reads_r1_file'] in ('', 'NA'):
					#amp_info['reads_r1_file'] = os.path.join(amp_file_dir, amp_name + '.r1.fq')
					amp_info['reads_r1_file'] = build_stage_filename(
				stage=STAGE_SPLIT,
				tag="reads_all_cells",
				amplicon=amp_name,
				read="r1",
				ext="fq",
				output_root=amp_file_dir,
			)
				if 'reads_r2_file' not in amp_info or amp_info['reads_r2_file'] in ('', 'NA'):
					#amp_info['reads_r2_file'] = os.path.join(amp_file_dir, amp_name + '.r2.fq')
					amp_info['reads_r2_file'] = build_stage_filename(
				stage=STAGE_SPLIT,
				tag="reads_all_cells",
				amplicon=amp_name,
				read="r2",
				ext="fq",
				output_root=amp_file_dir,
			)
				amplicon_information[amp_name] = amp_info
				amplicon_names.append(amp_name)

		logging.info ("Finished splitting reads")
		return amplicon_names,amplicon_information,info_file

	logging.info("Splitting reads to amplicons..")


	amplicon_names = [] #ordered list of amplicons
	amplicon_information = {}#a dict to store all amplicons information. Keys are amplicon names, values are dicts of values for those amplicons
	#first, find the locations of amplicons
	#create a fastq from the reads
	amplicon_fasta_file = output_root + ".amplicons.fa"
	with open(amplicon_file,'r') as amps_in, open(amplicon_fasta_file,'w') as amps_out:
		for line in amps_in:
			line_els = line.strip().split("\t")
			amp_name = line_els[0]
			amp_seqs = line_els[1]
			guide_seq = ""
			if len(line_els) > 2:
				guide_seq = line_els[2]
			input_ref_allele_counts = "2"
			if len(line_els) > 3:
				input_ref_allele_counts = line_els[3]

			amplicon_information[amp_name] = {
					'name':line_els[0],
					'amp_seqs':amp_seqs,
					'guide_seq':guide_seq,
					'input_amp_seqs':amp_seqs,
					'input_alternate_allele_seqs':'NA',
					'input_ref_allele_counts':input_ref_allele_counts,
					'aln_chr':'NA',
					'aln_start':'NA',
					'aln_end':'NA',
					'aln_score':'NA',
					'secondary_aln_chr':'NA',
					'secondary_aln_start':'NA',
					'secondary_aln_end':'NA',
					'secondary_aln_score':'NA',
					'aln_count':'0',
					'reads_r1_file':'NA',
					'reads_r2_file':'NA',
					}
			first_amp_seq = amp_seqs.split(",")[0]
			amps_out.write(">"+amp_name+"\n"+first_amp_seq+"\n")
			amplicon_names.append(amp_name)

	#parse alternate alleles
	if alt_alleles_file != "":
		read_alleles_count = 0
		with open(alt_alleles_file,'r') as fin:
			head = fin.readline()
			for line in fin:
				line_els = line.strip().split("\t")
				this_amp_name = line_els[0]
				this_given_alleles = line_els[1]
				amplicon_information[this_amp_name]['input_alternate_allele_seqs'] = line_els[2]
				amplicon_information[this_amp_name]['amp_seqs'] = line_els[2]
				seen_amp_seqs = {}
				for amp_seq in line_els[2].split(","):
					if amp_seq.lower() in seen_amp_seqs:
						logging.error('Amplicon is present twice for ' + this_amp_name + ' from alternate alleles file ' + alt_alleles_file + ': ' + amp_seq)
				read_alleles_count += 1
		logging.info('Read ' + str(read_alleles_count) + ' alternate alleles')

	#align amplicons to genome
	bowtie2_log = output_root+".amplicons.alignReads.log"
	aligned_amps_file = output_root + ".amplicons.aligned.sam"
	align_command = 'bowtie2 -k 2 -x %s -p %s -f -U %s 2>>%s > %s'%(bowtie2_index,n_processes,amplicon_fasta_file,bowtie2_log,aligned_amps_file)
	with open(bowtie2_log,'w') as bt2log:
		bt2log.write(align_command)
	run_command(align_command)

	#genome_amplicon_locs = {}#chr,start -> amplicon
	start_amplicon_locs = {}# chr,start -> amplicon
	end_amplicon_locs = {}# chr,end -> amplicon
	amplicon_alignment_details_file = output_root + ".amplicons.details.txt"
	amplicon_alignment_count = 0
	sec_count = 0
	with open(aligned_amps_file,'r') as fin:
		for line in fin:
			if line.startswith("@"):
				continue
			line_els = line.split("\t")
			amp_name = line_els[0]
			line_mapq = line_els[4]
			line_unmapped = int(line_els[1]) & 0x4
			line_secondary = int(line_els[1]) & 0x100
			line_chr = line_els[2]
			line_start = int(line_els[3]) - 1
			line_end = line_start + len(line_els[9])
			if not line_unmapped and not line_secondary:
				amplicon_information[amp_name]['aln_chr'] = line_chr
				amplicon_information[amp_name]['aln_start'] = str(line_start)
				amplicon_information[amp_name]['aln_end'] = str(line_end)
				amplicon_information[amp_name]['aln_score'] = line_mapq
				amplicon_alignment_count += 1

				start_key = line_chr+":"+str(line_start)
				end_key = line_chr+":"+str(line_end - 1)
				
				if start_key in start_amplicon_locs:
					raise Exception('Error: amplicons ' + amp_name + ' and ' + start_amplicon_locs[start_key] + ' align to the same location (' + start_key + ')')

				if end_key in end_amplicon_locs:
					raise Exception('Error: amplicons ' + amp_name + ' and ' + end_amplicon_locs[end_key] + ' align to the same location (' + end_key + ')')

				# Create a key range buffer to catch reads that do not perfectly start | end at the primer sites
				for i in range(line_start-3, line_start+4):
					start_amplicon_locs[f"{line_chr}:{str(i)}"] = amp_name
				for i in range(line_end-3, line_end+4):
					end_amplicon_locs[f"{line_chr}:{str(i)}"] = amp_name
			   
			if line_secondary:
				sec_count += 1
				amplicon_information[amp_name]['secondary_aln_chr'] = line_chr
				amplicon_information[amp_name]['secondary_aln_start'] = str(line_start)
				amplicon_information[amp_name]['secondary_aln_end'] = str(line_end)
				amplicon_information[amp_name]['secondary_aln_score'] = line_mapq
	logging.info('Found alignments for ' + str(amplicon_alignment_count) + ' amplicons')

	#primer_seqs: dict of primers -> amplicon for alignment
	primer_seqs = get_primer_seqs(amplicon_file,primer_lookup_len,adapter_DNA)

	amp_filehandles = {} #dict of amplicon => filehandles (tuple of (r1 and r2))
	amp_filenames = [] #array of names of created files
	unidentified_out1_name = os.path.join(amp_file_dir,"unidentified.r1.fq")
	unidentified_out2_name = os.path.join(amp_file_dir,"unidentified.r2.fq")
	unidentified_out1 = open(unidentified_out1_name,'w')
	unidentified_out2 = open(unidentified_out2_name,'w')

	logging.info('Splitting reads to amplicon files..')
	amplicon_count = defaultdict(int) #amplicon-> read count
	aln_count = defaultdict(int) #amp_r1,amp_r2,aln_amp_r1... -> read count
	unaln_barcode_count = defaultdict(int) #cell barcode -> count unaligned
	aln_barcode_count = defaultdict(int) # cell barcode -> count aligned
	tot_barcode_count = set()
	barcode_count_dict = defaultdict(int) #barcode -> barcode - > aligned / unaligned / percentage
	tot_reads_count = 0
	aln_reads_count = 0 #not secondary or not unaligned
	id_reads_count = 0
	chimeric_reads_count = 0 #primer1 != primer2 
	unidentified_reads_count = 0 #primers not match location or primers not match each other
	aligned_other_loc_reads_count = 0 # reads that have primers for an amplicon but are aligned to another location
	partial_alignment_rescue_count = 0 # primers agreed and one alignment-side check supported the amplicon with no contradiction
	unmapped_reads_count = 0 #unmapped by bowtie
	bam_iter = get_command_output('samtools view %s'%(aligned_bam))#read in the aligned bam file
	count = 0

	# set to keep track of failed cells
	failing_barcode = set()

	# The iterator keeps going past the end of samtools output
	for line1 in bam_iter:
		line2 = next(bam_iter)
		count += 1

		line1_els = line1.split("\t")
		line2_els = line2.split("\t")

		info1 = line1_els[0]
		info2 = line2_els[0]

		# Adding a break statement if the line is empty #
		if line1_els[0] == "":
			break

		seq1 = line1_els[9]
		seq2 = line2_els[9]
		qual1 = line1_els[10]
		qual2 = line2_els[10]

		tot_reads_count += 1

		line1_mapq = line1_els[5]
		line1_unmapped = int(line1_els[1]) & 0x4
		line1_rc = int(line1_els[1]) & 0x10
		line1_secondary = int(line1_els[1]) & 0x100

		line2_mapq = line2_els[5]
		line2_unmapped = int(line2_els[1]) & 0x4
		line2_rc = int(line2_els[1]) & 0x10
		line2_secondary = int(line2_els[1]) & 0x100
		line1_cigar = line1_els[5]

		line1_chr = line1_els[2]
		line1_start = int(line1_els[3])-1
		line2_chr = line2_els[2]
		line2_start = int(line2_els[3])-1
		line2_cigar = line2_els[5]

		if line1_secondary or line2_secondary:
			continue

		info1_first_bit = info1.split(' ')[0]
		info2_first_bit = info2.split(' ')[0]
		if info1_first_bit != info2_first_bit:
			raise Exception('Error, bam is not read name-sorted. Please sort using samtools sort -n -o out.bam {input.bam}.')
		barcode = info1_first_bit.split(":")[-1]
	
		# If barcode does not have enough reads, pass
		if reads_per_cell[barcode] < min_total_reads_per_barcode:
			failing_barcode.add(barcode)
			continue

		tot_barcode_count.add(barcode)
   
		if not(line1_unmapped or line2_unmapped):
			aln_reads_count += 1
	
		if barcode not in unaln_barcode_count.keys():
			unaln_barcode_count[barcode] = 0
		if barcode not in aln_barcode_count.keys():
			aln_barcode_count[barcode] = 0 
			
		seq1_to_write = seq1
		primer1 = seq1[0:primer_lookup_len]
		key1 = line1_chr + ":" + str(line1_start)
		if line1_rc:
			primer1 = seq1[(-1*primer_lookup_len):]
			key1 = line1_chr + ":" + str(alignment_end(line1_start, line1_cigar))
			seq1_to_write = reverse_complement(seq1)

		seq2_to_write = seq2
		primer2 = seq2[0:primer_lookup_len]
		key2 = line2_chr + ":" + str(line2_start)
		if line2_rc:
			primer2 = seq2[(-1*primer_lookup_len):]#after alignment, all reads are put on the forward strand so we don't need to reverse complement them
			key2 = line2_chr + ":" + str(alignment_end(line2_start, line2_cigar))
			seq2_to_write = reverse_complement(seq2)


		amp1 = "NA"
		amp2 = "NA"
		if primer1 in primer_seqs:
			amp1 = primer_seqs[primer1]
		if primer2 in primer_seqs:
			amp2 = primer_seqs[primer2]
		amp1_aln = "NA"
		amp2_aln = "NA"
		if key1 in start_amplicon_locs: # start_amplicon_locs[ChrN:10000000] = amp_name
			amp1_aln = start_amplicon_locs[key1]
		if key2 in end_amplicon_locs:
			amp2_aln = end_amplicon_locs[key2]

		amp_same_key = ""
		if assign_reads_to_all_possible_amplicons: # write all possible amplicons where this read could match - from the primer matching or the alignment
			if amp1 != "NA" and amp2 != "NA": # This block is focused on logging reads with unexpected primer seq combinations
				if amp1 == amp2: # primer seqs in agreement
					pass
				elif amp1_aln != "NA" and amp1_aln != amp1: # Amp determined by primer disagrees with amp from alignment
					aligned_other_loc_reads_count += 1
					unaln_barcode_count[barcode] += 1
				elif amp2_aln != "NA" and amp2_aln != amp2: 
					aligned_other_loc_reads_count += 1
					unaln_barcode_count[barcode] += 1
				elif amp1 != amp2: # Amplicon determined by primer sequences disagrees
					chimeric_reads_count += 1
					unaln_barcode_count[barcode] += 1
				else:
					unidentified_reads_count += 1
					unaln_barcode_count[barcode] += 1
			else:
				unidentified_reads_count += 1
				unaln_barcode_count[barcode] += 1

			candidate_amps = set([amp1,amp2,amp1_aln,amp2_aln]) # these are all possible amplicons this read could belong to
			candidate_amps.discard('NA') # We don't want to write out to an "NA" file
			for amp in candidate_amps:
				if amp not in amp_filehandles:
					# Construct filehandle to write reads for {amp}
					amp_filename_r1 = build_stage_filename(
						stage = STAGE_SPLIT,
						tag = "reads_all_cells",
						amplicon = amp,
						read = "r1",
						ext = "fq",
						output_root = amp_file_dir
					)
					
					amp_filename_r2 = build_stage_filename(
						stage = STAGE_SPLIT,
						tag = "reads_all_cells",
						amplicon = amp,
						read = "r2",
						ext = "fq",
						output_root = amp_file_dir
					)

					#print(f"Line 2329\n{amp}\t{amp_filename_r1}\t{amp_filename_r2}\n")

					#append a tuple of R1 and R2
					amp_filenames.append((amp_filename_r1,amp_filename_r2))
					amplicon_information[amp]['reads_r1_file'] = amp_filename_r1 + '.gz'
					amplicon_information[amp]['reads_r2_file'] = amp_filename_r2 + '.gz'

					fh_r1 = open(amp_filename_r1,'w')
					fh_r2 = open(amp_filename_r2,'w')
					amp_filehandles[amp] = (fh_r1,fh_r2)


				amp_filehandles[amp][0].write("@%s\n%s\n%s\n%s\n"%(info1,seq1_to_write,"+",qual1))
				amp_filehandles[amp][1].write("@%s\n%s\n%s\n%s\n"%(info2,seq2_to_write,"+",qual2))
				id_reads_count += 1
				amplicon_count[amp] += 1
				aln_barcode_count[barcode] += 1
		
		else: #if not assign to all possible amplicons, only write if the primers and the alignment match
			# Require primer agreement and at least one non-conflicting alignment-side check.
			# Some valid read pairs only recover one amplicon boundary from the genomic alignment;
			# previously those were discarded even when both primer calls agreed on the target.
			accepted_amplicon = None
			if amp1 == amp2 and amp1 != "NA":
				alignment_calls = [x for x in (amp1_aln, amp2_aln) if x != "NA"]
				if alignment_calls and all(x == amp1 for x in alignment_calls):
					accepted_amplicon = amp1
					if len(alignment_calls) < 2:
						partial_alignment_rescue_count += 1

			if accepted_amplicon is not None:
				if accepted_amplicon not in amp_filehandles:
					amp_filename_r1 = build_stage_filename(
						stage=STAGE_SPLIT,
						tag="reads_all_cells",
						amplicon=accepted_amplicon,
						read="r1",
						ext="fq",
						output_root=amp_file_dir,
					)
					amp_filename_r2 = build_stage_filename(
						stage=STAGE_SPLIT,
						tag="reads_all_cells",
						amplicon=accepted_amplicon,
						read="r2",
						ext="fq",
						output_root=amp_file_dir,
					)

					amp_filenames.append((amp_filename_r1, amp_filename_r2))
					amplicon_information[accepted_amplicon]['reads_r1_file'] = amp_filename_r1 + '.gz'
					amplicon_information[accepted_amplicon]['reads_r2_file'] = amp_filename_r2 + '.gz'

					fh_r1 = open(amp_filename_r1,'w')
					fh_r2 = open(amp_filename_r2,'w')
					amp_filehandles[accepted_amplicon] = (fh_r1,fh_r2)

				amp_filehandles[accepted_amplicon][0].write("@%s\n%s\n%s\n%s\n"%(info1,seq1_to_write,"+",qual1))
				amp_filehandles[accepted_amplicon][1].write("@%s\n%s\n%s\n%s\n"%(info2,seq2_to_write,"+",qual2))
				id_reads_count += 1
				amplicon_count[accepted_amplicon] += 1
			else:
				unidentified_out1.write("@%s\n%s\n%s\n%s\n"%(info1,seq1_to_write,"+\t"+"\t".join([amp1,amp2,amp1_aln,amp2_aln]),qual1))
				unidentified_out2.write("@%s\n%s\n%s\n%s\n"%(info2,seq2_to_write,"+",qual2))
				amp_same_key = "*"

				if amp1 != "NA" and amp2 != "NA":
					if amp1 != amp2:
						chimeric_reads_count += 1
					elif amp1_aln != "NA" and amp1_aln != amp1:
						aligned_other_loc_reads_count += 1
					elif amp2_aln != "NA" and amp2_aln != amp2:
						aligned_other_loc_reads_count += 1
					else:
						unidentified_reads_count += 1
				else:
					unidentified_reads_count += 1
		#create the dict based on the primer seqs
		#if there's a translocation, we don't expect the genome alignment to show that
		aln_count[(amp_same_key,amp1,amp2,amp1_aln,amp2_aln)] += 1

	unidentified_out1.close()
	unidentified_out2.close()
	for amp_name in amp_filehandles:
	   amp_filehandles[amp_name][0].close()
	   amp_filehandles[amp_name][1].close()

	identified_amplicon_file = output_root + ".splitReads.valid_amps.txt"
	with open(identified_amplicon_file,'w') as fout:
		for amp in amplicon_names:
			amp_aln_count = amplicon_count[amp]
			amplicon_information[amp]['aln_count'] = str(amp_aln_count)
			fout.write("%s\t%d\n"%(amp,amp_aln_count))
	

	logging.info("Cells not meeting the cell requirement: " + str(len(failing_barcode))) 
	logging.info("Total Barcode Count: " + str(len(tot_barcode_count)))
		   
	aligned_read_file = output_root + ".splitReads.aligned.txt"
	logging.info("Writing aligned reads to " + aligned_read_file)
	with open(aligned_read_file,'w') as fout:
		fout.write("Barcode\tAligned Count\n")
		for barcode in aln_barcode_count:
			fout.write("%s\t%d\n"%(barcode,aln_barcode_count[barcode]))
	logging.info("Finished writing aligned read file")
	
	unaligned_read_file = output_root + ".splitReads.unaligned.txt"
	logging.info('Writing unaligned reads to ' + unaligned_read_file)
	with open(unaligned_read_file,'w') as fout:
		fout.write("Barcode\tUnaligned Count\n")
		for barcode in unaln_barcode_count:
			fout.write("%s\t%d\n"%(barcode,unaln_barcode_count[barcode]))
	logging.info("Finished writing unaligned read file")

	identified_amplicon_aln_file = output_root + ".splitReads.amp_classification.txt"
	with open(identified_amplicon_aln_file,'w') as fout:
		fout.write("\t".join(['is_valid','amp1_from_seq','amp2_from_seq','amp1_from_align','amp2_from_align'])+"\n")
		for amp_key in sorted(aln_count.keys()):
			fout.write("%s\t%d\n"%("\t".join(amp_key),aln_count[amp_key]))

	amp_filenames.append((unidentified_out1_name,unidentified_out2_name))

	#zip output
	amp_commands = []
	for amp_filename_r1,amp_filename_r2 in amp_filenames:
		amp_commands.append('gzip -f ' + amp_filename_r1)
		amp_commands.append('gzip -f ' + amp_filename_r2)


	logging.info("gzipping output on "+ str(n_processes) + " threads..")

	#print(f"Line 2437:\n{amp_commands}")

	pool = mp.Pool(n_processes)
	pool.map_async(run_command, amp_commands).get(threading.TIMEOUT_MAX) 
	
	pool.close()
	pool.join()

	for amp in amplicon_names:
		r1_candidate = amplicon_information[amp].get('reads_r1_file')
		zipped_r1 = r1_candidate
		if zipped_r1 and not os.path.isfile(zipped_r1) and os.path.isfile(zipped_r1 + '.gz'):
			zipped_r1 = zipped_r1 + '.gz'
		if zipped_r1 and os.path.isfile(zipped_r1):
			#print(f"Writing\t{zipped_r1}")
			amplicon_information[amp]['reads_r1_file'] = zipped_r1
		else:
			#print(f"Line 2446\n{amp}\n{zipped_r1}\nExists?|{os.path.isfile(zipped_r1) if zipped_r1 else False}\n")
			amplicon_information[amp]['reads_r1_file'] = None

		r2_candidate = amplicon_information[amp].get('reads_r2_file')
		zipped_r2 = r2_candidate
		if zipped_r2 and not os.path.isfile(zipped_r2) and os.path.isfile(zipped_r2 + '.gz'):
			zipped_r2 = zipped_r2 + '.gz'
		if zipped_r2 and os.path.isfile(zipped_r2):
			amplicon_information[amp]['reads_r2_file'] = zipped_r2
		else:
			amplicon_information[amp]['reads_r2_file'] = None

	log_str = str(tot_reads_count) + " alignment pairs were processed (including multi-mapped and unaligned reads)\n" + \
			"  Of these, " + str(aln_reads_count) + " read pairs were aligned to the genome\n" + \
			"    Of these, " + str(id_reads_count) + " read pairs were aligned correctly and identified\n" + \
			"      " + str(partial_alignment_rescue_count) + " read pairs were rescued because both primer calls agreed and one alignment-side check supported the amplicon without contradiction\n" + \
			"    Of the unaligned reads, \n" + \
			"      " + str(chimeric_reads_count) + " read pairs were chimeric (contained primer sequences from different amplicons) \n" + \
			"      " + str(aligned_other_loc_reads_count) + " read pairs aligned to a different genomic location than the amplicon location \n" + \
			"      " + str(unidentified_reads_count) + " reads were otherwise unidentified\n"
	logging.info(log_str)
	log_file = output_root+".splitReads.log"
	with open(log_file,'w') as fout:
		fout.write(log_str)

	with open(info_file,'w') as fout:
		header_els = [
					'name',
					'amp_seqs',
					'guide_seq',
					'input_amp_seqs',
					'input_alternate_allele_seqs',
					'input_ref_allele_counts',
					'aln_chr',
					'aln_start',
					'aln_end',
					'aln_score',
					'secondary_aln_chr',
					'secondary_aln_start',
					'secondary_aln_end',
					'secondary_aln_score',
					'aln_count',
					'reads_r1_file',
					'reads_r2_file',
					]
		fout.write("\t".join(header_els)+"\n")
		for amplicon_name in amplicon_names:
			fout.write(
			"\t".join(
				str(amplicon_information[amplicon_name][x])
				if amplicon_information[amplicon_name][x] is not None
				else ""
				for x in header_els
				) + "\n"
			)

	if not keep_intermediate_files:
		logging.debug('Deleting intermediate amplicon files')
		delete_command = 'rm %s && rm %s'%(amplicon_fasta_file,aligned_amps_file)
		run_command(delete_command)

	return amplicon_names,amplicon_information,info_file

def run_crispresso_commands(amplicon_names,amplicon_information,output_root,crispresso_dir,suppress_sub_crispresso_plots,n_processes, alleles):
	"""
	Generate and execute CRISPResso2 commands for each amplicon.

	Depending on the `alleles` flag, this function runs CRISPResso2 on:
	  - Per-amplicon paired-end read FASTQs (standard mode), or
	  - Per-amplicon allele-only FASTQs (allele mode).

	Completed runs are detected via presence of a `.finished` file to
	support resumable execution.

	Parameters
	----------
	amplicon_names : list[str]
		List of amplicon identifiers.
	amplicon_information : dict
		Mapping amplicon_name -> metadata dict produced by
		`split_reads_by_amplicon`.
	output_root : str
		Base path for pipeline outputs.
	crispresso_dir : str
		Directory where CRISPResso outputs will be written.
	suppress_sub_crispresso_plots : bool
		If True, suppress CRISPResso2 report/plot generation.
	n_processes : int
		Number of worker processes to use for parallel execution.
	alleles : bool
		If True, run CRISPResso2 on allele-only FASTQs instead of
		paired-end reads.

	Returns
	-------
	dict
		Mapping amplicon_name -> result dictionary with keys:
			- 'name'
			- 'crispresso_command'
			- 'crispresso_run_folder'
			- 'finished_file'
			- 'log_file'
			- 'crispresso_result'
			- 'status' ('Completed', 'Failed', or 'Skipped')

	Notes
	-----
	- Commands are executed in parallel using multiprocessing.
	- Output metadata is written to:
		* `<output_root>.crispresso.info.txt`
		* `<output_root>.crispresso.filtered.info.txt` (allele mode)
	- Skips amplicons with zero aligned reads.
	"""
	# add prefix to info_file
	if alleles:
		info_file = output_root+".crispresso.filtered.info.txt"
		crispresso_dir = crispresso_dir + ".filtered"
		if not os.path.isdir(crispresso_dir):
			os.mkdir(crispresso_dir)
		
	else:
		info_file = output_root+".crispresso.info.txt"
	
	
	if os.path.isfile(info_file):
		with open(info_file,'r') as fin:
			head = fin.readline().strip()
			head_els = head.split("\t")
			crispresso_information = {}
			for line in fin:
				line_els = line.strip().split("\t")
				amp_info = dict(zip(head_els,line_els))
				crispresso_information[line_els[0]] = amp_info

		logging.info ("Finished running CRISPResso on targets")
		return crispresso_information

	crispresso_commands = []
	crispresso_information = {}

	not_run_count = 0
	finished_count = 0
	to_run_count = 0

	#print('Got to 2553')
	
	for amplicon_name in amplicon_names:
		crispresso_information[amplicon_name] = {}
		crispresso_information[amplicon_name]['name'] = amplicon_name
		if amplicon_information[amplicon_name]['aln_count'] == '0':
			#print(f"Got to 2559: {amplicon_name} within skipped block")
			crispresso_information[amplicon_name]['status'] = 'Skipped'
			crispresso_information[amplicon_name]['crispresso_command'] = 'NA'
			crispresso_information[amplicon_name]['crispresso_result'] = 'Skipped because had too few (%s) reads'%amplicon_information[amplicon_name]['aln_count']
			crispresso_information[amplicon_name]['finished_file'] = 'NA'
			crispresso_information[amplicon_name]['log_file'] = 'NA'
			crispresso_information[amplicon_name]['crispresso_run_folder'] = 'NA'
			not_run_count += 1
		else:

			#print(f"Running on {amplicon_name}\n{amplicon_information[amplicon_name]}\n")

			# These files are handled by updating crispresso_dir
			finished_file = os.path.join(crispresso_dir,amplicon_name+".finished")
			log_file = os.path.join(crispresso_dir,amplicon_name+".log")
			  
			#else:
			if not alleles:
				amp_filename_r1 = amplicon_information[amplicon_name].get('reads_r1_file')
				amp_filename_r2 = amplicon_information[amplicon_name].get('reads_r2_file')
			
			amplicon_seqs = amplicon_information[amplicon_name]['amp_seqs']
			guide = amplicon_information[amplicon_name]['guide_seq']
			guide_str = " -g " + guide + " "
			if guide.lower() == "na" or guide.lower() == "none":
				guide_str = ""
			suppress_sub_crispresso_plots_str = ""
			if suppress_sub_crispresso_plots:
				suppress_sub_crispresso_plots_str = " --suppress_report --suppress_plots"
			
			
			# Separate crispresso_cmd for alleles and non alleles
			if alleles:
				# set pass allele_seq file through crispresso
				allele_dir = os.path.join(output_root + ".seq_by_amplicon")
				#amp_filename = os.path.join(allele_dir, amplicon_name + "_filtered_allele.fq.gz")
				amp_filename = build_stage_filename(
					stage = STAGE_FILTER,
					tag = "alleles_qc_cells",
					amplicon = amplicon_name,
					ext = "fq.gz",
					output_root = allele_dir
				)
				#print(f"Got to line 2609\n{amp_filename}")

				if not amp_filename or not os.path.isfile(amp_filename):
					skip_reason = "Skipped because allele FASTQ input was missing after upstream filtering: %s" % amp_filename
					crispresso_information[amplicon_name]['status'] = 'Skipped'
					crispresso_information[amplicon_name]['crispresso_command'] = 'NA'
					crispresso_information[amplicon_name]['crispresso_result'] = skip_reason
					crispresso_information[amplicon_name]['finished_file'] = 'NA'
					crispresso_information[amplicon_name]['log_file'] = 'NA'
					crispresso_information[amplicon_name]['crispresso_run_folder'] = 'NA'
					logging.warning(skip_reason)
					not_run_count += 1
					continue
				
				if guide != "":
					crispresso_cmd = "CRISPResso -r1 " + amp_filename + " -a " + amplicon_seqs + " -g " + guide + suppress_sub_crispresso_plots_str + " -o "+crispresso_dir+" -n "+amplicon_name+" -w 2 --fastq_output --no_rerun --force_merge_pairs --exclude_bp_from_left 0 --exclude_bp_from_right 0 &> "+log_file +" && touch "+finished_file
				else:
					crispresso_cmd = "CRISPResso -r1 " + amp_filename + " -a " + amplicon_seqs + suppress_sub_crispresso_plots_str + " -o "+crispresso_dir+" -n "+amplicon_name+" -w 2 --fastq_output --no_rerun --force_merge_pairs --exclude_bp_from_left 0 --exclude_bp_from_right 0 &> "+log_file +" && touch "+finished_file    
				
			else:
				missing_inputs = []
				if not amp_filename_r1 or not os.path.isfile(amp_filename_r1):
					missing_inputs.append(amp_filename_r1 or 'missing_r1')
				if not amp_filename_r2 or not os.path.isfile(amp_filename_r2):
					missing_inputs.append(amp_filename_r2 or 'missing_r2')
				if missing_inputs:
					skip_reason = "Skipped because read FASTQ input(s) were missing after upstream filtering: %s" % ", ".join(missing_inputs)
					crispresso_information[amplicon_name]['status'] = 'Skipped'
					crispresso_information[amplicon_name]['crispresso_command'] = 'NA'
					crispresso_information[amplicon_name]['crispresso_result'] = skip_reason
					crispresso_information[amplicon_name]['finished_file'] = 'NA'
					crispresso_information[amplicon_name]['log_file'] = 'NA'
					crispresso_information[amplicon_name]['crispresso_run_folder'] = 'NA'
					logging.warning(skip_reason)
					not_run_count += 1
					continue
				if guide != "":
					crispresso_cmd = "CRISPResso -r1 " + amp_filename_r1 + " -r2 " + amp_filename_r2 + " -a "+amplicon_seqs+" -g "+guide+ suppress_sub_crispresso_plots_str + " -o "+crispresso_dir+" -n "+amplicon_name+" -w 2 --fastq_output --no_rerun --force_merge_pairs --exclude_bp_from_left 0 --exclude_bp_from_right 0 &> "+log_file +" && touch "+finished_file
				else:
					crispresso_cmd = "CRISPResso -r1 " + amp_filename_r1 + " -r2 " + amp_filename_r2 + " -a "+amplicon_seqs+ suppress_sub_crispresso_plots_str + " -o "+crispresso_dir+" -n "+amplicon_name+" -w 2 --fastq_output --no_rerun --force_merge_pairs --exclude_bp_from_left 0 --exclude_bp_from_right 0 &> "+log_file +" && touch "+finished_file    
			
			crispresso_information[amplicon_name]['crispresso_command'] = crispresso_cmd
			crispresso_information[amplicon_name]['finished_file'] = finished_file
			crispresso_information[amplicon_name]['log_file'] = log_file
			crispresso_information[amplicon_name]['crispresso_run_folder'] = os.path.join(crispresso_dir,'CRISPResso_on_'+amplicon_name)
			
			if os.path.isfile(finished_file):
				finished_count += 1
				continue
			else:
				crispresso_commands.append(crispresso_cmd)


	logging.info('Skipped CRISPResso analysis for ' + str(not_run_count) + ' amplicons, finished analysis for ' + str(finished_count) + ' amplicons')

	logging.info('Got ' + str(len(crispresso_commands)) + ' CRISPResso commands')

	if len(crispresso_commands) > 0:
		# start processes
		logging.info("Running on "+ str(n_processes) + " processes..")
		pool = mp.Pool(n_processes)
		result = pool.map_async(run_command, crispresso_commands).get(threading.TIMEOUT_MAX)
		pool.close()
		pool.join()

	for amplicon_name in amplicon_names:
		if 'status' in crispresso_information[amplicon_name] and crispresso_information[amplicon_name]['status'] == 'Skipped':
			pass
		else:
			finished_file = crispresso_information[amplicon_name]['finished_file']
			if os.path.isfile(finished_file):
				finished_file = crispresso_information[amplicon_name]['finished_file']
				crispresso_information[amplicon_name]['status'] = 'Completed'
				crispresso_information[amplicon_name]['crispresso_result'] = 'Completed'
			else:
				crispresso_information[amplicon_name]['status'] = 'Failed'
				log_file = crispresso_information[amplicon_name]['log_file']
				error_message = 'Failed, see ' + log_file
				if os.path.isfile(log_file):
					with open(log_file,'r') as lf:
						for line in lf:
							if 'ERROR:' in line:
								error_message = line.strip()
				crispresso_information[amplicon_name]['crispresso_result'] = error_message


	with open(info_file,'w') as fout:
		header_els = [
					'name',
					'crispresso_command',
					'crispresso_run_folder',
					'finished_file',
					'log_file',
					'crispresso_result',
					'status',
					]
		fout.write("\t".join(header_els)+"\n")
		for amplicon_name in amplicon_names:
			fout.write("\t".join([crispresso_information[amplicon_name][x] for x in header_els])+"\n")
	return crispresso_information

def run_command(cmd):
	"""
	Execute a shell command using subprocess.

	Parameters
	----------
	cmd : str
		Shell command to execute.

	Returns
	-------
	int or None
		Return code from subprocess call, or None if execution failed.

	Notes
	-----
	- Executes with `shell=True`.
	- Logs the command at DEBUG level.
	- Errors are logged but not raised.
	- Does not capture stdout.
	"""
	try:
		logging.debug('running: ' + cmd)
		return_value = sb.call(cmd,shell=True)
	except Exception as e:
		logging.error("error: %s on %s" % (e, cmd))

def get_command_output(command):
	"""
	Execute a shell command and return an iterator over stdout lines.

	Parameters
	----------
	command : str
		Shell command to execute.

	Returns
	-------
	iterator[str]
		Iterator yielding lines from command stdout.

	Notes
	-----
	- Uses subprocess.Popen with shell=True.
	- stderr is redirected to stdout.
	- Caller is responsible for consuming the iterator.
	""" 
	p = sb.Popen(command,
			stdout=sb.PIPE,
			stderr=sb.STDOUT,shell=True,
			universal_newlines=True,
			bufsize=-1)#bufsize system default
	return iter(p.stdout.readline, b'')

def parse_one_crispresso_output(this_args):
	"""
		Parse a single CRISPResso2 output folder into per-cell allele summaries.

	Reads the CRISPResso_output.fastq.gz file, extracts:
		- Cell barcode
		- Reference alignment
		- Indel/substitution status
		- Allele counts per reference

	Performs multinomial-based allele assignment to determine
	the most likely allele configuration per cell.

	Parameters
	----------
	this_args : dict
		Dictionary containing:
			- amplicon_name : str
			- amplicon_info_file : str
			- crispresso_run_folder : str
			- input_ref_allele_counts : str
			- min_num_reads_per_cell : int
			- min_allele_pct_cutoff : float
			- min_allele_count_cutoff : int
			- ignore_substitutions : bool
			- output_root : str
			- min_reads_per_amplicon_per_cell : int

	Returns
	-------
	None
		Results are written to disk:
			- <run_folder>.summ
			- <run_folder>.summ.finished
			- allele FASTQ files (if enabled)

	Notes
	-----
	- Intended for multiprocessing execution.
	- Assumes CRISPResso output structure is unchanged.
	- Multinomial modeling is used to distinguish signal from noise.
	"""
	amplicon_name = this_args['amplicon_name']
	amplicon_info_file = this_args['amplicon_info_file']
	crispresso_run_folder = this_args['crispresso_run_folder']
	input_ref_allele_counts = this_args['input_ref_allele_counts']
	min_num_reads_per_cell = this_args['min_num_reads_per_cell']
	min_allele_pct_cutoff = this_args['min_allele_pct_cutoff']
	min_allele_count_cutoff = this_args['min_allele_count_cutoff']
	ignore_substitutions = this_args['ignore_substitutions']
	amplicon_dir = os.path.join(this_args['output_root'] + ".seq_by_amplicon")
	min_reads_per_amplicon_per_cell = this_args['min_reads_per_amplicon_per_cell']
	

	folder_finished_file = crispresso_run_folder + ".summ.finished"
	crispresso_output_fastq = os.path.join(crispresso_run_folder, 'CRISPResso_output.fastq.gz')
	amp_arm_check_len = 30
	
	wildtype_allele = get_wildtype_allele(crispresso_run_folder)
	


	with open(amplicon_info_file,'r') as fin:
		head = fin.readline().strip()
		head_els = head.split("\t")
		amplicon_names = []
		amplicon_information = {}
		for line in fin:
			line_els = line.strip().split("\t")
			amp_info = dict(zip(head_els,line_els))
			amplicon_information[line_els[0]] = amp_info
			amplicon_names.append(line_els[0])

	this_amplicon_info = amplicon_information[amplicon_name]
	ok_left_sides = [x[0:amp_arm_check_len] for x in this_amplicon_info['amp_seqs'].split(",")]
	ok_right_sides = [x[-1*amp_arm_check_len:] for x in this_amplicon_info['amp_seqs'].split(",")]
	tot_count = 0
	crispresso2_aligned_count = 0
	data = {}
	alleles = {}
	allele_sequence_dict = {}
	seen_refs = []
	cell_read_counts = defaultdict(int)
	num_crispresso_references = 0
	num_references = len(input_ref_allele_counts.split(","))
	logging.debug('Parsing CRISPResso output for ' + amplicon_name)
	fastq_input_handle =  gzip.open(crispresso_output_fastq,'rt')
	next_fastq_id = fastq_input_handle.readline()
	while(next_fastq_id):
		#read through fastq in sets of 4
		fastq_id = next_fastq_id.split(" ")[0] #fastp adds ' merged_199_234' so trim that off
		fastq_seq = fastq_input_handle.readline().strip()
		fastq_plus = fastq_input_handle.readline().strip()
		fastq_qual = fastq_input_handle.readline()
		next_fastq_id = fastq_input_handle.readline()

		tot_count += 1
		if "ALN=NA " in fastq_plus: # Read did not align
			continue
		crispresso2_aligned_count += 1
		id_els = fastq_id.strip().split(":")
		cell = id_els[-1]

		known_amp_left = False
		known_amp_right = False
		if fastq_seq[0:amp_arm_check_len] in ok_left_sides:
			known_amp_left = True
		if fastq_seq[-1*amp_arm_check_len:] in ok_right_sides:
			known_amp_right = True

		if not ok_left_sides or not ok_right_sides:
			#print('mismatch: ' + fastq_seq[0:amp_arm_check_len] + ' with ' + str(ok_left_sides))
			#print('mismatch: ' + fastq_seq[-1*amp_arm_check_len:] + ' with ' + str(ok_right_sides))
			continue

		aln_ref = ""
		#match = re.search(" ALN=(\S+) ", fastq_plus)
		match = re.search(r" ALN=(\S+) ", fastq_plus)
		if match:
			aln_ref = match.group(1)
		#discard reads that align ambiguously
		if '&' in aln_ref:
			continue

		if cell not in data:
			data[cell] = {'mod':0,'unmod':0}
			alleles[cell] = {}
			allele_sequence_dict[cell] = {}

		if aln_ref not in data[cell]:
			data[cell][aln_ref] = {'mod':0,'unmod':0}
			if aln_ref not in seen_refs:
				seen_refs.append(aln_ref)


		allele = "NA"

		# Formation of allele_key should only consider the quant window in the gRNA
		if ignore_substitutions:
			match = re.search("(DEL=.* INS=.*) SUB=.* ALN_REF", fastq_plus)
			unmod_allele_str = "DEL= INS="
		else:
			match = re.search("(DEL=.* INS=.* SUB=.*) ALN_REF", fastq_plus)
			unmod_allele_str = "DEL= INS= SUB="

		if match:
			allele = match.group(1)
			
			# Check for gRNA input
			# if gRNA, where in the amplicon sequence?
			# Check for allele key values outside of gRNA
			# if outside of gRNA, convert to WT read
			# We don't expect CRISPR edits outside of gRNA region
			
			if allele == unmod_allele_str:
				data[cell]['unmod'] += 1
				data[cell][aln_ref]['unmod'] += 1
			else:
				data[cell]['mod'] += 1
			data[cell][aln_ref]['mod'] += 1
		allele_key = aln_ref + ":" + allele
		
		if allele_key not in alleles[cell]:
			alleles[cell][allele_key] = 0
			# new layer with sequence + count
			allele_sequence_dict[cell][allele_key] = {}
		
		if fastq_seq not in allele_sequence_dict[cell][allele_key]:
			allele_sequence_dict[cell][allele_key][fastq_seq] = 0
		
		alleles[cell][allele_key] += 1
		allele_sequence_dict[cell][allele_key][fastq_seq] += 1
		cell_read_counts[cell] += 1

	# Checking for proper allele_sequence_dict formation


	#if somehow CRISPResso is reporting more than the input number of references, throw this warning and reset the number of references to match that seen in CRISPResso
	if len(seen_refs) > num_references:
		logging.warning('WARNING - saw ' + str(seen_refs) + ' refs for cell ' + cell + ' and ' + amplicon_name + ' (Expecting only ' + str(num_references) + ')')
		num_references = len(seen_refs)

	#set up allele order
	input_ref_names = ['Reference']
	for i in range(1,num_references):
		input_ref_names[i] = 'Amplicon'+str(i)

	if 'NA' in input_ref_allele_counts:
		print('WARNING, NA input ref count!')
		input_ref_allele_counts = "2"
	input_ref_allele_counts = [int(x) for x in input_ref_allele_counts.split(",")]
	#if CRISPResso reports more alternate alleles, just assign them a presence of 1
	while len(input_ref_allele_counts) < num_references:
		input_ref_allele_counts.append(1)

	tot_allele_count = sum(input_ref_allele_counts)

	#create prob array for multinomial tests
	noise_prob = 0.01
	input_ref_allele_probs = []
	for this_allele_count in input_ref_allele_counts:
		alleles_prob = (1-noise_prob)/float(this_allele_count)
		prob_array = [alleles_prob]*this_allele_count
		prob_array.append(noise_prob)
		input_ref_allele_probs.append(prob_array)


	#os.path.join(amplicon_dir, amplicon_name + "_unfiltered_allele.fq")
	unfiltered_allele_file = build_stage_filename(
		stage = STAGE_SPLIT,
		tag = "alleles_all_cells",
		amplicon = amplicon_name,
		ext = "fq",
		output_root = amplicon_dir
	)


	#print(f"Got here 2868\n{unfiltered_allele_file=}")

	with open(crispresso_run_folder+".summ",'w') as fout, open(crispresso_run_folder+".summarize_indels.out",'w') as fsumm, open(crispresso_run_folder+".summarize_alleles.out",'w') as asumm, open(unfiltered_allele_file, "w") as aseq:
		fout.write("\t".join([str(x) for x in ['cell','all_cell_read_count','all_cell_mut_pct','all_cell_allele_string','final_cell_read_count','final_cell_mut_allele_pct','final_cell_allele_string','final_num_refs_covered','final_cell_allele_mod_string','final_cell_allele_mod_types_string','final_cell_allele_readcount_string','final_ref_read_count_string','final_ref_mut_allele_fracs_string']])+"\n")


		asumm.write("cell\tread_count\tmod_pct\t"+"\t".join(["allele_"+str(x) for x in range(tot_allele_count)]) + "\n")
		fsumm.write("cell\tread_count\tmod_pct\t"+"\t".join(["allele_"+str(x) for x in range(tot_allele_count)]) + "\n")
		for cell in sorted(data.keys()):
			if cell.strip() == "":
				continue
			#logging.debug('cell is ' + cell + ' with ' + str(cell_read_counts[cell]) + ' reads')
			mod_count = data[cell]['mod']
			unmod_count = data[cell]['unmod']
			all_cell_read_count = mod_count + unmod_count
			all_cell_mut_pct = round(100*mod_count/float(all_cell_read_count),2)

			final_alleles = []
			all_cell_alleles = sorted(alleles[cell].items(), key=lambda x:x[1],reverse=True)
			for idx, this_allele_count in enumerate(input_ref_allele_counts):
				this_allele_name = input_ref_names[idx] #reference name to look in CRISPResso output for
				this_prob_array = input_ref_allele_probs[idx] #probability array to use for multinomial
				this_cell_alleles = [x for x in alleles[cell].items() if x[0].split(":")[0] == this_allele_name] #alleles that match to this specific reference

				this_cell_alleles = sorted(this_cell_alleles, key=lambda x:x[1],reverse=True)

				this_cell_alleles_sum = sum([x[1] for x in this_cell_alleles])

				# i is the number of alleles from the final_alleles to test as real. The rest are noise.
				best_prob = None
				best_alleles = None
				for i in range(1,this_allele_count+1):
					# python indexing works in our favor here and will return nothing for accesses past the array length. e.g. d = [1,2]; d[0:5] = [1,2] and d[5:] = []
					alleles_real = this_cell_alleles[0:i]
					alleles_noise = this_cell_alleles[i:]

					if len(alleles_real) == 0:
						alleles_real = [('NA',0)]
					#distribute the chosen alleles over the num_max_alleles
					while len(alleles_real) < this_allele_count:
#                        #print('beginning ' + str(alleles_real))
						allele_to_halve = alleles_real.pop()
						half_count = int(allele_to_halve[1]/2)
						half_allele = (allele_to_halve[0],half_count)
						half_allele_2 = (allele_to_halve[0],allele_to_halve[1]-half_count)
						alleles_real = sorted(alleles_real+[half_allele,half_allele_2], key=lambda x:x[1],reverse=True)

					alleles_real_counts = [x[1] for x in alleles_real]
					alleles_noise_count = sum([x[1] for x in alleles_noise])

					#this is the array of counts for the multinomial
					prob_counts = alleles_real_counts + [alleles_noise_count]

					this_prob = multinomial.pmf(prob_counts,this_cell_alleles_sum,this_prob_array)
#                    print('pmf of ' + str(prob_counts) + ' cellall: ' + str(this_cell_alleles_sum) + ' prob array : ' + str(prob_array))
#                    print('this prob: ' + str(this_prob))
					if best_prob is None or this_prob > best_prob:
						best_prob = this_prob
						best_alleles = alleles_real



				for best_allele in best_alleles:
					final_alleles.append(best_allele)
			


			#done performing reference-specific assignment

			all_cell_allele_arr = []
			for idx,(allele,count) in enumerate(all_cell_alleles):
				#add all alleles to string (unfiltered)
				all_cell_allele_arr.append(allele+":"+str(count))
			all_cell_allele_string = ",".join(all_cell_allele_arr)

			final_ref_read_counts = defaultdict(int) #for each ref, how many final reads were there
			final_ref_mod_allele_counts = defaultdict(int)
			final_ref_unmod_allele_counts = defaultdict(int)
			final_cell_allele_arr = [] # list of alleles
			final_cell_allele_mod_arr = [] # U/M for modified
			final_cell_allele_mod_types_arr = [] #DIS for deletion, insertion, sub
			final_cell_allele_readcount_arr = [] #number of reads per allele

			final_mod_allele_count = 0 #how many mod alleles
			final_unmod_allele_count = 0 #how many unmod alleles
			final_cell_read_count = 0 # how many total final reads for this cell
			for idx,(allele,count) in enumerate(final_alleles):
				# add final alleles to string
				final_cell_read_count += count
				allele_ref,allele_status = allele.split(":")
				final_cell_allele_arr.append(allele)

				final_ref_read_counts[allele_ref] += count
				final_cell_allele_readcount_arr.append(str(count))
				this_allele_mut_type_str = ""
				if ignore_substitutions:
					if allele_status == 'DEL= INS=':
						final_unmod_allele_count += 1
						final_cell_allele_mod_arr.append("U")
						final_ref_unmod_allele_counts[allele_ref] += 1
						this_allele_mut_type_str = "U"
					else:
						final_mod_allele_count += 1
						final_cell_allele_mod_arr.append("M")
						final_ref_mod_allele_counts[allele_ref] += 1
						(del_str,ins_str) = [x.split("=")[1] for x in allele_status.split(" ")]
						if del_str != '':
							this_allele_mut_type_str += "D"
						if ins_str != '':
							this_allele_mut_type_str += "I"
				else: #include substitutions
					if allele_status == 'DEL= INS= SUB=':
						final_unmod_allele_count += 1
						final_cell_allele_mod_arr.append("U")
						final_ref_unmod_allele_counts[allele_ref] += 1
						this_allele_mut_type_str = "U"
					else:
						final_mod_allele_count += 1
						final_cell_allele_mod_arr.append("M")
						final_ref_mod_allele_counts[allele_ref] += 1
						(del_str,ins_str,sub_str) = [x.split("=")[1] for x in allele_status.split(" ")]
						if del_str != '':
							this_allele_mut_type_str += "D"
						if ins_str != '':
							this_allele_mut_type_str += "I"
						if sub_str != '':
							this_allele_mut_type_str += "S"
				final_cell_allele_mod_types_arr.append(this_allele_mut_type_str)


			#if final_cell_read_count >= read_count_per_amplicon_cutoff:                
			#    write_max_alleles(allele_sequence_dict,
			#                    cell,
			#                    final_cell_allele_arr,
			#                    amplicon_name,
			#                    amplicon_dir,
			#                    aseq,
			#                    wildtype_allele) 

			#print(f"At line 2999: write_max_alleles: {aseq=}")
			#asdf()
			write_max_alleles(allele_sequence_dict, cell, final_cell_allele_arr, amplicon_name, amplicon_dir, aseq, wildtype_allele)

			final_cell_allele_string = ",".join(final_cell_allele_arr)
			final_cell_allele_mod_string = ",".join(final_cell_allele_mod_arr)
			final_cell_allele_mod_types_string = ",".join(final_cell_allele_mod_types_arr)
			final_cell_allele_readcount_string = ",".join(final_cell_allele_readcount_arr)
			final_cell_mut_allele_pct = round(100*final_mod_allele_count/float(final_mod_allele_count + final_unmod_allele_count),2)

			#now compute for each reference
			final_num_refs_covered = len(final_ref_read_counts.keys())
			final_ref_mut_allele_fracs = ["NA"]*num_references
			final_ref_read_count = [0]*num_references
			for idx, this_ref_count in enumerate(input_ref_allele_counts):
				this_ref_name = input_ref_names[idx] #reference name to look in CRISPResso output for
				final_ref_read_count[idx] = final_ref_read_counts[this_ref_name]

				this_ref_mod = final_ref_mod_allele_counts[this_ref_name]
				this_ref_unmod = final_ref_unmod_allele_counts[this_ref_name]
				this_ref_tot = this_ref_mod + this_ref_unmod
				if this_ref_tot > 0:
					final_ref_mut_allele_fracs[idx] = round(100*this_ref_mod/float(this_ref_tot),2)
			final_ref_mut_allele_fracs_string = ",".join([str(x) for x in final_ref_mut_allele_fracs])
			final_ref_read_count_string = ",".join([str(x) for x in final_ref_read_count])
	
			# For each barcode, parse the allele dict and return the most sequence for the allele 
			fout.write("\t".join([str(x) for x in [cell,all_cell_read_count,all_cell_mut_pct,all_cell_allele_string,final_cell_read_count,final_cell_mut_allele_pct,final_cell_allele_string,final_num_refs_covered,final_cell_allele_mod_string,final_cell_allele_mod_types_string,final_cell_allele_readcount_string,final_ref_read_count_string,final_ref_mut_allele_fracs_string]])+"\n")
			if all_cell_read_count >= min_num_reads_per_cell:
				fsumm.write("\t".join([str(x) for x in [cell,final_cell_read_count,final_cell_mut_allele_pct]])+"\n")

				asumm.write("\t".join([str(x) for x in [cell,final_cell_read_count,final_cell_mut_allele_pct]+final_alleles])+"\n")

	with open (folder_finished_file,'w') as fout:
		fout.write("Total reads\t" + str(tot_count)+"\n")
		fout.write("CRISPResso2 aligned reads\t" + str(crispresso2_aligned_count)+"\n")
		fout.write(str(datetime.now()))

def write_max_alleles(allele_dict, barcode, allele_key, amplicon_name, amplicon_folder, allele_file, wildtype_allele):
	"""
	 Write consensus allele sequences for a barcode to a FASTQ file.

	For each allele key assigned to a barcode, selects the most common
	sequence. In case of ties:
		- If counts are 1 and wildtype allele is available, use wildtype.
		- Otherwise select the highest-frequency sequence.

	Parameters
	----------
	allele_dict : dict
		Nested mapping:
			barcode -> allele_key -> {sequence: count}
	barcode : str
		Cell barcode identifier.
	allele_key : str or list[str]
		Allele key(s) selected for the barcode.
	amplicon_name : str
		Amplicon identifier.
	amplicon_folder : str
		Path to amplicon output directory (unused but retained for interface consistency).
	allele_file : file-like object
		Open writable file handle for FASTQ output.
	wildtype_allele : str or None
		Wildtype allele sequence used for tie-breaking.

	Returns
	-------
	None

	Notes
	-----
	- Writes FASTQ-formatted records.
	- Does not return values; writes directly to `allele_file`.
	- Assumes allele_dict structure created by parse_one_crispresso_output.
	"""
	# Return the most common sequence for the alleles for each barcode
	allele_list = []
	barcode_allele_key = []
	adjusted_list = []


	if isinstance(allele_key, str):
		allele_key = [allele_key]


	if len(allele_key) == 1:
		allele_key.append(allele_key[0])

	if len(allele_key) == 0:
		return

	for allele in allele_key:
		adjusted = False
		allele_df = pd.DataFrame(list(allele_dict[barcode][allele].items()), columns = ['FASTA', 'Count'])
		max_value = allele_df['Count'].max()
		max_loc = allele_df['Count'] == max_value
		max_indices = allele_df.index[max_loc]
		if sum(max_loc) > 1: # If there are ties
			if max_value == 1: # If the tied values are 1, return the wildtype allele
				if wildtype_allele is None: # Check if a wildtype allele was found
					allele_sequence = allele_df.loc[max_indices[0], 'FASTA']
				else:
					allele_sequence = wildtype_allele
					adjusted = True # add a flag to report when we adjust a cell to the WT allele
					
			else: # If tied values are greater than 1, return the first one the path to the crispresso output for the amplicon
				allele_sequence = allele_df.loc[allele_df['Count'].idxmax(), 'FASTA']
		else: # There are no ties, choose highest frequency read
			allele_sequence = allele_df.loc[max_indices[0], 'FASTA']
		
		adjusted_list.append(adjusted)
		allele_list.append(allele_sequence)
		barcode_allele_key.append(allele)

	headers = []
	sequences = []
	qualities = []
	
	
	adjusted_list = ["Adjusted:" if x else "" for x in adjusted_list]
	
	for index, allele_key in enumerate(barcode_allele_key):
		headers.append(f"@{amplicon_name}:{adjusted_list[index]}{allele_key.replace(' ', '_')}:{barcode}:{index+1}")
		sequences.append(allele_list[index])
		qualities.append("a" * len(allele_list[index]))
	
	for h,s,q in zip(headers, sequences, qualities):
		allele_file.write(f"{h}\n{s}\n+\n{q}\n")
			
	return

def get_wildtype_allele(crispresso_run_folder):
	"""
	Parse the crispresso allele table and grab the 'wildtype' allele

	params:
		crispresso_run_folder: the path to the output directory of a crispresso run
		
	returns:
		max_allele: the sequence of the selected 'wildtype' allele for an amplicon
		
	"""
	crispresso2_info = CRISPRessoShared.load_crispresso_info(crispresso_run_folder)
	z = zipfile.ZipFile(os.path.join(crispresso_run_folder, crispresso2_info['running_info']['allele_frequency_table_zip_filename']))
	zf = z.open(crispresso2_info['running_info']['allele_frequency_table_filename'])
	df_alleles = pd.read_csv(zf, sep="\t")
	
	# Get the most common wild type allele
	df_alleles = df_alleles[(df_alleles['Read_Status'] == "UNMODIFIED") &
							(df_alleles['n_deleted'] == 0) &
							(df_alleles['n_inserted'] == 0) &
							(df_alleles['n_mutated'] == 0)]
	

	if len(df_alleles) == 0:
		return None
	
	max_allele = df_alleles.loc[df_alleles['#Reads'].idxmax()]
	
	max_allele = max_allele['Aligned_Sequence']
   
   
	return max_allele


def parse_crispresso_outputs(amplicon_names,amplicon_information,amplicon_info_file,crispresso_information,
							output_root, min_total_reads_per_barcode, min_reads_per_amplicon_per_cell, n_processes,num_max_alleles=2,num_references=1,
							min_num_reads_per_cell=5,min_allele_pct_cutoff=.1,min_allele_count_cutoff=2,
							ignore_substitutions=False,
							write_alleles=False):
	"""
	Generate and execute CRISPResso2 commands for each amplicon.

	Depending on the `alleles` flag, this function runs CRISPResso2 on:
	  - Per-amplicon paired-end read FASTQs (standard mode), or
	  - Per-amplicon allele-only FASTQs (allele mode).

	Completed runs are detected via presence of a `.finished` file to
	support resumable execution.

	Parameters
	----------
	amplicon_names : list[str]
		List of amplicon identifiers.
	amplicon_information : dict
		Mapping amplicon_name -> metadata dict produced by
		`split_reads_by_amplicon`.
	output_root : str
		Base path for pipeline outputs.
	crispresso_dir : str
		Directory where CRISPResso outputs will be written.
	suppress_sub_crispresso_plots : bool
		If True, suppress CRISPResso2 report/plot generation.
	n_processes : int
		Number of worker processes to use for parallel execution.
	alleles : bool
		If True, run CRISPResso2 on allele-only FASTQs instead of
		paired-end reads.

	Returns
	-------
	dict
		Mapping amplicon_name -> result dictionary with keys:
			- 'name'
			- 'crispresso_command'
			- 'crispresso_run_folder'
			- 'finished_file'
			- 'log_file'
			- 'crispresso_result'
			- 'status' ('Completed', 'Failed', or 'Skipped')

	Notes
	-----
	- Commands are executed in parallel using multiprocessing.
	- Output metadata is written to:
		* `<output_root>.crispresso.info.txt`
		* `<output_root>.crispresso.filtered.info.txt` (allele mode)
	- Skips amplicons with zero aligned reads.
	"""
	parse_output_args = []
	for name in amplicon_names:
		if crispresso_information[name]['status'] == 'Completed':
			crispresso_run_folder = crispresso_information[name]['crispresso_run_folder']
			crispresso_out = os.path.join(crispresso_run_folder,'CRISPResso_output.fastq.gz')
			input_ref_allele_counts = amplicon_information[name]['input_ref_allele_counts']
			folder_finished_file = crispresso_run_folder + ".summ.finished"

			if not os.path.isfile(folder_finished_file):
				this_args = {'amplicon_name':name,
							 'amplicon_info_file':amplicon_info_file,
							 'crispresso_run_folder':crispresso_run_folder,
							 'input_ref_allele_counts':input_ref_allele_counts,
							 'min_num_reads_per_cell':min_num_reads_per_cell,
							 'min_allele_pct_cutoff':min_allele_pct_cutoff,
							 'min_allele_count_cutoff':min_allele_count_cutoff,
							 'ignore_substitutions':ignore_substitutions,
							 'output_root': output_root,
							 'write_alleles': write_alleles,
							 'min_reads_per_amplicon_per_cell': min_reads_per_amplicon_per_cell
							 }
				parse_output_args.append(this_args)

	if len(parse_output_args) > 0:
		logging.info('Parsing ' + str(len(parse_output_args)) + ' CRISPResso folders on ' + str(n_processes) + ' threads..')
		# start processes
		if n_processes > 1 and len(parse_output_args) > 1:
			pool = mp.Pool(n_processes)
			pool.map_async(parse_one_crispresso_output, parse_output_args).get(threading.TIMEOUT_MAX)
			pool.close()
			pool.join()
		else:
			for this_args in parse_output_args:
				parse_one_crispresso_output(this_args)
	else:
		logging.info('Finished parsing CRISPResso folders')

	logging.info('Aggregating ' + str(len(amplicon_names)) + ' target summaries')
	data = {}
	for amplicon_name in amplicon_names:
		if crispresso_information[amplicon_name]['status'] != 'Completed':
			continue
		crispresso_run_folder = crispresso_information[amplicon_name]['crispresso_run_folder']
		summ_file = crispresso_run_folder + ".summ"
		if not os.path.isfile(summ_file):
			continue
		with open (summ_file,'r') as fin:
			head = fin.readline()
			for line in fin:
				line_els = line.strip().split("\t")
				if len(line_els) < 3:
					raise Exception('Unexpected line format: ' + line + ' in ' + summ_file)
				cell = line_els[0]
				all_cell_read_count = line_els[1]
				all_cell_mut_pct = line_els[2]
				final_cell_read_count = line_els[4]
				final_cell_mut_pct = line_els[5]

				if cell not in data:
					data[cell] = {}
				data[cell][amplicon_name]=(
						"\t"+all_cell_read_count+"\t"+all_cell_mut_pct,
						"\t"+final_cell_read_count+"\t"+final_cell_mut_pct)


	cells = sorted(data.keys())

	with open(output_root+".editingSummaryPseudobulk.txt",'w') as fout:
		header = "cell"
		for name in amplicon_names:
			header += "\ttotCount.%s\tmodPct.%s"%(name,name)
		fout.write(header+"\n")

		for cell in cells:
			line = cell
			for name in amplicon_names:
				val = "\t0\tNA"
				if name in data[cell]:
					val = data[cell][name][0]
				line += val
			fout.write(line+"\n")

	with open(output_root+".editingSummary.txt",'w') as fout:
		header = "cell"
		for name in amplicon_names:
			header += "\ttotCount.%s\tmodPct.%s"%(name,name)
		fout.write(header+"\n")

		for cell in cells:
			line = cell
			for name in amplicon_names:
				val = "\t0\tNA"
				if name in data[cell]:
					val = data[cell][name][1]
				line += val
			fout.write(line+"\n")

	logging.info("Finished reading and compiling summaries for %d cells"%len(cells))

	# Creating a formatted data frame object to use for plot generation
	df_colnames = []
	# Create colnames
	for name in amplicon_names:
		df_colnames.append("totCount.%s"%name)
		df_colnames.append("modPct.%s"%name)    
	# Create indices using cells
	df_index = [cell for cell in cells]

	# Empty prepped DataFrame
	summary_df = pd.DataFrame(columns = df_colnames, index = df_index)

	# Parse data dictionary and add to DataFrame
	for cell in cells:
		vals = []
		for amplicon in amplicon_names:
			val = [0, "NA"]
			if amplicon in data[cell]:
				val = data[cell][amplicon][0].strip().split('\t')
			val = [x if isinstance(x, int) else int(x) if x.isdigit() else x if x == "NA" else float(x) for x in val]
			vals += val
		summary_df.loc[cell] = vals
	
	
	# Filter to totCols
	totCols = summary_df.filter(like = "totCount")
	
	# Check if Amplicon Score File Exists
	amp_score_file = output_root + ".amplicon_score.txt"
	if os.path.isfile(amp_score_file): # if amp_score_file exists, read it in
		amp_score = pd.read_csv(amp_score_file, sep = "\t", index_col = 0)
	else: # Generate Amplicon Score and write it to a file
		amplicon_score_time = time.time()
		amp_score = generate_amplicon_score(totCols, min_reads_per_amplicon_per_cell = min_reads_per_amplicon_per_cell, min_total_reads_per_barcode = min_total_reads_per_barcode)
		amp_score.to_csv(amp_score_file, sep = "\t")
		end_amplicon_score_time = time.time() - amplicon_score_time
		logging.info("Generated amplicon score in %.2f seconds"%(end_amplicon_score_time))
	
	# writing out a filtered editingSummary file
	# this removes cells filtered out within generate_amplicon_score
	with open(output_root+".filteredEditingSummaryPseudobulk.txt",'w') as fout:
		header = "cell"
		for name in amplicon_names:
			header += "\ttotCount.%s\tmodPct.%s"%(name,name)
		fout.write(header+"\n")

		for cell in cells:
			if cell not in amp_score.index:
				continue
			line = cell
			for name in amplicon_names:
				val = "\t0\tNA"
				if name in data[cell]:
					val = data[cell][name][0]
				line += val
			fout.write(line+"\n")

	with open(output_root+".filteredEditingSummary.txt",'w') as fout:
		header = "cell"
		for name in amplicon_names:
			header += "\ttotCount.%s\tmodPct.%s"%(name,name)
		fout.write(header+"\n")

		for cell in cells:
			if cell not in amp_score.index:
				continue
			line = cell
			for name in amplicon_names:
				val = "\t0\tNA"
				if name in data[cell]:
					val = data[cell][name][1]
				line += val
			fout.write(line+"\n")

	logging.info("Finished reading and compiling summaries for %d filtered cells"%len(amp_score))


	summary_df = add_color_information(summary_df, amp_score)
	 
	return summary_df


def stratify_data(input_data):
	"""
	Assign quality category codes to barcodes.

	Barcodes are classified into one of:
		- 'HQ_HI'
		- 'HQ_LO'
		- 'LQ_HI'
		- 'LQ_LO'

	Classification is based on:
		- Amplicon Score cutoff (>= 1)
		- Barcode Rank cutoff (<= 10000)

	Parameters
	----------
	input_data : pandas.DataFrame
		Must contain:
			- 'Amplicon Score'
			- 'Barcode Rank'

	Returns
	-------
	pandas.DataFrame
		Same DataFrame with an added 'Color' column.

	Notes
	-----
	- Cutoffs are currently hard-coded.
	- Does not modify other columns.
	"""
	#amp_score_cutoff = 1
	depth_cutoff = 10000

	
	top_cells = input_data[input_data['Barcode Rank'] <= depth_cutoff]
	amp_score_cutoff = top_cells['Amplicon Score'].median()

	codes = []

	for i in range(len(input_data)):
		barcode_rank = input_data['Barcode Rank'].iloc[i]
		amp_score = input_data['Amplicon Score'].iloc[i]

		if barcode_rank <= depth_cutoff:
			if amp_score >= amp_score_cutoff:
				code = "HQ_HI"
			else:
				code = "LQ_HI"
		else:
			if amp_score >= amp_score_cutoff:
				code = "HQ_LO"
			else:
				code = "LQ_LO"

		codes.append(code)

	input_data['Color'] = codes
	return input_data 

def generate_amplicon_score(raw_tot_columns, min_reads_per_amplicon_per_cell, min_total_reads_per_barcode):
	"""
	Compute amplicon score and assign quality categories to barcodes.

	Amplicon score is computed by:
		- Evaluating percentile thresholds across amplicons.
		- Weighting contributions using predefined constants.
		- Summing across percentile tiers.

	Parameters
	----------
	raw_tot_columns : pandas.DataFrame
		DataFrame containing total read counts per amplicon per barcode.
	min_reads_per_amplicon_per_cell : int
		Minimum reads required at each amplicon.
	min_total_reads_per_barcode : int
		Minimum total reads required to retain a barcode.

	Returns
	-------
	pandas.DataFrame
		DataFrame indexed by barcode with columns:
			- 'Amplicon Score'
			- 'Read Count'
			- 'Barcode Rank'
			- 'Color' (quality category)

	Notes
	-----
	- Applies per-amplicon coverage filtering first.
	- Applies total-read filtering last.
	- Calls `stratify_data` to assign quality categories.
	"""
	percentile_cutoffs = [0.975, 0.99, 0.999, 0.9999]
	constant_values = [1, 10, 50, 100]
	
	# Filter to cells with 'min_reads_per_amplicon_per_cell' or more reads for all amplicons
	mask = (raw_tot_columns >= min_reads_per_amplicon_per_cell).all(axis = 1)
	#passing_barcodes = raw_tot_columns.index[mask]
	raw_tot_columns = raw_tot_columns.loc[mask]

	logging.info('Cells that did not pass the read count per amplicon cutoff:' + str(len(mask) - sum(mask)))
   
	# Filter
	# Create a list of the percentile values for each amplicon
	percentile_values = []
	for percentile in percentile_cutoffs:
		percentile_values.append(raw_tot_columns.quantile(percentile, axis = 0))
		
	# Create a dictionary to store the amplicon scores
	amplicon_dict = {}
	
	# Calculate the amplicon statistic
	for df_index, row in raw_tot_columns.iterrows():
		barcode_sum = 0
		for index in range(0, len(percentile_cutoffs)):
			# Get the percentile values for the current amplicon
			percentile_vals = percentile_values[index]
			# Get the constant value for the current percentile cutoff
			constant_val = constant_values[index]
			# Calculate the amplicon stat: 
			# percentage of amplicon values over the percentile value multiplied by a constant
			amplicon_stat = ((sum(row >= percentile_vals)) / len(percentile_vals)) * constant_val
			barcode_sum += amplicon_stat
		# Assign the amplicon score to the amplicon dictionary
		amplicon_dict[df_index] = barcode_sum
	
	amplicon_df = pd.DataFrame.from_dict(amplicon_dict, orient = 'index', columns = ['Amplicon Score'])
	raw_sum = raw_tot_columns.sum(axis = 1)
	amplicon_df['Read Count'] = raw_sum
	amplicon_df = amplicon_df.sort_values("Read Count", ascending = False)
	amplicon_df['Barcode Rank'] = range(1, len(amplicon_df) + 1)
	
	amplicon_df = stratify_data(amplicon_df)
	
	filtered_df = amplicon_df[amplicon_df['Read Count'] >= min_total_reads_per_barcode]
	#filtered_df = filtered_df[filtered_df.index.isin(passing_barcod
	
	
	return filtered_df

def add_color_information(editingSummary, color_df):
	"""
	Add quality category ('Color') column to editing summary DataFrame.

	Parameters
	----------
	editingSummary : pandas.DataFrame
		Editing summary DataFrame indexed by barcode.
	color_df : pandas.DataFrame
		DataFrame containing 'Color' classification indexed by barcode.

	Returns
	-------
	pandas.DataFrame
		Copy of `editingSummary` filtered to barcodes present in `color_df`,
		with an added 'Color' column.

	Notes
	-----
	- Barcodes not present in `color_df` are removed.
	- Does not modify original DataFrame.
	"""
	# Fillter editing summary to barcodes within color_df and add appropriate color
	# value to the editingSummary df    
	editingSummary = editingSummary[editingSummary.index.isin(color_df.index)]
	Color_col = [color_df.loc[barcode, 'Color'] for barcode in editingSummary.index]
	new_df = editingSummary.copy()
	new_df['Color'] = Color_col
	return new_df

def plot_amp_score(output_root):
	"""
	Generate a scatter plot of Amplicon Score versus Barcode Rank.

	Each barcode is plotted with:
		- X-axis: Barcode Rank (descending by total read count)
		- Y-axis: Amplicon Score

	Points are colored according to the 'Color' quality category
	assigned during amplicon score computation.

	Parameters
	----------
	output_root : str
		Base output prefix used to locate input summary files and
		write plot outputs.
		Required input file:
			- <output_root>.amplicon_score.txt

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Reads:
		- <output_root>.amplicon_score.txt

	Writes:
		- <output_root>.05_Amplicon_Score.pdf
		- <output_root>.05_Amplicon_Score.png

	Notes
	-----
	- Expected quality category short codes:
		{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.
	- The amplicon score is computed in `generate_amplicon_score`.
	- Barcode Rank is assigned after sorting by total read count.
	- Category counts are displayed directly on the plot.
	- Does not modify input data.
	"""
	plt.clf()
	plt.cla()

	_legacy_to_short = {
	"High_score_High_depth": "HQ_HI",
	"High_score_Low_depth":  "HQ_LO",
	"Low_score_High_depth":  "LQ_HI",
	"Low_score_Low_depth":   "LQ_LO",
	}
	
	_shortcode_to_color_and_label = {
		"HQ_HI": ("blue",   "HQ_HI"),  # High Score / High Depth
		"HQ_LO": ("green",  "HQ_LO"),  # High Score / Low Depth
		"LQ_HI": ("orange", "LQ_HI"),  # Low Score / High Depth
		"LQ_LO": ("red",    "LQ_LO"),  # Low Score / Low Depth
	}

	data = pd.read_csv(output_root + ".amplicon_score.txt", sep = "\t")
	value_counts = data['Color'].value_counts()

	colors = []
	categories = []
	for index in value_counts.index:
		short_code = _legacy_to_short.get(index,index)
		color_label = _shortcode_to_color_and_label.get(short_code, ("gray", short_code))
		colors.append(color_label[0])
		categories.append(color_label[1])
		  
	summary_stats = pd.DataFrame({"Count": value_counts, "Color": colors, "Category": categories})
	
	amp_score_max = data['Amplicon Score'].max()
   
	color_array = []
	for color_val in data['Color']:
		short_code = _legacy_to_short.get(color_val, color_val)
		color = _shortcode_to_color_and_label.get(short_code, ("gray", short_code))[0]
		color_array.append(color)
	
	x_loc = 0.9 * data['Barcode Rank'].max()
	
	# Set plot size
	plt.figure(figsize = (12,12))
	
	plt.scatter(data['Barcode Rank'],
				data['Amplicon Score'],
				color = color_array)
	plt.xlabel("Barcode Rank", fontsize = 22)
	plt.ylabel("Amplicon Score", fontsize = 22)
	plt.title("Amplicon Score", fontsize = 24)
	plt.yticks(fontsize = 16)
	plt.xticks(fontsize = 16)
	plt.tight_layout()
	
	for i in range(len(summary_stats)):
		plt.text(x_loc, (0.90 - (i * 0.05)) * amp_score_max,
				 str(summary_stats.iloc[i, 2]) + " Barcode Counts: " + str(summary_stats.iloc[i, 0]),
				 fontsize = 16, color = summary_stats.iloc[i, 1], ha = 'right')
	
	amp_plot_root = output_root + ".05_Amplicon_Score"
	plt.savefig(amp_plot_root + ".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(amp_plot_root + ".png", pad_inches = 1, bbox_inches = "tight")
	   
	logging.info("Finished amplicon score plot")
	summary_plot_obj = PlotObject(
		plot_name = amp_plot_root,
		plot_title = 'Amplicon Score Plot',
		plot_label = "Plot of the Barcode Rank (X) vs. Amplicon Score (Y) with color coding according to Amplicon Score category.",
		plot_datas = [
			("Amplicon Score", output_root + ".amplicon_score.txt")
		]
	) 
	return summary_plot_obj 
	 
def log_log_plot(parsed_information, output_root, cell_quality_to_analyze, filtered = True):
	"""
	Generate a log-log scatter plot of read count versus barcode rank.

	For each barcode:
		- Total read count is computed as the sum of all 'totCount.*' columns.
		- Barcodes are sorted in descending order of total read count.
		- Barcode rank is assigned accordingly (1 = highest read count).

	The plot displays:
		- X-axis: Barcode Rank (log scale)
		- Y-axis: Total Read Count (log scale)
		- Point color: Quality category ('Color' column)

	If `filtered` is True, only barcodes whose 'Color' value is in
	`cell_quality_to_analyze` are plotted.

	Parameters
	----------
	parsed_information : pandas.DataFrame
		Editing summary DataFrame containing:
			- 'totCount.*' columns
			- 'Color' column with quality category codes

	output_root : str
		Base output prefix used to write plot outputs.

	cell_quality_to_analyze : list[str]
		List of quality category short codes used for filtering when
		`filtered` is True.
		Expected values:
			{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.

	filtered : bool, optional
		If True, restrict plot to barcodes matching
		`cell_quality_to_analyze`.
		If False, plot all barcodes.

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Writes:
		- <output_root>.01_Log-Log.pdf
		- <output_root>.01_Log-Log.png

		or, if filtered is True:

		- <output_root>.01_Log-Log_filtered.pdf
		- <output_root>.01_Log-Log_filtered.png

	Notes
	-----
	- Log scaling is applied to both axes.
	- Quality categories are expected to use short codes:
		{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.
	- Does not modify input DataFrame.
	"""
	plt.cla()
	plt.clf()
	# Fix column names
	parsed_information.columns = parsed_information.columns.str.replace('[^a-zA-Z0-9]', '_')
	colors = parsed_information['Color']
	COLOR_DISPLAY_MAP = {
		"HQ_HI": "blue",
		"HQ_LO": "green",
		"LQ_HI": "yellow",
		"LQ_LO": "red"
	}
	# COLOR_DISPLAY_MAP = {
		# "HQ_HI": ("High Score / High Reads", "green"),
		# "HQ_LO": ("High Score / Low Reads", "blue"),
		# "LQ_HI": ("Low Score / High Reads", "orange"),
		# "LQ_LO": ("Low Score / Low Reads", "red"),
	# }
	
	
	colors = [COLOR_DISPLAY_MAP[color] for color in colors]
	
	totCols = parsed_information.filter(like = "totCount")
	rowsum_data = {"Barcode": totCols.index,
			"Read Count": totCols.sum(axis = 1)}

	rowsum_DF = pd.DataFrame(rowsum_data)

	rowsum_DF.index = range(1, len(totCols.index) + 1)
	rowsum_DF['Color'] = colors
	rowsum_DF = rowsum_DF.sort_values(by = "Read Count", ascending=False)
	rowsum_DF['Barcode Count'] = range(1, len(totCols.index) + 1)
   
	#barcodes = amplicon[amplicon['Color'].isin(cell_quality_to_analyze)].index.tolist()

	if filtered:
		color_picker = [COLOR_DISPLAY_MAP[key] for key in cell_quality_to_analyze]
		#print(f"{color_picker=}") 
		rowsum_DF = rowsum_DF[rowsum_DF['Color'].isin(color_picker)]
	
	# Main scatter plot
	plt.figure(figsize = (12,12))
	# Red points
	red_points = rowsum_DF[rowsum_DF['Color'] == 'red']
	plt.scatter(red_points['Barcode Count'], red_points['Read Count'], color='red')
	
	# Yellow points
	yellow_points = rowsum_DF[rowsum_DF['Color'] == 'yellow']
	plt.scatter(yellow_points['Barcode Count'], yellow_points['Read Count'], color='yellow')
	
	# Green Points
	green_points = rowsum_DF[rowsum_DF['Color'] == 'green']
	plt.scatter(green_points['Barcode Count'], green_points['Read Count'], color='green')
	
	# Blue Points
	blue_points = rowsum_DF[rowsum_DF['Color'] == 'blue']
	plt.scatter(blue_points['Barcode Count'], blue_points['Read Count'], color='blue')

	plt.xlabel("Barcode Rank", fontsize = 22)
	plt.ylabel("Read Count", fontsize = 22)
	plt.xscale("log")
	plt.yscale("log")
	plt.title("Barcode Rank vs. Read Count (log scale)", fontsize = 24)
	plt.xticks(fontsize = 16)
	plt.yticks(fontsize = 16)
	
	legend_elements = [mpatches.Patch(color=color, label=label + ": " + str(len(rowsum_DF[rowsum_DF['Color'] == color])) + " Barcodes") for label, color in COLOR_DISPLAY_MAP.items()]
	plt.legend(handles=legend_elements, fontsize = 16)
	
	if filtered:
		log_log_root = output_root + ".01_Log-Log_filtered"
	else:
		log_log_root = output_root + ".01_Log-Log"
		
	plt.savefig(log_log_root + ".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(log_log_root + ".png", pad_inches = 1, bbox_inches = "tight")
	   
	logging.info("Finished log-log plot")
	summary_plot_obj = PlotObject(
		plot_name = log_log_root,
		plot_title = 'Log-Log Plot',
		plot_label = "Log scale plot of the Barcode Rank (X) vs. Read Count (Y) with color coding according to Amplicon Score category.",
		plot_datas = [
			
		]
	)
	
	return summary_plot_obj 
		
def cell_per_amp_filtered(parsed_information, output_root, cell_quality_to_analyze):
	"""
	Generate a plot of the number of barcodes covering each amplicon
	at multiple read-count thresholds.

	For barcodes whose 'Color' value is in `cell_quality_to_analyze`,
	this function counts, for each amplicon, the number of barcodes
	with read counts greater than or equal to each threshold in:

		[1, 5, 10, 25, 50, 100]

	The result is plotted as a point plot showing how many cells
	sufficiently cover each amplicon under increasing coverage stringency.

	Parameters
	----------
	parsed_information : pandas.DataFrame
		Editing summary DataFrame containing:
			- 'totCount.*' columns
			- 'Color' column

	output_root : str
		Base output prefix used to write plot outputs.

	cell_quality_to_analyze : list[str]
		List of quality category short codes to include.
		Expected values:
			{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Writes:
		- <output_root>.02_CellCountPerAmplicon_filtered.pdf
		- <output_root>.02_CellCountPerAmplicon_filtered.png

	Notes
	-----
	- Only barcodes matching `cell_quality_to_analyze` are included.
	- Coverage is defined as totCount >= threshold.
	- Amplicons are ordered by mean coverage across thresholds.
	- Does not modify input DataFrame.
	"""
	parsed_information.columns = parsed_information.columns.str.replace('[^a-zA-Z0-9]', '_')

	parsed_information = parsed_information[parsed_information['Color'].isin(cell_quality_to_analyze)]

	# Filter for only high score / high reads and high score / low reads
	#parsed_information = parsed_information[parsed_information['Color'].isin(["High Score / High Reads", "High Score / Low Reads"])]

	# Grab total count columns
	totCols = parsed_information.filter(like = "totCount")

	# Get counts of cells with different read cutoffs
	read_cutoffs = [1, 5, 10, 25, 50, 100]
	counts = {f'{c} reads': (parsed_information[totCols.columns] >= c).sum() for c in read_cutoffs}

	# Create a DataFrame with counts
	tots = pd.DataFrame(counts).T
	tots.columns = tots.columns.str.replace('totCount.', '')

	# Sort columns by mean values
	amplicon_means = tots.mean().sort_values(ascending=False).index
	tots = tots[amplicon_means]

	# Reshape the data for plotting
	tots['Read_cutoff'] = tots.index
	tots2 = tots.melt(id_vars='Read_cutoff', var_name='Targets', value_name='Cell_count')
	tots2['Read_cutoff'] = pd.Categorical(tots2['Read_cutoff'], categories=tots.index)
	tots2['Targets'] = pd.Categorical(tots2['Targets'], categories=amplicon_means)

	# Create a line plot using Matplotlib
	plt.figure(figsize=(12, 12))
	sns.pointplot(data = tots2, x = 'Targets', y = 'Cell_count', hue = 'Read_cutoff', palette = 'bright')
	plt.title('Cell count per amplicon with minimum specified coverage', fontsize = 24)
	plt.yticks(fontsize = 16)
	plt.xticks(rotation=0 if len(totCols) < 50 else 90, fontsize = 16)
	plt.xlabel('Amplicons', fontsize = 22)
	plt.ylabel('Cell count', fontsize = 22)
	plt.tight_layout()
	cell_per_amp_root = output_root + ".02_CellCountPerAmplicon_filtered"
	plt.savefig(cell_per_amp_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(cell_per_amp_root+".png", pad_inches = 1, bbox_inches = "tight")

	logging.info("Finished cell count per amplicon plot")
	summary_plot_obj = PlotObject(
			plot_name = cell_per_amp_root,
			plot_title = 'Cell count per amplicon with minimum specified coverage',
			plot_label = 'Plotting the number of cells covering an amplicon at a given read cutoff. The log cell counts are on the Y-axis and the amplicon is on the X-axis.',
			plot_datas = [
				]
			)
	return summary_plot_obj
 

def amp_per_cell_filtered(parsed_information, output_root, cell_quality_to_analyze):
	"""
	Generate a plot of the number of amplicons covered per barcode
	at multiple read-count thresholds.

	For barcodes whose 'Color' value is in `cell_quality_to_analyze`,
	this function computes, for each barcode:

		- The number of amplicons with totCount >= threshold,
		  where threshold ∈ [1, 5, 10, 25, 50, 100].

	It then determines how many barcodes achieve coverage of
	100%, 90%, 80%, and 50% of total amplicons at each threshold,
	and visualizes the results as a line plot.

	Parameters
	----------
	parsed_information : pandas.DataFrame
		Editing summary DataFrame containing:
			- 'totCount.*' columns
			- 'Color' column

	output_root : str
		Base output prefix used to write plot outputs.

	cell_quality_to_analyze : list[str]
		List of quality category short codes to include.
		Expected values:
			{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Writes:
		- <output_root>.03_AmpliconCoveredPerCell_filtered.pdf
		- <output_root>.03_AmpliconCoveredPerCell_filtered.png

	Notes
	-----
	- Only barcodes matching `cell_quality_to_analyze` are included.
	- Coverage is defined as totCount >= threshold.
	- Coverage percentages are relative to total number of amplicons.
	- Does not modify input DataFrame.
	"""
	parsed_information.columns = parsed_information.columns.str.replace('[^a-zA-Z0-9]', '_')
	#parsed_information = parsed_information[parsed_information['Color'].isin(["High Score / High Reads", "High Score / Low Reads"])]
	parsed_information = parsed_information[parsed_information['Color'].isin(cell_quality_to_analyze)]    

	read_cutoffs = [1,5,10,25,50,100]
	# Grab total count columns
	cov_cols = [col for col in parsed_information.columns if col.startswith("totCount")]
	# Row sum of cells with amplicon coverage at specified cutoffs
	g1 = (parsed_information[cov_cols] >= 1).sum(axis=1) 
	g5 = (parsed_information[cov_cols] >= 5).sum(axis=1)
	g10 = (parsed_information[cov_cols] >= 10).sum(axis=1)
	g25 = (parsed_information[cov_cols] >= 25).sum(axis=1)
	g50 = (parsed_information[cov_cols] >= 50).sum(axis=1)
	g100 = (parsed_information[cov_cols] >= 100).sum(axis=1)

	# Determine number of amplicons for 90%, 80%, and 50% coverage of target amplicons
	cov_100pct = len(cov_cols)
	cov_90pct = round(len(cov_cols) * 0.9)
	cov_80pct = round(len(cov_cols) * 0.8)
	cov_50pct = round(len(cov_cols) * 0.5)
	break_vals = [cov_100pct, cov_90pct, cov_80pct, cov_50pct]

	vals1 = [len(g1[g1 > cov_100pct]), len(g1[g1 > cov_90pct]), len(g1[g1 > cov_80pct]), len(g1[g1 > cov_50pct])]
	vals5 = [len(g5[g5 > cov_100pct]), len(g5[g5 > cov_90pct]), len(g5[g5 > cov_80pct]), len(g5[g5 > cov_50pct])]
	vals10 = [len(g10[g10 > cov_100pct]), len(g10[g10 > cov_90pct]), len(g10[g10 > cov_80pct]), len(g10[g10 > cov_50pct])]
	vals25 = [len(g25[g25 > cov_100pct]), len(g25[g25 > cov_90pct]), len(g25[g25 > cov_80pct]), len(g25[g25 > cov_50pct])]
	vals50 = [len(g50[g50 > cov_100pct]), len(g50[g50 > cov_90pct]), len(g50[g50 > cov_80pct]), len(g50[g50 > cov_50pct])]
	vals100 = [len(g100[g100 > cov_100pct]), len(g100[g100 > cov_90pct]), len(g100[g100 > cov_80pct]), len(g100[g100 > cov_50pct])]

	# Create DataFrame structure
	tots = pd.DataFrame([vals1, vals5, vals10, vals25, vals50, vals100], 
					columns=["Target_" + str(val) for val in break_vals], 
					index=[str(x) + " reads" for x in read_cutoffs])

	# Preserve the intended X-axis order: highest number of covered targets
	# to lowest, rather than re-ordering buckets by the observed cell counts.
	target_order = list(tots.columns)

	# Reshape for plotting
	tots['Read_cutoff'] = [str(x) + " reads" for x in read_cutoffs]
	tots2 = tots.melt(id_vars = "Read_cutoff", var_name = "Targets", value_name = "Cell_count")

	tots2['Targets'] = pd.Categorical(tots2['Targets'], categories=target_order, ordered=True)
	tots2 = tots2.sort_values(by = "Targets")

	# Plot
	plt.figure(figsize=(12, 12))
	sns.lineplot(data=tots2, x='Targets', y='Cell_count', 
				hue='Read_cutoff', 
				markers=True, palette='bright')

	# Order the legend based on read cutoffs
	handles, labels = plt.gca().get_legend_handles_labels()
	sorted_handles_labels = sorted(zip(handles, labels), key=lambda x: int(x[1].split()[0]))
	handles, labels = zip(*sorted_handles_labels)
	plt.legend(handles, labels)
	
	plt.title('Amplicons covered per cell with minimum specified coverage', fontsize = 24)
	plt.xticks(ticks = tots2['Targets'], labels = tots2['Targets'].str.replace('Target_', ''), fontsize = 16)
	plt.yticks(fontsize = 16)
	plt.xlabel('Number of amplicon targets covered', fontsize = 16)
	plt.ylabel('Cell count', fontsize = 16)
	plt.tight_layout()
	amp_per_cell_obj_root = output_root + ".03_AmpliconCoveredPerCell_filtered"
	plt.savefig(amp_per_cell_obj_root+".png", pad_inches=1, bbox_inches='tight')
	plt.savefig(amp_per_cell_obj_root+".pdf", pad_inches=1, bbox_inches='tight') 

	summary_plot_obj = PlotObject(
		plot_name = amp_per_cell_obj_root,
		plot_title = "Amplicons covered per cell with minimum specified coverage",
		plot_label = "The number of cells with amplicon coverage at a given read cutoff. The number of amplicons covered is on the X-axis and the cell count is on the Y-axis.",
		plot_datas = [
		]
	)
	logging.info("Finished Amplicon coverage per cell plot.")
	
	return summary_plot_obj


def mod_per_amp_filtered(parsed_information, output_root, cell_quality_to_analyze):
	"""
	Generate a plot of average modification percentage per amplicon
	at multiple read-count thresholds.

	For barcodes whose 'Color' value is in `cell_quality_to_analyze`,
	this function computes, for each amplicon and each threshold
	in [1, 5, 10, 25, 50, 100]:

		- The mean modification percentage ('modPct.*')
		  among barcodes with totCount >= threshold.

	The results are visualized as a line plot showing how
	average modification varies with coverage stringency.

	Parameters
	----------
	parsed_information : pandas.DataFrame
		Editing summary DataFrame containing:
			- 'totCount.*' columns
			- 'modPct.*' columns
			- 'Color' column

	output_root : str
		Base output prefix used to write plot outputs.

	cell_quality_to_analyze : list[str]
		List of quality category short codes to include.
		Expected values:
			{'HQ_HI','HQ_LO','LQ_HI','LQ_LO'}.

	Returns
	-------
	PlotObject
		Metadata object describing generated plot and associated data files.

	Side Effects
	------------
	Writes:
		- <output_root>.04_ModPercentagePerAmp_filtered.pdf
		- <output_root>.04_ModPercentagePerAmp_filtered.png

	Notes
	-----
	- Only barcodes matching `cell_quality_to_analyze` are included.
	- Coverage threshold is applied before computing modification mean.
	- Does not modify input DataFrame.
	"""
	parsed_information.columns = parsed_information.columns.str.replace('[^a-zA-Z0-9]', '_')
	parsed_information = parsed_information[parsed_information['Color'].isin(cell_quality_to_analyze)]
	#parsed_information = parsed_information[parsed_information['Color'].isin(["High Score / High Reads", "High Score / Low Reads"])]
	
	read_cutoffs = [1,5,10,25,50,100]

	colors = cell_quality_to_analyze
	#colors = ['High Score / High Reads', 'High Score / Low Reads']
	PerMod_df = pd.DataFrame(columns=['Read_cutoff', 'Mod_average', 'Target', 'Color'])        
		
	for color in colors:
		for cutoff in read_cutoffs:
			totCols = parsed_information[parsed_information['Color'] == color].filter(like = 'totCount')
			modCols = parsed_information[parsed_information['Color'] == color].filter(like = 'modPct')
			for i in range(0, totCols.shape[1]):
				totCol = totCols.iloc[:, i]
				modCol = modCols.iloc[:, i]
				modAvg = modCol[totCol >= cutoff].mean()
				col_title = totCol.name.replace('totCount.', '')
				PerMod_df.loc[len(PerMod_df)] = [cutoff, modAvg, col_title, color]
				
	# Create a combined average df for ordering in the plot
	Modification_Average_df = PerMod_df.groupby(['Target'])['Mod_average'].mean()
	Modification_Average_df = Modification_Average_df.sort_values(ascending=False)
	Modification_Average_df = Modification_Average_df.reset_index()

	# Order PerMod_df amplicons according to overall modification average
	PerMod_df['Target'] = pd.Categorical(PerMod_df['Target'], categories= Modification_Average_df['Target'])
	PerMod_df = PerMod_df.sort_values(by = "Target")
	
	# Generate Plot
	sns.set_style("white")
	plt.figure(figsize=(12, 12)) # 6,6
	sns.lineplot(data=PerMod_df, x='Target', y='Mod_average', 
				hue='Read_cutoff', 
				markers=True, palette='bright', errorbar = None)
	plt.gca().patch.set_alpha(0)
	plt.title('Average modification percentage of an amplicon with minimum specified coverage', fontsize = 24)
	plt.xticks(rotation=90, fontsize = 16)
	plt.yticks(fontsize = 16)
	plt.xlabel('Amplicons', fontsize = 22)
	plt.ylabel('Modification Percentage', fontsize = 22)
	plt.tight_layout()
	
	mod_pct_plot_obj_root = output_root + ".04_ModPercentagePerAmp_filtered" 
	plt.savefig(mod_pct_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
	plt.savefig(mod_pct_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')
	logging.info("Finished modification percentage per amplicon plot.")
	summary_plot_obj = PlotObject(
			plot_name = mod_pct_plot_obj_root,
			plot_title = 'Average modification by target',
			plot_label = 'The average modification percentage (Y) of an amplicon (X) with a minimum specified read coverage.',
			plot_datas = [
				]
			)
	return summary_plot_obj



### EXAMPLE of creation of plot object

#    fig = plt.figure(figsize=(12,12))
#    ax = plt.subplot(111)
#    pie_values = []
#    pie_labels = []
#    for i in range(len(filter_labels)):
#        if values[i] > 0:
#            pie_values.append(values[i])
#            pie_labels.append(filter_labels[i]+"\n("+str(values[i])+")")
#    ax.pie(pie_values, labels=pie_labels, autopct="%1.2f%%")
#    ax.set_title('Read Assignment Summary')
#    plt.savefig(summary_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
#    plt.savefig(summary_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')
#
#    plot_count_str = "<br>".join(["%s N=%s"%x for x in zip(filter_labels,values)])
#    plot_label = 'Assignment summary for N=' + str(num_reads_input) + ' input reads:'
#    summary_plot_obj = PlotObject(
#            plot_name = summary_plot_obj_root,
#            plot_title = 'Read summary',
#            plot_label = plot_label + '<br>'+plot_count_str,
#            plot_datas = [
#                ('Read assignment summary',summary_plot_obj_root + ".txt")
#                ]
#            )

#### END OF EXAMPLE


class PlotObject:
	"""
	Holds information for plots for future output, namely:
		the plot name: root of plot (name.pdf and name.png should exist)
		the plot title: title to be shown to user
		the plot label: label to be shown under the plot
		the plot data: array of (tuple of display name and file name)
		the plot order: int specifying the order to display on the report (lower numbers are plotted first, followed by higher numbers)
	"""
	def __init__(self,plot_name,plot_title,plot_label,plot_datas,plot_order=50):
		self.name = plot_name
		self.title = plot_title
		self.label = plot_label
		self.datas = plot_datas
		self.order = plot_order

	def to_json(self):
		obj = {
				'plot_name':self.name,
				'plot_title':self.title,
				'plot_label':self.label,
				'plot_datas':self.datas,
				'plot_order':self.order
				}
		obj_str = json.dumps(obj,separators=(',',':'))
		return obj_str

	#construct from a json string
	@classmethod
	def from_json(cls, json_str):
		obj = json.loads(json_str)
		return cls(plot_name=obj['plot_name'],
				plot_title=obj['plot_title'],
				plot_label=obj['plot_label'],
				plot_datas=obj['plot_datas'],
				plot_order=obj['plot_order'])

	def __str__(self):
		return 'Plot object with name ' + self.name

	def __repr__(self):
		return f'PlotObject(name={self.name}, title={self.title}, label={self.label}, datas={self.datas} order={self.order})'



def make_report(report_file,report_name,results_folder,
			crispresso_run_names,crispresso_sub_html_files,
			summary_plot_objects=[]
		):
	"""
	Makes an HTML report for a CRISPRSCope run

	Parameters:
	report_file: path to the output report
	report_name: description of report type to be shown at top of report
	results_folder (string): absolute path to the CRISPRSCope output

	crispresso_run_names (arr of strings): names of crispresso runs
	crispresso_sub_html_files (dict): dict of run_name->file_loc

	summary_plot_objects (list): list of PlotObjects to plot
	"""

	logger = logging.getLogger()
	ordered_plot_objects = sorted(summary_plot_objects,key=lambda x: x.order)

	html_str = """
<!doctype html>
<html lang="en">
  <head>
	<meta charset="utf-8">
	<meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
	<title>"""+report_name+"""</title>

	<!-- Bootstrap core CSS -->
<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootswatch/4.5.2/flatly/bootstrap.min.css" integrity="sha384-qF/QmIAj5ZaYFAeQcrQ6bfVMAh4zZlrGwTPY7T/M+iTTLJqJBJjwwnsE5Y0mV7QK" crossorigin="anonymous">
  </head>

  <body>
<style>
html,
body {
  height: 100%;
}

body {
  padding-bottom: 40px;
  background-color: #f5f5f5;
}

.navbar-fixed-left {
  width: 200px;
  position: fixed;
  border-radius: 0;
  height: 100%;
  padding: 10px;
}

.navbar-fixed-left .navbar-nav > li {
  /*float: none;   Cancel default li float: left */
  width: 160px;
}

</style>

<nav class="navbar navbar-fixed-left navbar-dark bg-dark" style="overflow-y:auto">
	 <a class="navbar-brand" href="#">CRISPRSCope</a>
	  <ul class="nav navbar-nav me-auto">
"""
	for idx,plot_obj in enumerate(ordered_plot_objects):
		html_str += """        <li class="nav-item">
		  <a class="nav-link active" href="#plot"""+str(idx)+"""">"""+plot_obj.title+"""
		  </a>
		</li>"""
	if len(crispresso_run_names) > 0:
		html_str += """        <li class="nav-item">
		  <a class="nav-link active" href="#crispresso_output">CRISPResso Output
		  </a>
		</li>"""
	html_str += """      </ul>
</nav>
<div class='container'>
<div class='row justify-content-md-center'>
<div class='col-8'>
	<div class='text-center pb-4'>
	<h1 class='display-3 pt-5'>CRISPRSCope</h1><hr><h2>"""+report_name+"""</h2>
	</div>
"""

	data_path = ""
	for idx,plot_obj in enumerate(ordered_plot_objects):
		plot_path = plot_obj.name
		plot_path = os.path.basename(plot_path)
		plot_str = "<div class='card text-center mb-2' id='plot"+str(idx)+"'>\n\t<div class='card-header'>\n"
		plot_str += "<h5>"+plot_obj.title+"</h5>\n"
		plot_str += "</div>\n"
		plot_str += "<div class='card-body'>\n"
		plot_str += "<a href='"+data_path + plot_path+".pdf'><img src='"+data_path + plot_path + ".png' width='80%' ></a>\n"
		plot_str += "<label>"+plot_obj.label+"</label>\n"
		for (plot_data_label,plot_data_path) in plot_obj.datas:
			plot_data_path = os.path.basename(plot_data_path)
			plot_str += "<p class='m-0'><small>Data: <a href='"+data_path+plot_data_path+"'>" + plot_data_label + "</a></small></p>\n"
		plot_str += "</div></div>\n"
		html_str += plot_str

	if len(crispresso_run_names) > 0:
		run_string = """<div class='card text-center mb-2' id='crispresso_output'>
		  <div class='card-header'>
			<h5>CRISPResso Output</h5>
		  </div>
		  <div class='card-body p-0'>
			<div class="list-group list-group-flush">
			"""
		for crispresso_run_name in crispresso_run_names:
			run_string += "<a href='"+data_path+crispresso_sub_html_files[crispresso_run_name]+"' class='list-group-item list-group-item-action'>"+crispresso_run_name+"</a>\n"
		run_string += "</div></div></div>"
		html_str += run_string

	html_str += """
				</div>
			</div>
		</div>
	</body>
</html>
"""
	with open(report_file,'w') as fo:
		fo.write(html_str)
	logger.info('Wrote ' + report_file)
if __name__ == "__main__":
	main()
