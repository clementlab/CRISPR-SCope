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
from CRISPRSCope.io_utils import open_text_maybe_gzip
from CRISPRSCope.output_artifacts import OutputContext, OutputManifest

from .settings import CELL_QUALITY_CODES, AmpliconScoreConfig

def _numeric_tot_count_columns(df):
	"""
	Return total-count columns coerced to numeric values.
	"""
	return df.filter(like = "totCount").apply(pd.to_numeric, errors = "coerce")


def _numeric_mod_pct_columns(df):
	"""
	Return modification-percentage columns coerced to numeric values.
	"""
	return df.filter(like = "modPct").apply(pd.to_numeric, errors = "coerce")


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
	editingSummary_count = _numeric_tot_count_columns(editing)

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

	amplicon_cov_plot_root = OutputContext(output_root).plot_root("amplicon_coverage_plot")
	plt.savefig(amplicon_cov_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(amplicon_cov_plot_root+".png", pad_inches = 1, bbox_inches = "tight")

	logging.info("Finished generating the average amplicon coverage plot.")
	summary_plot_obj = declared_plot_object(output_root, "amplicon_coverage_plot")
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

	high_totals = _numeric_tot_count_columns(editing[editing.index.isin(high_barcodes)]).sum(axis=1)
	low_totals = _numeric_tot_count_columns(editing[editing.index.isin(low_barcodes)]).sum(axis=1)

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

	cell_cov_plot_root = OutputContext(output_root).plot_root("cell_coverage_boxplot")
	plt.savefig(cell_cov_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(cell_cov_plot_root+".png", pad_inches = 1, bbox_inches = "tight")

	logging.info("Finished generating the read counts per barcode barplot plot.")
	summary_plot_obj = declared_plot_object(output_root, "cell_coverage_boxplot")
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

	high_qual_edit = _numeric_tot_count_columns(editing[editing.index.isin(high_barcodes)])
	low_qual_edit = _numeric_tot_count_columns(editing[editing.index.isin(low_barcodes)])

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


	cell_cov_plot_root = OutputContext(output_root).plot_root("cell_coverage_plot")
	plt.savefig(cell_cov_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(cell_cov_plot_root+".png", pad_inches = 1, bbox_inches = "tight")

	logging.info("Finished generating the read counts per barcode barplot plot.")
	summary_plot_obj = declared_plot_object(output_root, "cell_coverage_plot")
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
	# Read in the final filtered editingSummary file
	editing = pd.read_csv(output_root + ".filteredEditingSummary.txt", sep = "\t", index_col = 0)
	# Read in amplicon score file
	amplicon = pd.read_csv(output_root + ".amplicon_score.txt", sep = "\t", index_col = 0)

	barcodes = amplicon[amplicon['Color'].isin(cell_quality_to_analyze)].index.tolist()
	editing = editing[editing.index.isin(barcodes)]

	editing = _numeric_mod_pct_columns(editing)

	editing = editing > 0

	row_sum = editing.sum(axis = 1)
	if row_sum.empty:
		logging.warning("Skipping edit histogram because no selected cells were found")
		return None

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

	hist_plot_root = OutputContext(output_root).plot_root("edit_histogram_plot")
	plt.savefig(hist_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(hist_plot_root+".png", pad_inches = 1, bbox_inches = "tight")

	logging.info("Finished generating the edit count histogram plot")
	summary_plot_obj = declared_plot_object(output_root, "edit_histogram_plot")
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
	# Read in the final filtered editingSummary file
	editing = pd.read_csv(output_root + ".filteredEditingSummary.txt", sep = "\t", index_col = 0)
	# Read in amplicon score file
	amplicon = pd.read_csv(output_root + ".amplicon_score.txt", sep = "\t", index_col = 0)

	# Filter barcodes if required
	barcodes = amplicon[amplicon['Color'].isin(cell_quality_to_analyze)].index.tolist()
	editing = editing[editing.index.isin(barcodes)]


	# Reduce to modified columns
	editing = _numeric_mod_pct_columns(editing)
	editing.columns = editing.columns.str.replace('modPct.', '')
	editing.fillna(0,inplace=True)

	# Convert to a boolean matrix
	editing = editing != 0

	# Find top 5 combinations of edit values
	counts = editing.apply(tuple, axis = 1).value_counts()
	top5 = counts.nlargest(5)

	# Filter editing matrix to top 5 combinations
	editing = editing[editing.apply(tuple, axis = 1).isin(top5.index)]

	# Remove columns where all modification values == False
	editing = editing.loc[:, ~(editing == False).all()]

	if editing.shape[1] == 0:
		logging.warning("Skipping edit combination upset plot because no edited sites were found")
		return None
	if editing.shape[1] == 1:
		logging.warning("Skipping edit combination upset plot because only one edited site was found")
		return None

	editing.index = pd.MultiIndex.from_frame(editing.astype(bool))

	upset = UpSet(editing, orientation = "horizontal", sort_by = "cardinality", show_counts = True)

	fig, ax = plt.subplots(figsize = (12, 12))
	upset.plot(fig = fig)

	for spine in ax.spines.values():
		spine.set_visible(False)
		ax.tick_params(left = False, bottom = False, labelleft = False, labelbottom = False, labelsize = 16)
	plt.title("Editing Sites and Intersections", fontsize=24)


	upset_plot_root = OutputContext(output_root).plot_root("edit_combinations_plot")
	plt.savefig(upset_plot_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(upset_plot_root+".png", pad_inches = 1, bbox_inches = "tight")

	logging.info("Finished generating the edit combination upset plot")
	summary_plot_obj = declared_plot_object(output_root, "edit_combinations_plot")
	return summary_plot_obj


def plot_amp_score(output_root, config=None):
	"""
	Generate a scatter plot of supported amplicon breadth versus barcode rank.

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
	- The amplicon score is the fraction of usable amplicons meeting the
	  configured support threshold.
	- Barcode Rank is assigned after sorting by total read count.
	- Category counts are displayed directly on the plot.
	- Does not modify input data.
	"""
	plt.clf()
	plt.cla()
	if config is None:
		config = AmpliconScoreConfig()

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

	amp_score_max = max(float(data['Amplicon Score'].max()), 1.0)

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
	plt.axhline(
		config.min_covered_fraction,
		color='black',
		linestyle='--',
		linewidth=1.5,
		label='Supported-breadth threshold',
	)
	plt.axvline(
		config.max_barcode_rank,
		color='gray',
		linestyle=':',
		linewidth=1.5,
		label='High-depth rank threshold',
	)
	plt.xlabel("Barcode Rank", fontsize = 22)
	plt.ylabel(
		f"Supported Amplicon Fraction (≥{config.min_reads_per_amplicon} reads)",
		fontsize = 22,
	)
	plt.title("Supported Amplicon Breadth", fontsize = 24)
	plt.ylim(-0.02, 1.02)
	plt.legend(loc='lower left', fontsize=12)
	plt.yticks(fontsize = 16)
	plt.xticks(fontsize = 16)
	plt.tight_layout()

	for i in range(len(summary_stats)):
		plt.text(x_loc, (0.90 - (i * 0.05)) * amp_score_max,
				 str(summary_stats.iloc[i, 2]) + " Barcode Counts: " + str(summary_stats.iloc[i, 0]),
				 fontsize = 16, color = summary_stats.iloc[i, 1], ha = 'right')

	amp_plot_root = OutputContext(output_root).plot_root("amplicon_score_plot")
	plt.savefig(amp_plot_root + ".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(amp_plot_root + ".png", pad_inches = 1, bbox_inches = "tight")

	logging.info("Finished amplicon score plot")
	summary_plot_obj = declared_plot_object(output_root, "amplicon_score_plot")
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

	totCols = _numeric_tot_count_columns(parsed_information)
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
		log_log_root = OutputContext(output_root).plot_root("log_log_filtered_plot")
		log_log_artifact_key = "log_log_filtered_plot"
	else:
		log_log_root = OutputContext(output_root).plot_root("log_log_plot")
		log_log_artifact_key = "log_log_plot"

	plt.savefig(log_log_root + ".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(log_log_root + ".png", pad_inches = 1, bbox_inches = "tight")

	logging.info("Finished log-log plot")
	summary_plot_obj = declared_plot_object(output_root, log_log_artifact_key)

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
	totCols = _numeric_tot_count_columns(parsed_information)

	# Get counts of cells with different read cutoffs
	read_cutoffs = [1, 5, 10, 25, 50, 100]
	counts = {f'{c} reads': (totCols >= c).sum() for c in read_cutoffs}

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
	cell_per_amp_root = OutputContext(output_root).plot_root("cell_count_per_amplicon_plot")
	plt.savefig(cell_per_amp_root+".pdf", pad_inches = 1, bbox_inches = "tight")
	plt.savefig(cell_per_amp_root+".png", pad_inches = 1, bbox_inches = "tight")

	logging.info("Finished cell count per amplicon plot")
	summary_plot_obj = declared_plot_object(output_root, "cell_count_per_amplicon_plot")
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
	totCols = _numeric_tot_count_columns(parsed_information)
	cov_cols = list(totCols.columns)
	# Row sum of cells with amplicon coverage at specified cutoffs
	g1 = (totCols >= 1).sum(axis=1)
	g5 = (totCols >= 5).sum(axis=1)
	g10 = (totCols >= 10).sum(axis=1)
	g25 = (totCols >= 25).sum(axis=1)
	g50 = (totCols >= 50).sum(axis=1)
	g100 = (totCols >= 100).sum(axis=1)

	# Determine number of amplicons for 90%, 80%, and 50% coverage of target amplicons
	cov_100pct = len(cov_cols)
	cov_90pct = round(len(cov_cols) * 0.9)
	cov_80pct = round(len(cov_cols) * 0.8)
	cov_50pct = round(len(cov_cols) * 0.5)
	break_vals = list(dict.fromkeys([cov_100pct, cov_90pct, cov_80pct, cov_50pct]))

	vals1 = [len(g1[g1 > val]) for val in break_vals]
	vals5 = [len(g5[g5 > val]) for val in break_vals]
	vals10 = [len(g10[g10 > val]) for val in break_vals]
	vals25 = [len(g25[g25 > val]) for val in break_vals]
	vals50 = [len(g50[g50 > val]) for val in break_vals]
	vals100 = [len(g100[g100 > val]) for val in break_vals]

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
	amp_per_cell_obj_root = OutputContext(output_root).plot_root("amplicon_covered_per_cell_plot")
	plt.savefig(amp_per_cell_obj_root+".png", pad_inches=1, bbox_inches='tight')
	plt.savefig(amp_per_cell_obj_root+".pdf", pad_inches=1, bbox_inches='tight')

	summary_plot_obj = declared_plot_object(output_root, "amplicon_covered_per_cell_plot")
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
		color_rows = parsed_information[parsed_information['Color'] == color]
		totCols = _numeric_tot_count_columns(color_rows)
		modCols = _numeric_mod_pct_columns(color_rows)
		mod_cols_by_target = {}
		for mod_col in modCols.columns:
			target = mod_col.replace('modPct.', '', 1).replace('modPct_', '', 1)
			mod_cols_by_target[target] = mod_col
		for cutoff in read_cutoffs:
			for tot_col in totCols.columns:
				totCol = totCols[tot_col]
				col_title = tot_col.replace('totCount.', '', 1).replace('totCount_', '', 1)
				if col_title not in mod_cols_by_target:
					continue
				modCol = modCols[mod_cols_by_target[col_title]]
				modAvg = modCol[totCol >= cutoff].mean()
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

	mod_pct_plot_obj_root = OutputContext(output_root).plot_root("modification_percentage_plot")
	plt.savefig(mod_pct_plot_obj_root+".pdf",pad_inches=1,bbox_inches='tight')
	plt.savefig(mod_pct_plot_obj_root+".png",pad_inches=1,bbox_inches='tight')
	logging.info("Finished modification percentage per amplicon plot.")
	summary_plot_obj = declared_plot_object(output_root, "modification_percentage_plot")
	return summary_plot_obj


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


def declared_plot_object(output_root, artifact_key):
	"""Build report metadata from a registered plot artifact."""
	metadata = OutputContext(output_root).plot_metadata(artifact_key)
	return PlotObject(
		plot_name=metadata["plot_name"],
		plot_title=metadata["plot_title"],
		plot_label=metadata["plot_label"],
		plot_datas=metadata["plot_datas"],
	)


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
