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

from .settings import (
    AMPLICON_SCORE_MAX_BARCODE_RANK_DEFAULT,
    AMPLICON_SCORE_MIN_COVERED_FRACTION_DEFAULT,
    AMPLICON_SCORE_MIN_READS_PER_AMPLICON_DEFAULT,
    AmpliconScoreConfig,
)

def stratify_data(input_data, config=None):
	"""
	Assign quality category codes to barcodes.

	Barcodes are classified into one of:
		- 'HQ_HI'
		- 'HQ_LO'
		- 'LQ_HI'
		- 'LQ_LO'

	Classification is based on the configured supported-breadth fraction and
	barcode-rank cutoff.

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

	The input must contain ``Supported Amplicons``, ``Usable Amplicons``, and
	``Barcode Rank``. The returned frame is a copy.
	"""
	if config is None:
		config = AmpliconScoreConfig()
	required_columns = {'Supported Amplicons', 'Usable Amplicons', 'Barcode Rank'}
	missing_columns = sorted(required_columns - set(input_data.columns))
	if missing_columns:
		raise ValueError(f"Missing amplicon-score columns: {missing_columns}")

	result = input_data.copy()
	required_supported = np.ceil(
		config.min_covered_fraction * result['Usable Amplicons'].astype(float)
	).astype(int)
	high_score = result['Supported Amplicons'].astype(int) >= required_supported
	high_depth = result['Barcode Rank'].astype(int) <= config.max_barcode_rank

	result['Color'] = 'LQ_LO'
	result.loc[high_score & high_depth, 'Color'] = 'HQ_HI'
	result.loc[high_score & ~high_depth, 'Color'] = 'HQ_LO'
	result.loc[~high_score & high_depth, 'Color'] = 'LQ_HI'
	return result


def generate_amplicon_score(
	raw_tot_columns,
	min_reads_per_amplicon_per_cell,
	min_total_reads_per_barcode,
	config=None,
):
	"""
	Compute supported-breadth scores and assign quality categories to barcodes.

	Each usable amplicon contributes at most one unit of support. The score is
	the fraction of usable amplicons with at least the configured read count.

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

	The existing all-amplicon minimum is applied before scoring, and the total
	read-count filter is applied after deterministic ranking.
	"""
	if config is None:
		config = AmpliconScoreConfig()
	output_columns = [
		'Amplicon Score',
		'Supported Amplicons',
		'Usable Amplicons',
		'Read Count',
		'Barcode Rank',
		'Color',
	]

	if raw_tot_columns.empty or raw_tot_columns.shape[1] == 0:
		return pd.DataFrame(columns = output_columns)

	if raw_tot_columns.columns.duplicated().any():
		dup_cols = raw_tot_columns.columns[raw_tot_columns.columns.duplicated(keep = False)].unique()
		raise ValueError(f"Duplicate amplicon names detected: {list(dup_cols)}")

	raw_tot_columns = raw_tot_columns.apply(pd.to_numeric, errors = "coerce")
	raw_tot_columns = raw_tot_columns.dropna(axis = 1, how = "all")
	if raw_tot_columns.empty or raw_tot_columns.shape[1] == 0:
		return pd.DataFrame(columns = output_columns)

	# Preserve the existing optional gate across every usable amplicon.
	mask = (raw_tot_columns >= min_reads_per_amplicon_per_cell).all(axis = 1)
	raw_tot_columns = raw_tot_columns.loc[mask]

	logging.info('Cells that did not pass the read count per amplicon cutoff:' + str(len(mask) - sum(mask)))

	if raw_tot_columns.empty:
		return pd.DataFrame(columns = output_columns)

	usable_amplicon_count = raw_tot_columns.shape[1]
	supported_amplicons = raw_tot_columns.ge(config.min_reads_per_amplicon).sum(axis=1)
	amplicon_df = pd.DataFrame({
		'Amplicon Score': supported_amplicons / float(usable_amplicon_count),
		'Supported Amplicons': supported_amplicons.astype(int),
		'Usable Amplicons': usable_amplicon_count,
		'Read Count': raw_tot_columns.sum(axis = 1),
	})
	amplicon_df['_Barcode Sort'] = amplicon_df.index.astype(str)
	amplicon_df = amplicon_df.sort_values(
		['Read Count', '_Barcode Sort'],
		ascending=[False, True],
		kind='stable',
	).drop(columns=['_Barcode Sort'])
	amplicon_df['Barcode Rank'] = range(1, len(amplicon_df) + 1)

	amplicon_df = stratify_data(amplicon_df, config=config)

	filtered_df = amplicon_df[amplicon_df['Read Count'] >= min_total_reads_per_barcode]
	return filtered_df[output_columns]


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
