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
from CRISPRSCope.cache import CacheConfig, tool_identity
from CRISPRSCope.io_utils import open_text_maybe_gzip
from CRISPRSCope.output_artifacts import OutputContext, OutputManifest
from .allele_calling import parse_allele_support

MIN_TOTAL_READS_PER_BARCODE_DEFAULT = 10
MIN_READS_PER_AMPLICON_PER_CELL_DEFAULT = 0
PARTIAL_RESCUE_MIN_MEAN_READ_QUALITY_DEFAULT = 30.0
AMPLICON_SCORE_MIN_READS_PER_AMPLICON_DEFAULT = 5
AMPLICON_SCORE_MIN_COVERED_FRACTION_DEFAULT = 2.0 / 3.0
AMPLICON_SCORE_MAX_BARCODE_RANK_DEFAULT = 10_000
H5AD_ZYGOSITY_DEFAULTS = {
    "wt_max_mod_pct": 20.0,
    "het_max_mod_pct": 80.0,
    "hom_min_mod_pct": 80.0,
    "compound_het_min_allele2_pct": 20.0,
}
CELL_QUALITY_CODES = {
    "HQ_HI": {"label": "HighScore_HighDepth", "display": "High score / High depth", "color": "#1f77b4"},
    "HQ_LO": {"label": "HighScore_LowDepth", "display": "High score / Low depth", "color": "#2ca02c"},
    "LQ_HI": {"label": "LowScore_HighDepth", "display": "Low score / High depth", "color": "#ff7f0e"},
    "LQ_LO": {"label": "LowScore_LowDepth", "display": "Low score / Low depth", "color": "#d62728"},
}

@dataclass(frozen=True)
class AmpliconScoreConfig:
	"""Settings for supported-breadth barcode classification."""

	min_reads_per_amplicon: int = AMPLICON_SCORE_MIN_READS_PER_AMPLICON_DEFAULT
	min_covered_fraction: float = AMPLICON_SCORE_MIN_COVERED_FRACTION_DEFAULT
	max_barcode_rank: int = AMPLICON_SCORE_MAX_BARCODE_RANK_DEFAULT

	def __post_init__(self):
		if self.min_reads_per_amplicon < 1:
			raise ValueError("amplicon_score_min_reads_per_amplicon must be >= 1")
		if not np.isfinite(self.min_covered_fraction) or not 0 < self.min_covered_fraction <= 1:
			raise ValueError("amplicon_score_min_covered_fraction must be > 0 and <= 1")
		if self.max_barcode_rank < 1:
			raise ValueError("amplicon_score_max_barcode_rank must be >= 1")


def _parse_cache_config(settings_file):
	"""Parse the cache policy without changing the legacy settings tuple."""
	settings = _parse_settings_file(settings_file)
	return CacheConfig.from_value(settings.get("cache_mode", "auto"))


def _resolve_settings_path(path_value: str, settings_dir: str) -> str:
	"""
	Resolve a path from the settings file relative to the settings file directory.
	"""
	path_value = str(path_value).strip()
	if os.path.isabs(path_value):
		return os.path.abspath(path_value)
	return os.path.abspath(os.path.join(settings_dir, path_value))


def _resolve_settings_path_list(path_values: str, settings_dir: str) -> str:
	"""
	Resolve a comma-separated path list from the settings file.
	"""
	return ",".join(
		_resolve_settings_path(path_value, settings_dir)
		for path_value in str(path_values).split(",")
	)


def _parse_bool_setting(settings, key, default=False):
	if key not in settings:
		return default
	val = str(settings[key]).strip().lower()
	if val in ("true", "yes", "1"):
		return True
	if val in ("false", "no", "0"):
		return False
	raise ValueError(f"Invalid value for {key}: {settings[key]!r}. Expected True or False.")


def _parse_int_setting(settings, key, default, minimum=None):
	raw_value = settings.get(key, default)
	try:
		value = int(raw_value)
	except Exception as e:
		raise ValueError(f"Invalid value for {key}: {raw_value!r} ({e})")
	if minimum is not None and value < minimum:
		raise ValueError(f"{key} must be >= {minimum}")
	return value


def _parse_float_setting(settings, key, default, minimum=None):
	raw_value = settings.get(key, default)
	try:
		value = float(raw_value)
	except Exception as e:
		raise ValueError(f"Invalid value for {key}: {raw_value!r} ({e})")
	if minimum is not None and value < minimum:
		raise ValueError(f"{key} must be >= {minimum}")
	return value


def _parse_amplicon_score_config(settings_file):
	"""Parse and validate supported-breadth amplicon-score settings."""
	settings = _parse_settings_file(settings_file)
	min_reads_per_amplicon = _parse_int_setting(
		settings,
		'amplicon_score_min_reads_per_amplicon',
		AMPLICON_SCORE_MIN_READS_PER_AMPLICON_DEFAULT,
		minimum=1,
	)
	min_covered_fraction = _parse_float_setting(
		settings,
		'amplicon_score_min_covered_fraction',
		AMPLICON_SCORE_MIN_COVERED_FRACTION_DEFAULT,
	)
	if not np.isfinite(min_covered_fraction) or not 0 < min_covered_fraction <= 1:
		raise ValueError("amplicon_score_min_covered_fraction must be > 0 and <= 1")
	max_barcode_rank = _parse_int_setting(
		settings,
		'amplicon_score_max_barcode_rank',
		AMPLICON_SCORE_MAX_BARCODE_RANK_DEFAULT,
		minimum=1,
	)
	return AmpliconScoreConfig(
		min_reads_per_amplicon=min_reads_per_amplicon,
		min_covered_fraction=min_covered_fraction,
		max_barcode_rank=max_barcode_rank,
	)


def _parse_allele_calling_config(settings_file):
	"""Return the genotype depth and the typed allele support threshold."""
	settings = _parse_settings_file(settings_file)
	for legacy in ('min_allele_pct_cutoff', 'min_allele_count_cutoff'):
		if legacy in settings:
			raise ValueError(f"{legacy} was unused and is retired; use min_allele_support instead")
	depth = _parse_int_setting(settings, 'min_reads_per_amplicon_for_genotype', 8, minimum=0)
	support = parse_allele_support(settings.get('min_allele_support', '2'))
	return depth, support.value


def _parse_editing_rate_ci_config(settings_file):
	"""Parse and validate editing-rate confidence interval settings."""
	from CRISPRSCope.editing_rate_ci import EditingRateCIConfig

	settings = _parse_settings_file(settings_file)
	enabled = _parse_bool_setting(settings, 'write_editing_rate_ci', default=True)
	bootstrap_iterations = _parse_int_setting(
		settings,
		'editing_rate_ci_bootstrap_iterations',
		10_000,
		minimum=100,
	)
	permutation_iterations = _parse_int_setting(
		settings,
		'editing_rate_ci_permutation_iterations',
		10_000,
		minimum=100,
	)
	confidence_level = _parse_float_setting(
		settings,
		'editing_rate_ci_confidence_level',
		0.95,
	)
	if not 0 < confidence_level < 1:
		raise ValueError("editing_rate_ci_confidence_level must be greater than 0 and less than 1")
	seed = _parse_int_setting(settings, 'editing_rate_ci_seed', 42, minimum=0)
	coverage_exact_max_reads = _parse_int_setting(
		settings,
		'editing_rate_ci_coverage_exact_max_reads',
		10,
		minimum=0,
	)
	coverage_bin_width_reads = _parse_int_setting(
		settings,
		'editing_rate_ci_coverage_bin_width_reads',
		5,
		minimum=1,
	)
	return EditingRateCIConfig(
		enabled=enabled,
		bootstrap_iterations=bootstrap_iterations,
		permutation_iterations=permutation_iterations,
		confidence_level=confidence_level,
		seed=seed,
		coverage_exact_max_reads=coverage_exact_max_reads,
		coverage_bin_width_reads=coverage_bin_width_reads,
	)


def _parse_editing_rate_depth_stability_config(settings_file):
	"""Parse and validate editing-rate depth-stability settings."""
	from CRISPRSCope.editing_rate_ci import EditingRateDepthStabilityConfig

	settings = _parse_settings_file(settings_file)
	enabled = _parse_bool_setting(
		settings,
		'write_editing_rate_depth_stability',
		default=False,
	)
	iterations = _parse_int_setting(
		settings,
		'editing_rate_depth_stability_iterations',
		1_000,
		minimum=100,
	)
	raw_percentages = settings.get(
		'editing_rate_depth_stability_percentages',
		'10,25,50,75,90',
	)
	try:
		percentage_parts = [part.strip() for part in str(raw_percentages).split(',')]
		if not percentage_parts or any(not part for part in percentage_parts):
			raise ValueError("expected a non-empty comma-separated list")
		percentages = tuple(float(part) for part in percentage_parts)
	except Exception as e:
		raise ValueError(
			"Invalid value for editing_rate_depth_stability_percentages: "
			f"{raw_percentages!r} ({e})"
		)
	if any(not np.isfinite(value) or not 0 < value < 100 for value in percentages):
		raise ValueError(
			"editing_rate_depth_stability_percentages values must be greater than 0 and less than 100"
		)
	if any(left >= right for left, right in zip(percentages, percentages[1:])):
		raise ValueError(
			"editing_rate_depth_stability_percentages must be strictly increasing and unique"
		)
	confidence_level = _parse_float_setting(
		settings,
		'editing_rate_ci_confidence_level',
		0.95,
	)
	if not 0 < confidence_level < 1:
		raise ValueError("editing_rate_ci_confidence_level must be greater than 0 and less than 1")
	seed = _parse_int_setting(settings, 'editing_rate_ci_seed', 42, minimum=0)
	return EditingRateDepthStabilityConfig(
		enabled=enabled,
		iterations=iterations,
		percentages=percentages,
		confidence_level=confidence_level,
		seed=seed,
	)


def _settings_value_is_path(value: str) -> bool:
	return str(value).strip().lower() not in ("", "true", "yes", "1", "false", "no", "0", "none")


def _normalize_optional_guide(guide):
	if guide is None:
		return ""
	guide = str(guide).strip()
	if guide.lower() in ("", "na", "none"):
		return ""
	return guide


def _build_input_ref_names(num_references):
	return ['Reference'] + ['Amplicon' + str(i) for i in range(1, num_references)]


def _resolve_existing_fastq_path(path):
	if not path or path in ("NA", "None"):
		return None
	if os.path.isfile(path):
		return path
	if not str(path).endswith(".gz") and os.path.isfile(path + ".gz"):
		return path + ".gz"
	return None


def _parse_settings_file(settings_file):
	settings = {}
	with open(settings_file, 'r') as sin:
		for line_number, line in enumerate(sin, start=1):
			line = line.rstrip("\n")
			stripped = line.strip()
			if stripped == "" or stripped.startswith("#"):
				continue
			if "\t" not in line:
				raise ValueError(
					f"Invalid settings line {line_number}: line must be key<TAB>value"
				)
			key, value = line.split("\t", 1)
			key = key.strip()
			value = value.strip()
			if not key or value == "":
				raise ValueError(
					f"Invalid settings line {line_number}: line must be key<TAB>value"
				)
			settings[key] = value
	return settings


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
			debug_rescued_reads_bam : str,
			debug_rejected_rescue_reads_bam : str,
			debug_require_strict_amplicon_alignment : bool,
			partial_rescue_min_mean_read_quality : float,
			write_output_manifest : bool,
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
	settings_dir = os.path.dirname(settings_file)

	# Checking if a debug parameter has been provided to the function call #
	logging_level = logging.INFO
	if len(args) > 2 and 'debug' in args[2].lower():
		logging_level=logging.DEBUG

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
	logging.getLogger('matplotlib').setLevel(logging.WARNING)
	logging.getLogger('fontTools').setLevel(logging.WARNING)
	logging.getLogger('fontTools.subset').setLevel(logging.WARNING)

	logging.info('Parsing settings file..')

	# Parse the settings file and write to a dictionary {key '\t' value}
	settings = _parse_settings_file(settings_file)
	_parse_allele_calling_config(settings_file)

	# Checking for various required settings and raising exceptions for missing values
	if 'r1' not in settings:
		raise Exception('Settings file must contain an entry for r1')
	r1 = _resolve_settings_path_list(settings['r1'], settings_dir)
	# Verifying that all of the read files exist
	for r1_file in r1.split(","):
		if not os.path.isfile(r1_file):
			raise Exception('Input r1 ' + r1_file + ' does not exist')

	# Same as above for r1
	if 'r2' not in settings:
		raise Exception('Settings file must contain an entry for r2')
	r2 = _resolve_settings_path_list(settings['r2'], settings_dir)
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

	allow_barcode_mismatches = _parse_bool_setting(settings, 'allowBarcodeMismatches', default=False)

	# Checking for inclusion and existence of a barcode file
	if 'barcodes' not in settings:
		raise Exception('Settings file must contain an entry for barcodes')
	barcode_file = _resolve_settings_path(settings['barcodes'], settings_dir)
	if not os.path.isfile(barcode_file):
		raise Exception('Barcode file ' + barcode_file + ' does not exist')


	# Checking for existence and inclusion of an amplicon file
	if 'amplicons' not in settings:
		raise Exception('Settings file must contain an entry for amplicons')
	amplicon_file = _resolve_settings_path(settings['amplicons'], settings_dir)
	if not os.path.isfile(amplicon_file):
		raise Exception('Amplicon file does not exist at ' + amplicon_file)

	# Settings primer lookup length
	if 'primerLookupLen' in settings:
		raise ValueError("primerLookupLen is no longer supported. Use primer_lookup_len instead.")
	primer_lookup_len = _parse_int_setting(settings, 'primer_lookup_len', 18, minimum=1)

	# Settings adapter DNA seq
	adapter_DNA = "TGTCTCTTATACACATCTCCGAGCCCACGAG"
	if 'adapter_DNA' in settings:
		adapter_DNA = settings['adapter_DNA'].strip()
	else:
		for legacy_key in ('plamsid_DNA', 'plasmid_DNA'):
			if legacy_key in settings:
				adapter_DNA = settings[legacy_key].strip()
				logging.warning("%s is deprecated; use adapter_DNA instead.", legacy_key)
				break

	n_processes = mp.cpu_count()
	if 'processes' in settings and settings['processes'] != 'max':
		n_processes=int(settings['processes'])
	if not n_processes or n_processes < 1:
		raise ValueError("n_processes must be >= 1")

	keep_intermediate_files = _parse_bool_setting(settings, 'keep_intermediate_files', default=False)

	ignore_substitutions = _parse_bool_setting(settings, 'ignore_substitutions', default=False)

	assign_reads_to_all_possible_amplicons = _parse_bool_setting(settings, 'assign_reads_to_all_possible_amplicons', default=False)

	debug_require_strict_amplicon_alignment = _parse_bool_setting(settings, 'debug_require_strict_amplicon_alignment', default=False)
	if assign_reads_to_all_possible_amplicons and debug_require_strict_amplicon_alignment:
		raise ValueError("debug_require_strict_amplicon_alignment cannot be used with assign_reads_to_all_possible_amplicons")

	suppress_sub_crispresso_plots = _parse_bool_setting(settings, 'suppress_sub_crispresso_plots', default=False)

	write_h5ad = _parse_bool_setting(settings, 'write_h5ad', default=True)
	write_output_manifest = _parse_bool_setting(settings, 'write_output_manifest', default=False)

	# --- normalized cutoffs (use canonical defaults and validate) ---
	min_total_reads_per_barcode = _parse_int_setting(settings, 'min_total_reads_per_barcode', MIN_TOTAL_READS_PER_BARCODE_DEFAULT, minimum=0)
	min_reads_per_amplicon_per_cell = _parse_int_setting(settings, 'min_reads_per_amplicon_per_cell', MIN_READS_PER_AMPLICON_PER_CELL_DEFAULT, minimum=0)

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
			if _parse_bool_setting(settings, key):
				cell_quality_to_analyze.add(code)

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
		bowtie2_index = _resolve_settings_path(settings['bowtie2_index'].replace(".fa",""), settings_dir)
	if 'genome' in settings:
		bowtie2_index = _resolve_settings_path(settings['genome'].replace(".fa",""), settings_dir)
	if not os.path.isfile(bowtie2_index+".1.bt2") and not os.path.isfile(bowtie2_index+".1.bt2l"):
		raise Exception('bowtie2_index file does not exist at ' + bowtie2_index + ".bt2 or " + bowtie2_index + ".bt2l")

	# Checking for existence and inclusion of an alternate alleles file
	alt_alleles_file = ""
	if 'alt_alleles_file' in settings:
		alt_alleles_file = _resolve_settings_path(settings['alt_alleles_file'], settings_dir)
		if not os.path.isfile(alt_alleles_file):
			raise Exception('Alt alleles file does not exist at ' + alt_alleles_file)

	# Checking for a user provided output_root
		# if not provided, an suffix is appended onto the settings file name
	output_root = settings_file
	if 'output_root' in settings:
		output_root = _resolve_settings_path(settings['output_root'], settings_dir)

	h5ad_output = output_root + ".h5ad"
	if 'h5ad_output' in settings:
		h5ad_output = _resolve_settings_path(settings['h5ad_output'], settings_dir)

	debug_rescued_reads_bam = ""
	if 'debug_rescued_reads_bam' in settings:
		debug_rescued_reads_bam = settings['debug_rescued_reads_bam'].strip()
		if debug_rescued_reads_bam.lower() in ("true", "yes", "1"):
			debug_rescued_reads_bam = output_root + ".splitReads.rescued.bam"
		elif debug_rescued_reads_bam.lower() in ("false", "no", "0", "none"):
			debug_rescued_reads_bam = ""
		elif _settings_value_is_path(debug_rescued_reads_bam):
			debug_rescued_reads_bam = _resolve_settings_path(debug_rescued_reads_bam, settings_dir)

	debug_rejected_rescue_reads_bam = ""
	if 'debug_rejected_rescue_reads_bam' in settings:
		debug_rejected_rescue_reads_bam = settings['debug_rejected_rescue_reads_bam'].strip()
		if debug_rejected_rescue_reads_bam.lower() in ("true", "yes", "1"):
			debug_rejected_rescue_reads_bam = output_root + ".splitReads.rejected_rescue_candidates.bam"
		elif debug_rejected_rescue_reads_bam.lower() in ("false", "no", "0", "none"):
			debug_rejected_rescue_reads_bam = ""
		elif _settings_value_is_path(debug_rejected_rescue_reads_bam):
			debug_rejected_rescue_reads_bam = _resolve_settings_path(debug_rejected_rescue_reads_bam, settings_dir)

	partial_rescue_min_mean_read_quality = _parse_float_setting(
		settings,
		'partial_rescue_min_mean_read_quality',
		PARTIAL_RESCUE_MIN_MEAN_READ_QUALITY_DEFAULT,
		minimum=0,
	)

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
		os.makedirs(amp_file_dir, exist_ok=True)

	crispresso_dir = output_root + ".crispresso"
	if not os.path.isdir(crispresso_dir):
		os.makedirs(crispresso_dir, exist_ok=True)


	# Resolve required tools once. Cache keys reuse these memoized identities.
	for command, error_message in (
		(("bowtie2", "--version"), "Error: bowtie2 is required"),
		(("samtools", "--version"), "Error: samtools is required"),
		(("CRISPResso", "--version"), "Error: CRISPResso2 is required"),
	):
		try:
			tool_identity(command)
		except Exception as error:
			raise Exception(error_message) from error


	return (r1, r2, constant1, constant2, allow_barcode_mismatches,barcode_file, amplicon_file, primer_lookup_len, adapter_DNA, amp_file_dir, alt_alleles_file, bowtie2_index, crispresso_dir, output_root, n_processes, keep_intermediate_files, ignore_substitutions, assign_reads_to_all_possible_amplicons, suppress_sub_crispresso_plots, min_total_reads_per_barcode, min_reads_per_amplicon_per_cell, cell_quality_to_analyze, write_h5ad, h5ad_output, h5ad_export_config, debug_rescued_reads_bam, debug_rejected_rescue_reads_bam, debug_require_strict_amplicon_alignment, partial_rescue_min_mean_read_quality, write_output_manifest, settings_file)
