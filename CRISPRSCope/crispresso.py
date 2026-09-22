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
from CRISPRSCope.cache import (
	OutputRequirement,
	large_file_fingerprint,
	safe_remove_owned,
	tool_identity,
)
from CRISPRSCope.io_utils import open_text_maybe_gzip
from CRISPRSCope.output_artifacts import OutputContext, OutputManifest

from .paths import (
	STAGE_FILTER,
	STAGE_SPLIT,
	_command_to_string,
    _raise_command_error,
    build_stage_filename,
    safe_remove,
)
from .settings import _build_input_ref_names, _normalize_optional_guide
from .summaries import add_color_information, generate_amplicon_score

def _filter_to_hq_reads_at_single_amplicon(args):
	"""
	Filter paired per-amplicon FASTQ files (R1/R2) to only reads whose barcode
	is in `barcodes`.

	Parameters
	----------
	args : tuple
		(amp, amplicon_dir, barcodes, regenerate)

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

	amp, amplicon_dir, barcodes, regenerate = args
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
	if not regenerate and os.path.isfile(out1_name) and os.path.isfile(out2_name):
		return {
			"Amplicon": amp,
			"Status": "Success",
			"R1": out1_name,
			"R2": out2_name,
			"Reason": "MatchingBarcodeCache",
		}
	safe_remove(out1_name, silent=True)
	safe_remove(out2_name, silent=True)


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
		A 4-tuple ``(amp, amplicon_dir, barcodes, regenerate)``:
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
	amp_dict, amplicon_dir, barcodes, regenerate = args
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
	if not regenerate and os.path.isfile(amplicon_out_file):
		return {
			"Amplicon": amp,
			"Status": "Success",
			"In": amplicon_in_file,
			"Out": amplicon_out_file,
			"Reason": "MatchingBarcodeCache",
		}
	safe_remove(amplicon_out_file, silent=True)


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


def _barcode_set_sha256(barcodes):
	"""Return a deterministic digest for a selected barcode set."""
	digest = hashlib.sha256()
	for barcode in sorted(str(barcode) for barcode in barcodes):
		digest.update(barcode.encode('utf-8'))
		digest.update(b'\n')
	return digest.hexdigest()


def filter_amplicon_reads(output_root, parsed_information, amplicon_names,
							  cell_quality_to_analyze, n_processes, cache_manager=None):
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
	barcode_cache_file = output_root + ".filtered_barcodes.sha256"
	barcode_hash = _barcode_set_sha256(barcodes)
	cached_barcode_hash = None
	if os.path.isfile(barcode_cache_file):
		with open(barcode_cache_file, 'r') as handle:
			cached_barcode_hash = handle.read().strip()
	regenerate = cached_barcode_hash != barcode_hash

	logging.info("Filtering for reads from %d high quality barcodes...", len(barcodes))
	if regenerate:
		logging.info("Selected barcode set changed; regenerating filtered FASTQs")

	worker_args = [(amp, amplicon_dir, barcodes, regenerate) for amp in amplicon_names]

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
		safe_remove(barcode_cache_file, silent=True)
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

	allele_worker_args = [(amp, amplicon_dir, barcodes, regenerate) for amp in amp_success]

	with mp.Pool(n_processes) as pool:
		allele_results = pool.map(_filter_to_hq_alleles_at_single_amplicon, allele_worker_args)

	allele_failures = [r for r in allele_results if r.get('Status') != 'Success']
	if amp_failure or allele_failures:
		safe_remove(barcode_cache_file, silent=True)
	else:
		with open(barcode_cache_file, 'w') as handle:
			handle.write(barcode_hash + "\n")

	return {
		'read_filter_successes': amp_success,
		'read_filter_failures': amp_failure,
		'allele_filter_successes': [r for r in allele_results if r.get('Status') == 'Success'],
		'allele_filter_failures': allele_failures,
	}


def _decompressed_fastq_sha256(path):
	"""Return a stable SHA-256 digest of decompressed FASTQ content."""
	digest = hashlib.sha256()
	with open(path, 'rb') as raw_handle:
		is_gzip = raw_handle.read(2) == b'\x1f\x8b'
	open_fastq = gzip.open if is_gzip else open
	with open_fastq(path, 'rb') as handle:
		for chunk in iter(lambda: handle.read(1024 * 1024), b''):
			digest.update(chunk)
	return digest.hexdigest()


def _filtered_allele_fastq_path(output_root, amplicon_name):
	return build_stage_filename(
		stage=STAGE_FILTER,
		tag="alleles_qc_cells",
		amplicon=amplicon_name,
		ext="fq.gz",
		output_root=output_root + ".seq_by_amplicon",
	)


def _crispresso_cache_requirements(finished_file, crispresso_run_folder):
	return (
		OutputRequirement(
			"finished", finished_file, strategy="sha256", allow_empty=True,
		),
		OutputRequirement(
			"crispresso_info",
			os.path.join(crispresso_run_folder, "CRISPResso2_info.json"),
			strategy="sha256", validator="json",
		),
		OutputRequirement(
			"crispresso_fastq",
			os.path.join(crispresso_run_folder, "CRISPResso_output.fastq.gz"),
			strategy="stat",
		),
	)


def _build_crispresso_cache_record(
	cache_manager, amplicon_name, amplicon_info, suppress_sub_crispresso_plots,
	alleles, input_paths,
):
	stage = "crispresso_filtered" if alleles else "crispresso_reads"
	upstream_stage = "filter_selected" if alleles else "split_reads"
	upstream_scope = amplicon_name if alleles else "run"
	upstream = cache_manager.load(upstream_stage, upstream_scope)
	dependencies = [cache_manager.dependency(upstream)] if upstream is not None else []
	return cache_manager.new_record(
		stage,
		amplicon_name,
		algorithm_version=1,
		dependencies=dependencies,
		inputs={
			"fastqs": [large_file_fingerprint(path) for path in input_paths],
			"amplicon_sequences": amplicon_info["amp_seqs"],
			"guide_sequence": _normalize_optional_guide(amplicon_info.get("guide_seq", "")),
		},
		parameters={
			"mode": "filtered_alleles" if alleles else "reads",
			"suppress_sub_crispresso_plots": bool(suppress_sub_crispresso_plots),
		},
		tools={"CRISPResso": tool_identity(("CRISPResso", "--version"))},
	)


def run_crispresso_commands(amplicon_names,amplicon_information,output_root,crispresso_dir,suppress_sub_crispresso_plots,n_processes, alleles, cache_manager=None):
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
			os.makedirs(crispresso_dir, exist_ok=True)

	else:
		info_file = output_root+".crispresso.info.txt"


	cached_information = {}
	current_allele_input_hashes = {}
	if alleles and cache_manager is None:
		for amplicon_name in amplicon_names:
			allele_input = _filtered_allele_fastq_path(output_root, amplicon_name)
			if allele_input and os.path.isfile(allele_input):
				current_allele_input_hashes[amplicon_name] = _decompressed_fastq_sha256(allele_input)

	if os.path.isfile(info_file) and cache_manager is None:
		with open(info_file,'r') as fin:
			head = fin.readline().strip()
			head_els = head.split("\t")
			for line in fin:
				line_els = line.strip().split("\t")
				amp_info = dict(zip(head_els,line_els))
				cached_information[line_els[0]] = amp_info

		cache_is_valid = True
		for amplicon_name in amplicon_names:
			amp_info = cached_information.get(amplicon_name)
			if amp_info is None:
				cache_is_valid = False
				break
			if amp_info.get('status') == 'Failed':
				cache_is_valid = False
				logging.warning("Ignoring stale CRISPResso info cache because %s previously failed", amplicon_name)
				break
			if amp_info.get('status') != 'Completed':
				if alleles and current_allele_input_hashes.get(amplicon_name):
					cache_is_valid = False
					break
				continue
			finished_file = amp_info.get('finished_file')
			crispresso_run_folder = amp_info.get('crispresso_run_folder')
			crispresso_info_file = os.path.join(crispresso_run_folder, 'CRISPResso2_info.json')
			if not finished_file or not os.path.isfile(finished_file) or not os.path.isfile(crispresso_info_file):
				cache_is_valid = False
				logging.warning("Ignoring stale CRISPResso info cache because %s is missing completion outputs", amplicon_name)
				break
			if alleles:
				current_hash = current_allele_input_hashes.get(amplicon_name)
				if not current_hash or amp_info.get('input_sha256') != current_hash:
					cache_is_valid = False
					logging.info(
						"Invalidating filtered CRISPResso cache for %s because its allele input changed",
						amplicon_name,
					)
					break
		if cache_is_valid:
			logging.info ("Finished running CRISPResso on targets")
			return cached_information
		safe_remove(info_file, silent=True)

	crispresso_commands = []
	crispresso_information = {}
	cache_records = {}
	cache_requirements = {}

	not_run_count = 0
	finished_count = 0
	to_run_count = 0

	#print('Got to 2553')

	for amplicon_name in amplicon_names:
		crispresso_information[amplicon_name] = {}
		crispresso_information[amplicon_name]['name'] = amplicon_name
		crispresso_information[amplicon_name]['input_sha256'] = (
			current_allele_input_hashes.get(amplicon_name, '') if alleles else ''
		)
		if amplicon_information[amplicon_name]['aln_count'] == '0':
			#print(f"Got to 2559: {amplicon_name} within skipped block")
			crispresso_information[amplicon_name]['status'] = 'Skipped'
			crispresso_information[amplicon_name]['crispresso_command'] = 'NA'
			crispresso_information[amplicon_name]['crispresso_result'] = 'Skipped because had too few (%s) reads'%amplicon_information[amplicon_name]['aln_count']
			crispresso_information[amplicon_name]['finished_file'] = 'NA'
			crispresso_information[amplicon_name]['log_file'] = 'NA'
			crispresso_information[amplicon_name]['crispresso_run_folder'] = 'NA'
			not_run_count += 1
			if cache_manager is not None:
				stage = "crispresso_filtered" if alleles else "crispresso_reads"
				skip_record = cache_manager.new_record(
					stage, amplicon_name, algorithm_version=1,
					inputs={"aln_count": "0"},
					parameters={"mode": "filtered_alleles" if alleles else "reads"},
				)
				decision = cache_manager.evaluate(skip_record, ())
				if not decision.is_hit:
					cache_manager.commit(
						skip_record, (),
						result={"status": "skipped", "reason": "zero_aligned_reads"},
					)
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
			guide = _normalize_optional_guide(amplicon_information[amplicon_name].get('guide_seq', ''))


			# Separate crispresso_cmd for alleles and non alleles
			if alleles:
				# set pass allele_seq file through crispresso
				amp_filename = _filtered_allele_fastq_path(output_root, amplicon_name)
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

				crispresso_args = [
					"CRISPResso",
					"-r1", amp_filename,
					"-a", amplicon_seqs,
				]
				if guide:
					crispresso_args.extend(["-g", guide])

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
				crispresso_args = [
					"CRISPResso",
					"-r1", amp_filename_r1,
					"-r2", amp_filename_r2,
					"-a", amplicon_seqs,
				]
				if guide:
					crispresso_args.extend(["-g", guide])

			if suppress_sub_crispresso_plots:
				crispresso_args.extend(["--suppress_report", "--suppress_plots"])
				crispresso_args.extend([
				"-o", crispresso_dir,
				"-n", amplicon_name,
				"-w", "2",
				"--fastq_output",
					"--exclude_bp_from_left", "0",
				"--exclude_bp_from_right", "0",
			])
			if not alleles:
				crispresso_args.append("--crispresso_merge")
			crispresso_run_folder = os.path.join(crispresso_dir,'CRISPResso_on_'+amplicon_name)
			crispresso_cmd = shlex.join(crispresso_args) + " > " + shlex.quote(log_file) + " 2>&1 && touch " + shlex.quote(finished_file)
			crispresso_information[amplicon_name]['crispresso_command'] = crispresso_cmd
			crispresso_information[amplicon_name]['finished_file'] = finished_file
			crispresso_information[amplicon_name]['log_file'] = log_file
			crispresso_information[amplicon_name]['crispresso_run_folder'] = crispresso_run_folder

			crispresso_info_file = os.path.join(crispresso_run_folder, 'CRISPResso2_info.json')
			if cache_manager is not None:
				input_paths = [amp_filename] if alleles else [amp_filename_r1, amp_filename_r2]
				cache_record = _build_crispresso_cache_record(
					cache_manager, amplicon_name, amplicon_information[amplicon_name],
					suppress_sub_crispresso_plots, alleles, input_paths,
				)
				requirements = _crispresso_cache_requirements(finished_file, crispresso_run_folder)
				decision = cache_manager.evaluate(cache_record, requirements)
				if decision.is_hit:
					finished_count += 1
					crispresso_information[amplicon_name]['status'] = 'Completed'
					crispresso_information[amplicon_name]['crispresso_result'] = 'Completed'
					continue
				cache_records[amplicon_name] = cache_record
				cache_requirements[amplicon_name] = requirements
				safe_remove_owned(
					finished_file, allowed_root=crispresso_dir,
					expected_name=amplicon_name + ".finished",
				)
				safe_remove_owned(
					crispresso_run_folder, allowed_root=crispresso_dir,
					expected_name="CRISPResso_on_" + amplicon_name,
				)
			elif alleles and (os.path.isfile(finished_file) or os.path.isdir(crispresso_run_folder)):
				cached_hash = cached_information.get(amplicon_name, {}).get('input_sha256')
				current_hash = current_allele_input_hashes.get(amplicon_name)
				if not current_hash or cached_hash != current_hash:
					logging.info(
						"Removing stale filtered CRISPResso outputs for %s",
						amplicon_name,
					)
					safe_remove(finished_file, silent=True)
					if os.path.isdir(crispresso_run_folder):
						shutil.rmtree(crispresso_run_folder)
			if cache_manager is None and os.path.isfile(finished_file) and os.path.isfile(crispresso_info_file):
				finished_count += 1
				continue
			else:
				if os.path.isfile(finished_file) and not os.path.isfile(crispresso_info_file):
					os.remove(finished_file)
				crispresso_commands.append({
					'amplicon_name': amplicon_name,
					'args': crispresso_args,
					'command': crispresso_cmd,
					'log_file': log_file,
					'finished_file': finished_file,
					'crispresso_run_folder': crispresso_run_folder,
					'input_sha256': crispresso_information[amplicon_name]['input_sha256'],
				})


	logging.info('Skipped CRISPResso analysis for ' + str(not_run_count) + ' amplicons, finished analysis for ' + str(finished_count) + ' amplicons')

	logging.info('Got ' + str(len(crispresso_commands)) + ' CRISPResso commands')

	command_errors = []
	if len(crispresso_commands) > 0:
		# start processes
		logging.info("Running on "+ str(n_processes) + " processes..")
		pool = mp.Pool(n_processes)
		result = pool.map_async(run_crispresso_command, crispresso_commands).get(threading.TIMEOUT_MAX)
		pool.close()
		pool.join()
		for job, completed_job in zip(crispresso_commands, result):
			if completed_job.get('error'):
				command_errors.append(completed_job)
			elif cache_manager is not None:
				amplicon_name = job['amplicon_name']
				cache_manager.commit(
					cache_records[amplicon_name], cache_requirements[amplicon_name]
				)

	for amplicon_name in amplicon_names:
		if 'status' in crispresso_information[amplicon_name] and crispresso_information[amplicon_name]['status'] == 'Skipped':
			pass
		else:
			finished_file = crispresso_information[amplicon_name]['finished_file']
			crispresso_info_file = os.path.join(crispresso_information[amplicon_name]['crispresso_run_folder'], 'CRISPResso2_info.json')
			if os.path.isfile(finished_file) and os.path.isfile(crispresso_info_file):
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
				if os.path.isfile(finished_file) and not os.path.isfile(crispresso_info_file):
					error_message = 'Failed, missing CRISPResso2_info.json after CRISPResso exit; see ' + log_file
				crispresso_information[amplicon_name]['crispresso_result'] = error_message


	with open(info_file,'w') as fout:
		header_els = [
					'name',
					'input_sha256',
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
	if command_errors:
		completed_job = command_errors[0]
		_raise_command_error(
			completed_job.get('command'), completed_job.get('returncode'),
			context=completed_job.get('error'),
		)
	return crispresso_information


def run_crispresso_command(job):
	"""
	Run one CRISPResso command without a shell and mark it complete only
	after CRISPResso exits successfully and writes CRISPResso2_info.json.
	"""
	args = job['args']
	command = job['command']
	log_file = job['log_file']
	finished_file = job['finished_file']
	crispresso_run_folder = job['crispresso_run_folder']
	crispresso_info_file = os.path.join(crispresso_run_folder, 'CRISPResso2_info.json')
	try:
		logging.debug('running: ' + command)
		os.makedirs(os.path.dirname(log_file), exist_ok=True)
		with open(log_file, 'w') as log_handle:
			completed = sb.run(args, stdout=log_handle, stderr=sb.STDOUT, shell=False)
		if completed.returncode == 0 and os.path.isfile(crispresso_info_file):
			with open(finished_file, 'w'):
				pass
			return {'returncode': completed.returncode, 'error': None, 'command': command}
		if os.path.isfile(finished_file):
			os.remove(finished_file)
		error = 'CRISPResso exited with return code %s' % completed.returncode
		if completed.returncode == 0:
			error = 'CRISPResso exited successfully but did not write %s' % crispresso_info_file
		with open(log_file, 'a') as log_handle:
			log_handle.write('\nERROR: %s\n' % error)
		logging.error("%s on %s", error, command)
		return {'returncode': completed.returncode, 'error': error, 'command': command}
	except Exception as e:
		if os.path.isfile(finished_file):
			os.remove(finished_file)
		try:
			with open(log_file, 'a') as log_handle:
				log_handle.write('\nERROR: %s\n' % e)
		except Exception:
			pass
		logging.error("error: %s on %s" % (e, command))
		return {'returncode': None, 'error': str(e), 'command': command}


def run_command(cmd):
	"""
	Execute a shell command using subprocess.

	Parameters
	----------
	cmd : str
		Shell command to execute.

	Returns
	-------
	dict
		Return code and error message from subprocess call.

	Notes
	-----
	- Logs the command at DEBUG level.
	- Raises ExternalCommandError on failures.
	"""
	try:
		logging.debug('running: ' + _command_to_string(cmd))
		completed = sb.run(cmd, shell=isinstance(cmd, str), stdout=sb.PIPE, stderr=sb.PIPE)
	except Exception as e:
		logging.error("error: %s on %s" % (e, _command_to_string(cmd)))
		_raise_command_error(cmd, context=str(e))
	if completed.returncode != 0:
		stderr = completed.stderr.decode(errors="replace") if isinstance(completed.stderr, bytes) else completed.stderr
		logging.error("return code %s on %s", completed.returncode, _command_to_string(cmd))
		_raise_command_error(cmd, completed.returncode, stderr=stderr)
	return {'returncode': completed.returncode, 'error': None}


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
	- Uses subprocess.Popen without a shell when an argument list is provided.
	- stderr is redirected to stdout.
	- Caller is responsible for consuming the iterator.
	"""
	p = sb.Popen(command,
			stdout=sb.PIPE,
			stderr=sb.STDOUT,
			shell=isinstance(command, str),
			universal_newlines=True,
			bufsize=-1)#bufsize system default
	return iter(p.stdout.readline, '')


def parse_one_crispresso_output(this_args):
	"""
		Parse a single CRISPResso2 output folder into per-cell allele summaries.

	Reads the CRISPResso output FASTQ file, extracts:
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
	with open_text_maybe_gzip(crispresso_output_fastq,'rt') as fastq_input_handle:
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
	input_ref_names = _build_input_ref_names(num_references)

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
		fout.write("Ignore substitutions\t" + str(bool(ignore_substitutions))+"\n")
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


def _parse_cache_matches_ignore_substitutions(folder_finished_file, ignore_substitutions):
	if not os.path.isfile(folder_finished_file):
		return False

	expected_value = str(bool(ignore_substitutions))
	observed_value = None
	with open(folder_finished_file, 'r') as fin:
		for line in fin:
			if line.startswith("Ignore substitutions\t"):
				observed_value = line.rstrip("\n").split("\t", 1)[1]
				break

	if observed_value is None:
		return False
	return observed_value == expected_value


def parse_crispresso_outputs(amplicon_names,amplicon_information,amplicon_info_file,crispresso_information,
								output_root, min_total_reads_per_barcode, min_reads_per_amplicon_per_cell, n_processes,num_max_alleles=2,num_references=1,
								min_num_reads_per_cell=5,min_allele_pct_cutoff=.1,min_allele_count_cutoff=2,
								ignore_substitutions=False,
								write_alleles=False,
								amplicon_score_config=None,
								cache_manager=None):
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

			if not _parse_cache_matches_ignore_substitutions(folder_finished_file, ignore_substitutions):
				if os.path.isfile(folder_finished_file):
					logging.info(
						"Reparsing %s because ignore_substitutions changed to %s",
						name,
						ignore_substitutions,
					)
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

	usable_amplicon_names = []
	for amplicon_name in amplicon_names:
		if crispresso_information[amplicon_name]['status'] != 'Completed':
			continue
		crispresso_run_folder = crispresso_information[amplicon_name]['crispresso_run_folder']
		summ_file = crispresso_run_folder + ".summ"
		if os.path.isfile(summ_file):
			usable_amplicon_names.append(amplicon_name)

	logging.info(
		'Aggregating %d target summaries (%d with usable CRISPResso output)' % (
			len(amplicon_names),
			len(usable_amplicon_names),
		)
	)
	data = {}
	for amplicon_name in usable_amplicon_names:
		crispresso_run_folder = crispresso_information[amplicon_name]['crispresso_run_folder']
		summ_file = crispresso_run_folder + ".summ"
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

	outputs = OutputContext(output_root)
	with open(outputs.path("editing_summary_pseudobulk"),'w') as fout:
		header = "cell"
		for name in amplicon_names:
			header += "\ttotCount.%s\tmodPct.%s"%(name,name)
		fout.write(header+"\n")

		for cell in cells:
			line = cell
			for name in amplicon_names:
				val = "\tNA\tNA"
				if name in usable_amplicon_names:
					val = "\t0\tNA"
				if name in data[cell]:
					val = data[cell][name][0]
				line += val
			fout.write(line+"\n")

	with open(outputs.path("editing_summary"),'w') as fout:
		header = "cell"
		for name in amplicon_names:
			header += "\ttotCount.%s\tmodPct.%s"%(name,name)
		fout.write(header+"\n")

		for cell in cells:
			line = cell
			for name in amplicon_names:
				val = "\tNA\tNA"
				if name in usable_amplicon_names:
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
			val = ["NA", "NA"]
			if amplicon in usable_amplicon_names:
				val = [0, "NA"]
			if amplicon in data[cell]:
				val = data[cell][amplicon][0].strip().split('\t')
			val = [x if isinstance(x, int) else int(x) if x.isdigit() else x if x == "NA" else float(x) for x in val]
			vals += val
		summary_df.loc[cell] = vals


	usable_tot_cols = ["totCount.%s" % name for name in usable_amplicon_names]
	totCols = summary_df[usable_tot_cols].apply(pd.to_numeric, errors = 'coerce') if usable_tot_cols else pd.DataFrame(index = summary_df.index)

	# Always regenerate so scoring-code or setting changes cannot leave stale
	# classifications in downstream filtered outputs.
	amp_score_file = outputs.path("amplicon_score")
	amplicon_score_time = time.time()
	amp_score = generate_amplicon_score(
		totCols,
		min_reads_per_amplicon_per_cell=min_reads_per_amplicon_per_cell,
		min_total_reads_per_barcode=min_total_reads_per_barcode,
		config=amplicon_score_config,
	)
	amp_score.to_csv(amp_score_file, sep = "\t")
	end_amplicon_score_time = time.time() - amplicon_score_time
	logging.info("Generated supported-breadth amplicon score in %.2f seconds", end_amplicon_score_time)

	with open(outputs.path("filtered_editing_summary_pseudobulk"),'w') as fout:
		header = "cell"
		for name in amplicon_names:
			header += "\ttotCount.%s\tmodPct.%s"%(name,name)
		fout.write(header+"\n")

		for cell in cells:
			if cell not in amp_score.index:
				continue
			line = cell
			for name in amplicon_names:
				val = "\tNA\tNA"
				if name in usable_amplicon_names:
					val = "\t0\tNA"
				if name in data[cell]:
					val = data[cell][name][0]
				line += val
			fout.write(line+"\n")

	summary_df = add_color_information(summary_df, amp_score)

	return summary_df


def _crispresso_annotation_is_modified(annotation_line, ignore_substitutions=False):
	"""
	Classify one CRISPResso output FASTQ annotation line as modified/unmodified.
	"""
	if "ALN=NA" in annotation_line:
		return None
	field_values = {}
	for field in ("DEL", "INS", "SUB"):
		match = re.search(r"(?:^|\s)" + field + r"=([^\s]*)", annotation_line)
		if match is None:
			field_values[field] = None
		else:
			field_values[field] = match.group(1)
	if field_values["DEL"] is None or field_values["INS"] is None:
		return None
	if not ignore_substitutions and field_values["SUB"] is None:
		return None
	mod_fields = ["DEL", "INS"]
	if not ignore_substitutions:
		mod_fields.append("SUB")
	return any(field_values[field] not in (None, "") for field in mod_fields)


def _load_final_allele_read_support(summ_file):
	"""
	Read first-pass per-cell final allele read support from a .summ file.
	"""
	read_support = defaultdict(dict)
	if not summ_file or not os.path.isfile(summ_file):
		return read_support
	with open(summ_file, "r") as fin:
		header = fin.readline().strip().split("\t")
		header_idx = {name: idx for idx, name in enumerate(header)}
		if "cell" not in header_idx or "final_cell_allele_readcount_string" not in header_idx:
			return read_support
		for line in fin:
			line_els = line.rstrip("\n").split("\t")
			if len(line_els) <= header_idx["final_cell_allele_readcount_string"]:
				continue
			cell = line_els[header_idx["cell"]]
			readcount_string = line_els[header_idx["final_cell_allele_readcount_string"]]
			if readcount_string in ("", "NA"):
				continue
			for allele_idx, read_count in enumerate(readcount_string.split(","), start=1):
				try:
					read_support[cell][allele_idx] = int(read_count)
				except ValueError:
					read_support[cell][allele_idx] = 0
	return read_support


def _parse_filtered_crispresso_allele_output(crispresso_output_fastq, read_support, valid_barcodes, ignore_substitutions=False):
	"""
	Aggregate filtered CRISPResso allele classifications by barcode.
	"""
	results = defaultdict(lambda: {"support": 0, "modified": 0, "total": 0})
	if not crispresso_output_fastq or not os.path.isfile(crispresso_output_fastq):
		return results

	with open_text_maybe_gzip(crispresso_output_fastq, "rt") as fin:
		while True:
			header = fin.readline()
			if not header:
				break
			sequence = fin.readline()
			annotation = fin.readline()
			quality = fin.readline()
			if not quality:
				break

			header_token = header.strip().split(" ")[0]
			if header_token.startswith("@"):
				header_token = header_token[1:]
			header_els = header_token.split(":")
			if len(header_els) < 3:
				continue
			barcode = header_els[-2]
			try:
				allele_idx = int(header_els[-1])
			except ValueError:
				continue
			if barcode not in valid_barcodes:
				continue

			is_modified = _crispresso_annotation_is_modified(annotation, ignore_substitutions=ignore_substitutions)
			if is_modified is None:
				continue

			results[barcode]["total"] += 1
			if is_modified:
				results[barcode]["modified"] += 1
			results[barcode]["support"] += read_support.get(barcode, {}).get(allele_idx, 0)
	return results


def _write_filtered_summary_table(path, amplicon_names, cells, usable_amplicon_names, filtered_data):
	"""
	Write a filtered editing summary table with CRISPResso-derived genotype calls.
	"""
	with open(path, "w") as fout:
		header = "cell"
		for name in amplicon_names:
			header += "\ttotCount.%s\tmodPct.%s" % (name, name)
		fout.write(header + "\n")

		for cell in cells:
			line = cell
			for name in amplicon_names:
				val = "\tNA\tNA"
				if name in usable_amplicon_names:
					val = "\t0\tNA"
				if cell in filtered_data and name in filtered_data[cell]:
					this_data = filtered_data[cell][name]
					if this_data["total"] > 0:
						mod_pct = round(100 * this_data["modified"] / float(this_data["total"]), 2)
						val = "\t%s\t%s" % (this_data["support"], mod_pct)
				line += val
			fout.write(line + "\n")


def write_filtered_editing_summary_from_filtered_crispresso(
	amplicon_names,
	crispresso_information,
	crispresso_filtered_information,
	output_root,
	ignore_substitutions=False,
):
	"""
	Write filteredEditingSummary from filtered CRISPResso allele classifications.
	"""
	outputs = OutputContext(output_root)
	amp_score_file = outputs.path("amplicon_score")
	if not os.path.isfile(amp_score_file):
		raise FileNotFoundError("Amplicon score file does not exist: " + amp_score_file)
	amp_score = pd.read_csv(amp_score_file, sep="\t", index_col=0)
	cells = list(amp_score.index)
	valid_barcodes = set(cells)

	filtered_data = defaultdict(dict)
	usable_amplicon_names = []
	for amplicon_name in amplicon_names:
		filtered_info = crispresso_filtered_information.get(amplicon_name, {})
		if filtered_info.get("status") != "Completed":
			continue
		crispresso_run_folder = filtered_info.get("crispresso_run_folder")
		crispresso_output_fastq = os.path.join(crispresso_run_folder, "CRISPResso_output.fastq.gz") if crispresso_run_folder else None
		if not crispresso_output_fastq or not os.path.isfile(crispresso_output_fastq):
			continue

		usable_amplicon_names.append(amplicon_name)
		first_pass_info = crispresso_information.get(amplicon_name, {})
		first_pass_folder = first_pass_info.get("crispresso_run_folder")
		first_pass_summ = first_pass_folder + ".summ" if first_pass_folder else None
		read_support = _load_final_allele_read_support(first_pass_summ)
		amplicon_calls = _parse_filtered_crispresso_allele_output(
			crispresso_output_fastq,
			read_support,
			valid_barcodes,
			ignore_substitutions=ignore_substitutions,
		)
		for cell, call_data in amplicon_calls.items():
			filtered_data[cell][amplicon_name] = call_data

	_write_filtered_summary_table(
		outputs.path("filtered_editing_summary"),
		amplicon_names,
		cells,
		set(usable_amplicon_names),
		filtered_data,
	)
	logging.info("Finished writing filtered editing summary for %d filtered cells", len(cells))

	filtered_summary = pd.read_csv(outputs.path("filtered_editing_summary"), sep="\t", index_col=0)
	return add_color_information(filtered_summary, amp_score)
