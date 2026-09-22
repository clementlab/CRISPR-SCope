"""Shared imports for CRISPRSCope pipeline-stage modules."""
from __future__ import annotations

import errno
import glob
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
	optional_file_fingerprint,
	safe_remove_owned,
	small_file_fingerprint,
	tool_identity,
)
from CRISPRSCope.io_utils import open_text_maybe_gzip
from CRISPRSCope.output_artifacts import OutputContext, OutputManifest

from .crispresso import get_command_output, run_command
from .fastq_processing import reverse_complement
from .paths import (
	STAGE_SPLIT,
	_command_to_string,
	_raise_command_error,
	build_stage_filename,
	safe_remove,
	safe_write_path,
)
from .settings import PARTIAL_RESCUE_MIN_MEAN_READ_QUALITY_DEFAULT, _resolve_existing_fastq_path


def _split_cache_requirements(
	output_root, amp_file_dir, info_file, amplicon_information=None,
	debug_rescued_reads_bam="", debug_rejected_rescue_reads_bam="",
):
	outputs = OutputContext(output_root)
	requirements = [
		OutputRequirement("amplicon_info", info_file, strategy="sha256", validator="tsv", required_header=("name",)),
		OutputRequirement("valid_amplicons", outputs.path("valid_amplicons"), strategy="sha256", allow_empty=True),
		OutputRequirement("aligned_read_counts", outputs.path("aligned_read_counts"), strategy="sha256", validator="tsv", required_header=("Barcode", "Aligned Count")),
		OutputRequirement("unaligned_read_counts", outputs.path("unaligned_read_counts"), strategy="sha256", validator="tsv", required_header=("Barcode", "Unaligned Count")),
		OutputRequirement("amplicon_classification", outputs.path("amplicon_classification"), strategy="sha256", validator="tsv", required_header=("is_valid", "amp1_from_seq")),
	]
	for amp_name, amp_info in (amplicon_information or {}).items():
		if str(amp_info.get("aln_count", "0")) == "0":
			continue
		for read in ("r1", "r2"):
			path = build_stage_filename(
				stage=STAGE_SPLIT, tag="reads_all_cells", amplicon=amp_name,
				read=read, ext="fq.gz", output_root=amp_file_dir,
			)
			requirements.append(OutputRequirement(
				f"reads:{amp_name}:{read}", path, strategy="stat", allow_empty=True,
				validator="gzip",
			))
	if debug_rescued_reads_bam:
		requirements.append(OutputRequirement(
			"debug_rescued_reads", debug_rescued_reads_bam, strategy="stat",
			allow_empty=True, validator="bam",
		))
	if debug_rejected_rescue_reads_bam:
		requirements.append(OutputRequirement(
			"debug_rejected_rescue_reads", debug_rejected_rescue_reads_bam,
			strategy="stat", allow_empty=True, validator="bam",
		))
	return tuple(requirements)


def _prune_obsolete_split_outputs(
	cache_manager, previous_record, current_requirements, amp_file_dir,
):
	"""Remove old per-amplicon FASTQs only when a trusted split record owned them."""
	if previous_record is None or cache_manager.config.mode.value == "disabled":
		return []
	current_paths = {
		requirement.normalized_path() for requirement in current_requirements
	}
	removed = []
	for output in previous_record.outputs:
		key = str(output.get("key", ""))
		if not key.startswith("reads:"):
			continue
		try:
			amplicon_name, read = key[len("reads:"):].rsplit(":", 1)
		except ValueError:
			logging.warning("Ignoring malformed split cache output key during cleanup: %s", key)
			continue
		if read not in {"r1", "r2"}:
			logging.warning("Ignoring unexpected split cache read key during cleanup: %s", key)
			continue
		expected = build_stage_filename(
			stage=STAGE_SPLIT, tag="reads_all_cells", amplicon=amplicon_name,
			read=read, ext="fq.gz", output_root=amp_file_dir,
		)
		if os.path.abspath(str(output.get("path", ""))) != os.path.abspath(expected):
			logging.warning("Ignoring split cache output with unexpected path during cleanup: %s", key)
			continue
		if os.path.abspath(expected) in current_paths:
			continue
		try:
			if safe_remove_owned(
				expected, allowed_root=amp_file_dir,
				expected_name=os.path.basename(expected),
			):
				removed.append(expected)
		except ValueError as error:
			logging.warning("Refusing unsafe stale split cleanup for %s: %s", key, error)
	return removed


def _build_split_cache_record(
	cache_manager, aligned_bam, amplicon_file, alt_alleles_file, bowtie2_index,
	primer_lookup_len, adapter_DNA, min_total_reads_per_barcode,
	assign_reads_to_all_possible_amplicons, debug_rescued_reads_bam,
	debug_require_strict_amplicon_alignment, debug_rejected_rescue_reads_bam,
	partial_rescue_min_mean_read_quality,
):
	parse_record = cache_manager.load("parse_align")
	dependencies = []
	if parse_record is not None:
		dependencies.append(cache_manager.dependency(parse_record, ("aligned_bam",)))
	index_files = sorted(
		path for path in glob.glob(bowtie2_index + ".*")
		if path.endswith((".bt2", ".bt2l"))
	)
	return cache_manager.new_record(
		"split_reads",
		algorithm_version=1,
		dependencies=dependencies,
		inputs={
			"aligned_bam": large_file_fingerprint(aligned_bam),
			"amplicons": small_file_fingerprint(amplicon_file),
			"alternate_alleles": optional_file_fingerprint(alt_alleles_file, strategy="sha256"),
			"bowtie2_index": [large_file_fingerprint(path) for path in index_files],
		},
		parameters={
			"primer_lookup_len": int(primer_lookup_len),
			"adapter_DNA": adapter_DNA,
			"min_total_reads_per_barcode": int(min_total_reads_per_barcode),
			"assign_reads_to_all_possible_amplicons": bool(assign_reads_to_all_possible_amplicons),
			"debug_rescued_reads_bam": os.path.abspath(debug_rescued_reads_bam) if debug_rescued_reads_bam else "",
			"debug_require_strict_amplicon_alignment": bool(debug_require_strict_amplicon_alignment),
			"debug_rejected_rescue_reads_bam": os.path.abspath(debug_rejected_rescue_reads_bam) if debug_rejected_rescue_reads_bam else "",
			"partial_rescue_min_mean_read_quality": float(partial_rescue_min_mean_read_quality),
		},
		tools={
			"bowtie2": tool_identity(("bowtie2", "--version")),
			"samtools": tool_identity(("samtools", "--version")),
		},
	)

def _load_split_read_cache(info_file, amp_file_dir):
	amplicon_names = []
	amplicon_information = {}
	cache_is_valid = True
	with open(info_file, 'r') as fin:
		head = fin.readline().strip()
		head_els = head.split("\t")
		for line in fin:
			line_els = line.strip().split("\t")
			if not line_els or line_els[0] == "":
				continue
			amp_info = dict(zip(head_els,line_els))
			amp_name = line_els[0]
			if 'reads_r1_file' not in amp_info or amp_info['reads_r1_file'] in ('', 'NA'):
				amp_info['reads_r1_file'] = build_stage_filename(
					stage=STAGE_SPLIT,
					tag="reads_all_cells",
					amplicon=amp_name,
					read="r1",
					ext="fq",
					output_root=amp_file_dir,
				)
			if 'reads_r2_file' not in amp_info or amp_info['reads_r2_file'] in ('', 'NA'):
				amp_info['reads_r2_file'] = build_stage_filename(
					stage=STAGE_SPLIT,
					tag="reads_all_cells",
					amplicon=amp_name,
					read="r2",
					ext="fq",
					output_root=amp_file_dir,
				)
			if amp_info.get('aln_count') not in ('0', 0, None):
				r1_cached = _resolve_existing_fastq_path(amp_info.get('reads_r1_file'))
				r2_cached = _resolve_existing_fastq_path(amp_info.get('reads_r2_file'))
				if not r1_cached or not r2_cached:
					cache_is_valid = False
					logging.warning(
						"Ignoring stale split-read cache for %s because referenced FASTQ files are missing",
						amp_name,
					)
				else:
					amp_info['reads_r1_file'] = r1_cached
					amp_info['reads_r2_file'] = r2_cached
			amplicon_information[amp_name] = amp_info
			amplicon_names.append(amp_name)
	return cache_is_valid, amplicon_names, amplicon_information


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


def alignment_end(pos, cigar):
	if cigar == "*" or cigar is None:
		return pos
	ref_len = 0
	for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
		length = int(length)
		if op in ('M', 'D', 'N', '=', 'X'):
			ref_len += length
	return pos + ref_len - 1


def mean_phred_quality(qual_string):
	"""
	Return the mean Phred+33 base quality for one SAM quality string.
	"""
	if qual_string is None:
		return None
	qual_string = qual_string.strip()
	if qual_string == "" or qual_string == "*":
		return None
	return sum(ord(ch) - 33 for ch in qual_string) / len(qual_string)


def alignment_boundary_key(read_chr, read_start0, cigar, is_reverse):
	"""
	Return the inward-facing amplicon-boundary key supported by an alignment.

	Forward-strand reads support the amplicon start boundary. Reverse-strand
	reads support the amplicon end boundary. Outward-facing boundary hits are
	therefore ignored by construction.
	"""
	if read_chr in (None, "", "*") or cigar in (None, "*") or read_start0 is None:
		return None
	read_pos = alignment_end(read_start0, cigar) if is_reverse else read_start0
	return read_chr + ":" + str(read_pos)


def inward_alignment_amplicon(read_chr, read_start0, cigar, is_reverse, start_lookup, end_lookup):
	"""
	Call an amplicon from an inward-facing alignment boundary, or NA.
	"""
	key = alignment_boundary_key(read_chr, read_start0, cigar, is_reverse)
	if key is None:
		return "NA"
	if is_reverse:
		return end_lookup.get(key, "NA")
	return start_lookup.get(key, "NA")


def _open_rescued_reads_writer(aligned_bam, rescued_reads_path):
	"""
	Open a SAM/BAM writer for read pairs accepted by partial alignment rescue.
	"""
	if not rescued_reads_path:
		return (None, None)

	safe_write_path(rescued_reads_path)
	header = sb.check_output(
		["samtools", "view", "-H", aligned_bam],
		universal_newlines=True,
	)

	if rescued_reads_path.lower().endswith(".sam"):
		writer = open(rescued_reads_path, "w")
		writer.write(header)
		return (writer, None)

	proc = sb.Popen(
		["samtools", "view", "-b", "-h", "-o", rescued_reads_path, "-"],
		stdin=sb.PIPE,
		universal_newlines=True,
	)
	proc.stdin.write(header)
	return (proc.stdin, proc)


def _close_rescued_reads_writer(writer, proc, rescued_reads_path):
	"""
	Close the rescued-read SAM/BAM writer and verify samtools completed.
	"""
	if writer is None:
		return

	writer.close()
	if proc is not None:
		return_code = proc.wait()
		if return_code != 0:
			raise Exception(
				"samtools failed while writing rescued read BAM "
				+ rescued_reads_path
			)


def _quality_is_below_threshold(mean_quality, min_mean_quality):
	if mean_quality is None:
		return False
	try:
		if np.isnan(mean_quality):
			return False
	except TypeError:
		pass
	return mean_quality < min_mean_quality


def _classify_amplicon_assignment(
	amp1,
	amp2,
	amp1_aln,
	amp2_aln,
	require_strict_amplicon_alignment=False,
	r1_mean_quality=None,
	r2_mean_quality=None,
	partial_rescue_min_mean_read_quality=PARTIAL_RESCUE_MIN_MEAN_READ_QUALITY_DEFAULT,
):
	"""
	Classify one read-pair amplicon assignment from primer and alignment calls.
	"""
	result = {
		"accepted_amplicon": None,
		"rescued_by_partial_alignment": False,
		"would_rescue_under_strict": False,
		"reject_reason": None,
	}

	if amp1 != amp2 or amp1 == "NA":
		result["reject_reason"] = "primer_disagreement_or_missing"
		return result

	alignment_calls = [x for x in (amp1_aln, amp2_aln) if x != "NA"]
	if any(x != amp1 for x in alignment_calls):
		result["reject_reason"] = "contradictory_alignment"
		return result

	if not alignment_calls:
		result["reject_reason"] = "no_inward_boundary_support"
		return result

	current_rescue_used = len(alignment_calls) < 2
	if current_rescue_used and partial_rescue_min_mean_read_quality is not None:
		if (
			_quality_is_below_threshold(r1_mean_quality, partial_rescue_min_mean_read_quality)
			or _quality_is_below_threshold(r2_mean_quality, partial_rescue_min_mean_read_quality)
		):
			result["reject_reason"] = "low_mean_quality"
			return result

	if require_strict_amplicon_alignment:
		strict_amplicon = None
		if amp1 == amp2 == amp1_aln == amp2_aln and amp1 != "NA":
			strict_amplicon = amp1
		result["accepted_amplicon"] = strict_amplicon
		result["would_rescue_under_strict"] = current_rescue_used and strict_amplicon is None
		if strict_amplicon is None:
			result["reject_reason"] = "strict_alignment_required"
		return result

	result["accepted_amplicon"] = amp1
	result["rescued_by_partial_alignment"] = current_rescue_used
	return result


def split_reads_by_amplicon(aligned_bam, output_root,amplicon_file,alt_alleles_file,primer_lookup_len,amp_file_dir,bowtie2_index,adapter_DNA,n_processes,keep_intermediate_files, reads_per_cell, min_total_reads_per_barcode, assign_reads_to_all_possible_amplicons=False, debug_rescued_reads_bam="", debug_require_strict_amplicon_alignment=False, debug_rejected_rescue_reads_bam="", partial_rescue_min_mean_read_quality=PARTIAL_RESCUE_MIN_MEAN_READ_QUALITY_DEFAULT, cache_manager=None):
	"""
	Split reads from a name-sorted aligned BAM into per-amplicon FASTQ files.

	Behavior summary
	----------------
	- Build primer lookup tables and optionally align amplicons to the genome.
	- Iterate through the name-sorted BAM and assign each read pair to one or
	  more amplicons based on primer matches and inward-facing amplicon-boundary
	  alignment support.
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
	debug_rescued_reads_bam : str
		Optional path to write read pairs accepted by partial alignment rescue.
		Use a .sam suffix for SAM output; any other suffix writes BAM.
	debug_require_strict_amplicon_alignment : bool
		If True, disable partial-alignment rescue and require both primer calls
		and both alignment-side calls to agree on the same amplicon.
	debug_rejected_rescue_reads_bam : str
		Optional path to write primer-agreeing rescue candidates rejected by the
		quality gate, contradictory alignment evidence, or missing inward
		boundary support.
	partial_rescue_min_mean_read_quality : float
		Minimum mean Phred quality required for each mate in partial rescues.
		Fully alignment-supported assignments are not gated by this setting.

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
	if assign_reads_to_all_possible_amplicons and debug_require_strict_amplicon_alignment:
		raise ValueError("debug_require_strict_amplicon_alignment cannot be used with assign_reads_to_all_possible_amplicons")

	info_file = output_root+".splitReads.ampliconInfo.txt"
	cache_record = None
	previous_split_record = None
	if cache_manager is not None:
		if cache_manager.config.mode.value != "disabled":
			previous_split_record = cache_manager.load("split_reads")
		cache_record = _build_split_cache_record(
			cache_manager, aligned_bam, amplicon_file, alt_alleles_file,
			bowtie2_index, primer_lookup_len, adapter_DNA,
			min_total_reads_per_barcode, assign_reads_to_all_possible_amplicons,
			debug_rescued_reads_bam, debug_require_strict_amplicon_alignment,
			debug_rejected_rescue_reads_bam, partial_rescue_min_mean_read_quality,
		)
		cached_information = None
		if os.path.isfile(info_file):
			try:
				_cache_valid, cached_names, cached_information = _load_split_read_cache(info_file, amp_file_dir)
			except (OSError, ValueError, IndexError):
				cached_information = None
		cache_requirements = _split_cache_requirements(
			output_root, amp_file_dir, info_file, cached_information,
			debug_rescued_reads_bam, debug_rejected_rescue_reads_bam,
		)
		if cache_manager.evaluate(cache_record, cache_requirements).is_hit:
			logging.info("Finished splitting reads from validated cache")
			return cached_names, cached_information, info_file
	elif os.path.isfile(info_file):
		cache_is_valid, amplicon_names, amplicon_information = _load_split_read_cache(info_file, amp_file_dir)
		if cache_is_valid:
			logging.info ("Finished splitting reads")
			return amplicon_names,amplicon_information,info_file
		os.remove(info_file)

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
	align_command = ["bowtie2", "-k", "2", "-x", bowtie2_index, "-p", str(n_processes), "-f", "-U", amplicon_fasta_file]
	with open(bowtie2_log,'w') as bt2log:
		bt2log.write(_command_to_string(align_command) + "\n")
		with open(aligned_amps_file, "w") as aligned_amps:
			completed = sb.run(align_command, stdout=aligned_amps, stderr=bt2log)
	if completed.returncode != 0:
		_raise_command_error(align_command, completed.returncode, context=f"see {bowtie2_log}")

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
	would_rescue_strict_alignment_count = 0 # strict debug mode only: reads rejected that current rescue would have accepted
	rejected_rescue_low_quality_count = 0
	rejected_rescue_contradictory_alignment_count = 0
	rejected_rescue_no_inward_boundary_support_count = 0
	unmapped_reads_count = 0 #unmapped by bowtie
	bam_iter = get_command_output(["samtools", "view", aligned_bam])#read in the aligned bam file
	rescued_reads_writer, rescued_reads_proc = _open_rescued_reads_writer(aligned_bam, debug_rescued_reads_bam)
	if debug_rescued_reads_bam:
		logging.info("Writing partial-alignment rescued reads to " + debug_rescued_reads_bam)
	rejected_rescue_reads_writer, rejected_rescue_reads_proc = _open_rescued_reads_writer(aligned_bam, debug_rejected_rescue_reads_bam)
	if debug_rejected_rescue_reads_bam:
		logging.info("Writing rejected rescue candidate reads to " + debug_rejected_rescue_reads_bam)
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
		if line1_rc:
			primer1 = seq1[(-1*primer_lookup_len):]
			seq1_to_write = reverse_complement(seq1)

		seq2_to_write = seq2
		primer2 = seq2[0:primer_lookup_len]
		if line2_rc:
			primer2 = seq2[(-1*primer_lookup_len):]#after alignment, all reads are put on the forward strand so we don't need to reverse complement them
			seq2_to_write = reverse_complement(seq2)


		amp1 = "NA"
		amp2 = "NA"
		if primer1 in primer_seqs:
			amp1 = primer_seqs[primer1]
		if primer2 in primer_seqs:
			amp2 = primer_seqs[primer2]
		amp1_aln = inward_alignment_amplicon(
			line1_chr,
			line1_start,
			line1_cigar,
			bool(line1_rc),
			start_amplicon_locs,
			end_amplicon_locs,
		)
		amp2_aln = inward_alignment_amplicon(
			line2_chr,
			line2_start,
			line2_cigar,
			bool(line2_rc),
			start_amplicon_locs,
			end_amplicon_locs,
		)

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
			assignment = _classify_amplicon_assignment(
				amp1,
				amp2,
				amp1_aln,
				amp2_aln,
				require_strict_amplicon_alignment = debug_require_strict_amplicon_alignment,
				r1_mean_quality = mean_phred_quality(qual1),
				r2_mean_quality = mean_phred_quality(qual2),
				partial_rescue_min_mean_read_quality = partial_rescue_min_mean_read_quality,
			)
			accepted_amplicon = assignment["accepted_amplicon"]
			reject_reason = assignment["reject_reason"]
			if assignment["rescued_by_partial_alignment"]:
				partial_alignment_rescue_count += 1
				if rescued_reads_writer is not None:
					rescued_reads_writer.write(line1)
					rescued_reads_writer.write(line2)
			if assignment["would_rescue_under_strict"]:
				would_rescue_strict_alignment_count += 1
			if reject_reason == "low_mean_quality":
				rejected_rescue_low_quality_count += 1
			elif reject_reason == "contradictory_alignment":
				rejected_rescue_contradictory_alignment_count += 1
			elif reject_reason == "no_inward_boundary_support":
				rejected_rescue_no_inward_boundary_support_count += 1
			if reject_reason in ("low_mean_quality", "contradictory_alignment", "no_inward_boundary_support"):
				if rejected_rescue_reads_writer is not None:
					rejected_rescue_reads_writer.write(line1)
					rejected_rescue_reads_writer.write(line2)

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
				aln_barcode_count[barcode] += 1
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
	_close_rescued_reads_writer(rescued_reads_writer, rescued_reads_proc, debug_rescued_reads_bam)
	_close_rescued_reads_writer(rejected_rescue_reads_writer, rejected_rescue_reads_proc, debug_rejected_rescue_reads_bam)
	for amp_name in amp_filehandles:
		amp_filehandles[amp_name][0].close()
		amp_filehandles[amp_name][1].close()

	outputs = OutputContext(output_root)
	identified_amplicon_file = outputs.path("valid_amplicons")
	with open(identified_amplicon_file,'w') as fout:
		for amp in amplicon_names:
			amp_aln_count = amplicon_count[amp]
			amplicon_information[amp]['aln_count'] = str(amp_aln_count)
			fout.write("%s\t%d\n"%(amp,amp_aln_count))


	logging.info("Cells not meeting the cell requirement: " + str(len(failing_barcode)))
	logging.info("Total Barcode Count: " + str(len(tot_barcode_count)))

	aligned_read_file = outputs.path("aligned_read_counts")
	logging.info("Writing aligned reads to " + aligned_read_file)
	with open(aligned_read_file,'w') as fout:
		fout.write("Barcode\tAligned Count\n")
		for barcode in aln_barcode_count:
			fout.write("%s\t%d\n"%(barcode,aln_barcode_count[barcode]))
	logging.info("Finished writing aligned read file")

	unaligned_read_file = outputs.path("unaligned_read_counts")
	logging.info('Writing unaligned reads to ' + unaligned_read_file)
	with open(unaligned_read_file,'w') as fout:
		fout.write("Barcode\tUnaligned Count\n")
		for barcode in unaln_barcode_count:
			fout.write("%s\t%d\n"%(barcode,unaln_barcode_count[barcode]))
	logging.info("Finished writing unaligned read file")

	identified_amplicon_aln_file = outputs.path("amplicon_classification")
	with open(identified_amplicon_aln_file,'w') as fout:
		fout.write("\t".join(['is_valid','amp1_from_seq','amp2_from_seq','amp1_from_align','amp2_from_align'])+"\n")
		for amp_key in sorted(aln_count.keys()):
			fout.write("%s\t%d\n"%("\t".join(amp_key),aln_count[amp_key]))

	amp_filenames.append((unidentified_out1_name,unidentified_out2_name))

	#zip output
	amp_commands = []
	for amp_filename_r1,amp_filename_r2 in amp_filenames:
		amp_commands.append(["gzip", "-f", amp_filename_r1])
		amp_commands.append(["gzip", "-f", amp_filename_r2])


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

		if debug_require_strict_amplicon_alignment:
			rescue_log_line = "      " + str(would_rescue_strict_alignment_count) + " read pairs would have been rescued by partial-alignment rescue but were excluded by strict alignment mode\n"
		else:
			rescue_log_line = "      " + str(partial_alignment_rescue_count) + " read pairs were rescued because both primer calls agreed and one inward-facing alignment boundary supported the amplicon without contradiction\n"
		rejected_rescue_log_line = \
				"      " + str(rejected_rescue_low_quality_count) + " rescue candidate read pairs were rejected due to low mean read quality\n" + \
				"      " + str(rejected_rescue_contradictory_alignment_count) + " rescue candidate read pairs were rejected due to contradictory alignment evidence\n" + \
				"      " + str(rejected_rescue_no_inward_boundary_support_count) + " rescue candidate read pairs had primer agreement but no inward-facing boundary support\n"
		log_str = str(tot_reads_count) + " alignment pairs were processed (including multi-mapped and unaligned reads)\n" + \
				"  Of these, " + str(aln_reads_count) + " read pairs were aligned to the genome\n" + \
				"    Of these, " + str(id_reads_count) + " read pairs were aligned correctly and identified\n" + \
				rescue_log_line + \
				rejected_rescue_log_line + \
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
		safe_remove(amplicon_fasta_file, silent=True)
		safe_remove(aligned_amps_file, silent=True)

	if cache_manager is not None:
		final_requirements = _split_cache_requirements(
				output_root, amp_file_dir, info_file, amplicon_information,
				debug_rescued_reads_bam, debug_rejected_rescue_reads_bam,
			)
		_prune_obsolete_split_outputs(
			cache_manager, previous_split_record, final_requirements, amp_file_dir,
		)
		cache_manager.commit(cache_record, final_requirements)

	return amplicon_names,amplicon_information,info_file
