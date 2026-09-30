"""Shared imports for CRISPRSCope pipeline-stage modules."""
from __future__ import annotations

import errno
import glob
import gzip
import hashlib
from itertools import zip_longest
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
import tempfile
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
	small_file_fingerprint,
	tool_identity,
)
from CRISPRSCope.io_utils import open_text_maybe_gzip
from CRISPRSCope.output_artifacts import OutputContext, OutputManifest

from .crispresso import run_command
from .paths import (
	STAGE_ALIGN,
	STAGE_PARSE,
	_command_to_string,
	_raise_command_error,
	build_stage_filename,
	safe_remove,
)

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
		for rec1, rec2 in zip_longest(f1, f2):
			if rec1 is None or rec2 is None:
				raise ValueError(
					"Paired FASTQ files contain different numbers of records: "
					f"{r1_path} and {r2_path}"
				)
			metrics.tot_reads += 1

			info1 = rec1.name.strip()
			info1_first_bit = info1.split(' ')[0]
			seq1 = str(rec1.sequence)
			plus1 = '+'
			qual1 = rec1.qualities

			info2 = rec2.name.strip()
			info2_first_bit = info2.split(' ')[0]
			if info1_first_bit != info2_first_bit:
				raise ValueError(
					"Paired FASTQ records have different read identifiers: "
					f"{info1_first_bit!r} != {info2_first_bit!r}"
				)
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
	log_file = out_name + ".bowtie2.log"
	bowtie_cmd = ["bowtie2", "-x", bowtie2_index, "-p", str(threads), "-1", r1_path, "-2", r2_path]
	samtools_cmd = ["samtools", "view", "-bS", "-"]
	logging.debug("running: %s | %s", _command_to_string(bowtie_cmd), _command_to_string(samtools_cmd))
	with open(log_file, "a") as bowtie_log, open(out_name, "wb") as bam_out:
		bowtie_proc = sb.Popen(bowtie_cmd, stdout=sb.PIPE, stderr=bowtie_log)
		samtools_proc = sb.Popen(samtools_cmd, stdin=bowtie_proc.stdout, stdout=bam_out, stderr=sb.PIPE)
		bowtie_proc.stdout.close()
		_, samtools_stderr = samtools_proc.communicate()
		bowtie_return = bowtie_proc.wait()
		if bowtie_return != 0:
			_raise_command_error(bowtie_cmd, bowtie_return, context=f"see {log_file}")
		if samtools_proc.returncode != 0:
			_raise_command_error(
				samtools_cmd,
				samtools_proc.returncode,
				context=f"creating {out_name}",
				stderr=samtools_stderr.decode(errors="replace") if samtools_stderr else None,
			)


def parse_and_align_reads(r1_fastqs,r2_fastqs,constant1,constant2,
						  output_root,barcode_file,allow_barcode_mismatches,
						  adapter_DNA,bowtie2_index,n_processes,keep_intermediate_files=False,
						  cache_manager=None):
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

	cache_record = None
	cache_requirements = (
		OutputRequirement("aligned_bam", aligned_bam, strategy="stat", validator="bam"),
		OutputRequirement(
			"cell_counts", cell_file, strategy="sha256", allow_empty=True,
			validator="cell_counts",
		),
	)
	if cache_manager is not None:
		index_files = sorted(
			path for path in glob.glob(bowtie2_index + ".*")
			if path.endswith((".bt2", ".bt2l"))
		)
		cache_record = cache_manager.new_record(
			"parse_align",
			algorithm_version=1,
			inputs={
				"r1": [large_file_fingerprint(path) for path in r1_fastqs.split(",")],
				"r2": [large_file_fingerprint(path) for path in r2_fastqs.split(",")],
				"barcodes": small_file_fingerprint(barcode_file),
				"bowtie2_index": [large_file_fingerprint(path) for path in index_files],
			},
			parameters={
				"constant1": constant1,
				"constant2": constant2,
				"allow_barcode_mismatches": bool(allow_barcode_mismatches),
				"adapter_DNA": adapter_DNA,
			},
			tools={
				"bowtie2": tool_identity(("bowtie2", "--version")),
				"samtools": tool_identity(("samtools", "--version")),
			},
		)
		if cache_manager.evaluate(cache_record, cache_requirements).is_hit:
			reads_per_cell = {}
			with open(cell_file, 'r') as fin:
				for line in fin:
					reads_per_cell[line.split('\t')[0]] = int(line.split('\t')[1].strip())
			logging.info("Finished parsing reads from validated cache")
			return (aligned_bam, reads_per_cell)
	elif os.path.isfile(info_file) and os.path.isfile(aligned_bam) and os.path.isfile(cell_file):
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
		input_bams = [result[2] for result in parsed_results]
		#inter_bam = output_root+".merged.bam"
		inter_bam = build_stage_filename(STAGE_ALIGN, "intermediate_merged", ext="bam", output_root=output_root)
		start_bam_cat = time.time()
		run_command(["samtools", "cat", "-o", inter_bam] + input_bams)
		end_bam_cat = time.time() - start_bam_cat
		logging.info("BAM cat ended in %.3f seconds", end_bam_cat)
		logging.info("Inside parsed_results > 1")
	if len(parsed_results) == 1:
		inter_bam = parsed_results[0][2]
		logging.info("Inside parsed_results == 1")

	bam_threads = min(12, n_processes)
	temporary_aligned_bam = aligned_bam + ".tmp.%d" % os.getpid()
	sort_bam_cmd = ["samtools", "sort", "-n", "-@", str(bam_threads), "-o", temporary_aligned_bam, inter_bam]

	start_bam_sort = time.time()
	try:
		run_command(sort_bam_cmd)
		os.replace(temporary_aligned_bam, aligned_bam)
	except BaseException:
		safe_remove(temporary_aligned_bam, silent=True)
		raise
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
	cell_parent = os.path.dirname(cell_file) or os.getcwd()
	fd, temporary_cell_file = tempfile.mkstemp(
		prefix="." + os.path.basename(cell_file) + ".", suffix=".tmp", dir=cell_parent
	)
	try:
		with os.fdopen(fd, 'w') as fout:
			for cell in merged_metrics.reads_per_cell:
				fout.write(f"{cell}\t{merged_metrics.reads_per_cell[cell]}\n")
		os.replace(temporary_cell_file, cell_file)
	except BaseException:
		safe_remove(temporary_cell_file, silent=True)
		raise

	if cache_manager is not None:
		cache_manager.commit(cache_record, cache_requirements)

	return (aligned_bam, merged_metrics.reads_per_cell)
