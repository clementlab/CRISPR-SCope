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


class ExternalCommandError(RuntimeError):
	"""Raised when an external command exits unsuccessfully."""


def _command_to_string(command):
	if isinstance(command, (list, tuple)):
		return shlex.join([str(x) for x in command])
	return str(command)


def _raise_command_error(command, returncode=None, context=None, stderr=None):
	message = "External command failed"
	if returncode is not None:
		message += f" with return code {returncode}"
	message += f": {_command_to_string(command)}"
	if context:
		message += f" ({context})"
	if stderr:
		message += f"\n{stderr.strip()}"
	raise ExternalCommandError(message)


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
