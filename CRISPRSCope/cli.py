"""Stable command-line entry point and compatibility facade.

Implementation lives in focused pipeline-stage modules. Existing imports from
``CRISPRSCope.cli`` remain available during this compatibility-preserving refactor.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import subprocess as sb
import sys
from pathlib import Path

import matplotlib.pyplot as plt

# ``python path/to/CRISPRSCope/cli.py settings.txt`` is a supported legacy
# invocation.  Direct-file execution has no package context, so establish one
# before importing the façade's relative modules.  Module execution
# (``python -m CRISPRSCope.cli``) already provides this context.
if __package__ in (None, ""):
    package_parent = str(Path(__file__).resolve().parent.parent)
    if package_parent not in sys.path:
        sys.path.insert(0, package_parent)
    __package__ = "CRISPRSCope"

from CRISPRSCope import __version__
from . import amplicon_assignment as _amplicon_assignment
from . import crispresso as _crispresso
from . import pipeline as _pipeline
from .amplicon_assignment import *
from .crispresso import *
from .fastq_processing import *
from .paths import *
from .plots_and_report import *
from .settings import *
from .summaries import *
from .pipeline import (
    _require_selected_barcodes,
    write_editing_rate_ci_output,
    write_editing_rate_depth_stability_output,
    write_h5ad_output,
)
from .amplicon_assignment import (
    _classify_amplicon_assignment,
    _close_rescued_reads_writer,
    _load_split_read_cache,
    _open_rescued_reads_writer,
    _quality_is_below_threshold,
)
from .crispresso import (
    _barcode_set_sha256,
    _crispresso_annotation_is_modified,
    _decompressed_fastq_sha256,
    _filtered_allele_fastq_path,
    _load_final_allele_read_support,
    _parse_cache_matches_ignore_substitutions,
    _parse_filtered_crispresso_allele_output,
    _write_filtered_summary_table,
)
from .plots_and_report import _numeric_mod_pct_columns, _numeric_tot_count_columns
from .settings import (
    _build_input_ref_names,
    _normalize_optional_guide,
    _parse_amplicon_score_config,
    _parse_bool_setting,
    _parse_editing_rate_ci_config,
    _parse_editing_rate_depth_stability_config,
    _parse_float_setting,
    _parse_int_setting,
    _parse_settings_file,
    _resolve_existing_fastq_path,
    _resolve_settings_path,
    _resolve_settings_path_list,
    _settings_value_is_path,
)

_ACTIVE_OUTPUT_MANIFEST = None

def _set_active_manifest(manifest):
    global _ACTIVE_OUTPUT_MANIFEST
    _ACTIVE_OUTPUT_MANIFEST = manifest

def _main_impl():
    return _pipeline.run_pipeline(manifest_observer=_set_active_manifest)


def parse_crispresso_outputs(*args, **kwargs):
    """Compatibility wrapper honoring legacy ``cli`` monkeypatches."""
    _crispresso.parse_one_crispresso_output = parse_one_crispresso_output
    return _crispresso.parse_crispresso_outputs(*args, **kwargs)


def run_crispresso_commands(*args, **kwargs):
    """Compatibility wrapper honoring legacy ``cli`` monkeypatches."""
    _crispresso.run_crispresso_command = run_crispresso_command
    return _crispresso.run_crispresso_commands(*args, **kwargs)


def split_reads_by_amplicon(*args, **kwargs):
    """Compatibility wrapper honoring legacy ``cli`` monkeypatches."""
    _amplicon_assignment.get_command_output = get_command_output
    return _amplicon_assignment.split_reads_by_amplicon(*args, **kwargs)

def main():
    """Run the pipeline and, when requested, write its diagnostic manifest."""
    global _ACTIVE_OUTPUT_MANIFEST
    _ACTIVE_OUTPUT_MANIFEST = None
    try:
        result = _main_impl()
    except BaseException as error:
        manifest = _ACTIVE_OUTPUT_MANIFEST
        if manifest is not None:
            manifest.fail(manifest.active_stage or "initialization", error)
            try:
                manifest.write()
            except BaseException:
                import logging
                logging.exception("Failed to write output manifest after pipeline failure")
        raise
    else:
        manifest = _ACTIVE_OUTPUT_MANIFEST
        if manifest is not None:
            manifest.complete()
            manifest.write()
        return result
    finally:
        _ACTIVE_OUTPUT_MANIFEST = None

if __name__ == "__main__":
    main()
