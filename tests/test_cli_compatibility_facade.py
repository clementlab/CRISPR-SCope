"""Compatibility guarantees for legacy imports from ``CRISPRSCope.cli``."""

import multiprocessing
import os
import subprocess
import sys

import matplotlib.pyplot as plt

from CRISPRSCope import __version__
from CRISPRSCope import cli
from CRISPRSCope import amplicon_assignment, crispresso, fastq_processing, paths
from CRISPRSCope import pipeline, plots_and_report, settings, summaries


def test_cli_reexports_stage_helpers_and_runtime_module_aliases():
    assert cli.AmpliconScoreConfig is settings.AmpliconScoreConfig
    assert cli.build_stage_filename is paths.build_stage_filename
    assert cli.parse_fq_file_pair is fastq_processing.parse_fq_file_pair
    assert callable(cli.split_reads_by_amplicon)
    assert cli.parse_one_crispresso_output is crispresso.parse_one_crispresso_output
    assert cli.generate_amplicon_score is summaries.generate_amplicon_score
    assert cli.make_report is plots_and_report.make_report
    assert cli.write_editing_rate_ci_output is pipeline.write_editing_rate_ci_output

    # Stage-to-stage helpers retain the exact owning implementation after extraction.
    assert fastq_processing.run_command is crispresso.run_command
    assert crispresso._build_input_ref_names is settings._build_input_ref_names
    assert crispresso.STAGE_SPLIT == paths.STAGE_SPLIT
    assert amplicon_assignment._raise_command_error is paths._raise_command_error
    assert amplicon_assignment.safe_write_path is paths.safe_write_path

    assert cli.os is os
    assert cli.sb is subprocess
    assert cli.mp is multiprocessing
    assert cli.plt is plt


def test_cli_keeps_main_proxy_for_monkeypatch_based_callers(monkeypatch):
    called = []
    monkeypatch.setattr(cli, "_main_impl", lambda: called.append("ran"))

    assert cli.main() is None
    assert called == ["ran"]


def test_cli_can_be_invoked_by_legacy_direct_file_path(tmp_path):
    result = subprocess.run(
        [sys.executable, cli.__file__, "--version"],
        cwd=tmp_path,
        check=True,
        text=True,
        capture_output=True,
    )

    assert __version__ in result.stdout
