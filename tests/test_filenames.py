import os
from CRISPRSCope.cli import build_stage_filename, validate_output_root
import pytest


def test_build_stage_filename_basic_prefix_mode(tmp_path):
    """
    Prefix mode: output_root does NOT end with '/' and the directory does not exist as a directory.
    Expect files to be written into the parent (tmp_path) with basename(output_root) as prefix.
    """
    root_prefix = str(tmp_path / "demo_prefix")
    validated = validate_output_root(root_prefix)  # returns abs path (parent is tmp_path)
    fname = build_stage_filename(3, "amp", amplicon="AMP:1", read="r1", ext="fastq.gz", output_root=validated)
    assert os.path.dirname(fname) == str(tmp_path)
    bn = os.path.basename(fname)
    assert bn.startswith("demo_prefix_03_amp")
    assert ".AMP_1." in bn
    assert bn.endswith(".fastq.gz")


def test_build_stage_filename_dir_mode(tmp_path):
    """
    Directory mode: output_root ends with '/' or is an existing directory.
    Expect files to be written *inside* that directory and file base to start with stage tag.
    """
    d = str(tmp_path / "demo_dir")
    os.makedirs(d, exist_ok=True)
    # Pass a trailing slash to explicitly request directory mode
    validated = validate_output_root(d + os.sep)
    fname = build_stage_filename(3, "amp", amplicon="AMP:1", read="r1", ext="fastq.gz", output_root=validated)
    assert os.path.dirname(fname) == d
    bn = os.path.basename(fname)
    # in dir mode the filename starts with the stage tag ("03_amp...")
    assert bn.startswith("03_amp")
    assert ".AMP_1." in bn
    assert bn.endswith(".fastq.gz")


def test_build_stage_filename_invalid_stage(tmp_path):
    root = str(tmp_path / "demo_invalid")
    validated = validate_output_root(root)
    with pytest.raises(ValueError):
        build_stage_filename(-1, "raw", output_root=validated)


def test_build_stage_filename_requires_output_root():
    with pytest.raises(ValueError):
        build_stage_filename(1, "raw")  # missing required output_root
