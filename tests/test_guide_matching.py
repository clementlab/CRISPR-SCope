"""Guide preflight tests use the CRISPResso2 installed with this project."""
import logging
import gzip
from types import SimpleNamespace

import pytest

from CRISPRSCope import crispresso
from CRISPRSCope import pipeline, settings as settings_module
from CRISPRSCope.guide_matching import GuideResolution, validate_guides
from CRISPRSCope.settings import _parse_flexiguide_homology


GUIDE = "ACGTACGTACGTACGTACGT"
MISMATCHED_LOCUS = "ACGTACGTACGTACGTACGA"


def _amplicons(tmp_path, *rows):
    path = tmp_path / "amplicons.txt"
    path.write_text("\n".join(rows) + "\n")
    return str(path)


def test_exact_guides_match_forward_and_reverse_strands(tmp_path):
    from CRISPResso2 import CRISPRessoShared

    path = _amplicons(tmp_path,
        f"forward\tTTTT{GUIDE}TTTT\t{GUIDE}\t2",
        f"reverse\tTTTT{CRISPRessoShared.reverse_complement(GUIDE)}TTTT\t{GUIDE}\t2\texact",
    )
    resolved = validate_guides(path)
    assert resolved["forward"].exact == (GUIDE,)
    assert resolved["reverse"].exact == (GUIDE,)


def test_exact_mismatches_are_collected_with_actionable_locations(tmp_path):
    path = _amplicons(tmp_path,
        f"targetA\tTTTT{MISMATCHED_LOCUS}TTTT\t{GUIDE}\t2",
        f"targetB\tTTTT{MISMATCHED_LOCUS}TTTT\t{GUIDE}\t2\texact",
    )
    with pytest.raises(ValueError, match="Guide validation failed before read processing") as error:
        validate_guides(path)
    assert "targetA (line 1)" in str(error.value)
    assert "targetB (line 2)" in str(error.value)
    assert GUIDE in str(error.value)
    assert "guide_match_mode to flexible" in str(error.value)


def test_flexible_mismatch_threshold_and_omission(tmp_path, caplog):
    path = _amplicons(tmp_path,
        f"offtarget\tTTTT{MISMATCHED_LOCUS}TTTT\t{GUIDE}\t2\tflexible",
    )
    assert validate_guides(path, homology=94)["offtarget"].flexible == (GUIDE,)
    with caplog.at_level(logging.WARNING):
        result = validate_guides(path, homology=95)["offtarget"]
    assert result.omitted == (GUIDE,)
    assert not result.exact and not result.flexible
    assert "quantify the full amplicon" in caplog.text


def test_multiple_guides_and_reference_sequences_share_row_mode(tmp_path, caplog):
    path = _amplicons(tmp_path,
        f"mixed\tTTTT{GUIDE}TTTT,TTTT{MISMATCHED_LOCUS}TTTT\t{GUIDE},TTTTTTTTTTTTTTTTTTTT\t2\tflexible",
    )
    with caplog.at_level(logging.WARNING):
        result = validate_guides(path)["mixed"]
    assert result.exact == (GUIDE,)
    assert result.omitted == ("TTTTTTTTTTTTTTTTTTTT",)
    assert "omitting this guide" in caplog.text
    assert "quantify the full amplicon" not in caplog.text


def test_alternate_allele_references_and_blank_guides(tmp_path):
    path = _amplicons(tmp_path,
        f"alt\tTTTT{MISMATCHED_LOCUS}TTTT\t{GUIDE}\t2",
        "blank\tACGTACGT\tNA\t2",
    )
    alt = tmp_path / "alternate.txt"
    alt.write_text(f"name\talleles\treferences\nalt\t0\tTTTT{GUIDE}TTTT\n")
    result = validate_guides(path, str(alt))
    assert result["alt"].exact == (GUIDE,)
    assert not result["blank"].exact


@pytest.mark.parametrize("value", ["-1", "101", "abc", "80.5"])
def test_invalid_homology_setting(tmp_path, value):
    settings = tmp_path / "settings.txt"
    settings.write_text(f"flexiguide_homology\t{value}\n")
    with pytest.raises(ValueError, match="flexiguide_homology"):
        _parse_flexiguide_homology(str(settings))


def test_homology_default_and_override(tmp_path):
    settings = tmp_path / "settings.txt"
    settings.write_text("# empty settings\n")
    assert _parse_flexiguide_homology(str(settings)) == 80
    settings.write_text("flexiguide_homology\t90\n")
    assert _parse_flexiguide_homology(str(settings)) == 90


@pytest.mark.parametrize("alleles", [False, True])
def test_resolved_guide_arguments_are_used_in_both_passes(tmp_path, alleles, monkeypatch):
    output_root = str(tmp_path / "run")
    base_dir = tmp_path / "run.crispresso"
    actual_dir = tmp_path / ("run.crispresso.filtered" if alleles else "run.crispresso")
    actual_dir.mkdir()
    seq_dir = tmp_path / "run.seq_by_amplicon"
    seq_dir.mkdir()
    resolution = GuideResolution("flexible", (GUIDE,), (MISMATCHED_LOCUS,), (), 91)
    info = {
        "amp_seqs": f"TTTT{GUIDE}TTTT,TTTT{MISMATCHED_LOCUS}TTTT",
        "guide_seq": f"{GUIDE},{MISMATCHED_LOCUS}",
        "guide_resolution": resolution,
        "aln_count": "1",
    }
    if alleles:
        allele_fastq = seq_dir / "04_alleles_qc_cells.ampA.fq.gz"
        with gzip.open(allele_fastq, "wt") as handle:
            handle.write("@r\nACGT\n+\nIIII\n")
    else:
        for read in ("r1", "r2"):
            path = seq_dir / f"ampA.{read}.fq.gz"
            path.write_text("input\n")
            info[f"reads_{read}_file"] = str(path)
    class Pool:
        def __init__(self, *_args):
            pass

        def map_async(self, _function, jobs):
            for job in jobs:
                folder = actual_dir / "CRISPResso_on_ampA"
                folder.mkdir(exist_ok=True)
                (folder / "CRISPResso2_info.json").write_text("{}\n")
                (actual_dir / "ampA.finished").write_text("")
            return SimpleNamespace(get=lambda *_args: [
                {"returncode": 0, "error": None, "command": job["command"]}
                for job in jobs
            ])

        def close(self):
            pass

        def join(self):
            pass

    monkeypatch.setattr(crispresso.mp, "Pool", Pool)
    result = crispresso.run_crispresso_commands(
        ["ampA"], {"ampA": info}, output_root, str(base_dir), False, 1,
        alleles=alleles,
    )
    command = result["ampA"]["crispresso_command"]
    assert f"-g {GUIDE}" in command
    assert f"--flexiguide_seq {MISMATCHED_LOCUS}" in command
    assert "--flexiguide_homology 91" in command


def test_crispresso_cache_key_changes_with_mode_and_threshold(tmp_path, monkeypatch):
    from CRISPRSCope.cache import CacheManager

    fastq = tmp_path / "reads.fastq.gz"
    with gzip.open(fastq, "wt") as handle:
        handle.write("@r\nACGT\n+\nIIII\n")
    monkeypatch.setattr(crispresso, "tool_identity", lambda *_args: {"path": "/tools/CRISPResso", "version": "test"})
    manager = CacheManager(str(tmp_path / "run"))
    info = {"amp_seqs": f"TTTT{GUIDE}TTTT", "guide_seq": GUIDE}

    def key():
        return crispresso._build_crispresso_cache_record(
            manager, "ampA", info, False, False, [str(fastq)]
        ).cache_key

    original = key()
    info["guide_resolution"] = GuideResolution("exact", (GUIDE,), (), (), 80)
    exact = key()
    info["guide_resolution"] = GuideResolution("flexible", (GUIDE,), (), (), 80)
    flexible = key()
    info["guide_resolution"] = GuideResolution("flexible", (GUIDE,), (), (), 90)
    changed_threshold = key()
    assert len({original, exact, flexible, changed_threshold}) == 4


def test_exact_mismatch_stops_pipeline_before_read_parsing(tmp_path, monkeypatch):
    import sys

    for filename in ("r1.fastq", "r2.fastq", "barcodes.txt", "genome.1.bt2"):
        (tmp_path / filename).write_text("input\n")
    amplicons = _amplicons(tmp_path,
        f"target\tTTTT{MISMATCHED_LOCUS}TTTT\t{GUIDE}\t2\texact"
    )
    settings_file = tmp_path / "settings.txt"
    settings_file.write_text("\n".join([
        f"r1\t{tmp_path / 'r1.fastq'}",
        f"r2\t{tmp_path / 'r2.fastq'}",
        f"barcodes\t{tmp_path / 'barcodes.txt'}",
        f"amplicons\t{amplicons}",
        f"bowtie2_index\t{tmp_path / 'genome'}",
        "constant1\tACGT",
        "constant2\tTGCA",
        "processes\t1",
    ]) + "\n")
    monkeypatch.setattr(sys, "argv", ["CRISPRSCope", str(settings_file)])
    monkeypatch.setattr(settings_module, "tool_identity", lambda *_args: {"path": "/tools/fake"})
    monkeypatch.setattr(
        pipeline, "parse_and_align_reads",
        lambda *_args, **_kwargs: pytest.fail("Read parsing began before guide validation"),
    )
    with pytest.raises(ValueError, match="target.*guide RNA"):
        pipeline.run_pipeline()


def test_omitted_flexible_guide_is_not_passed_to_crispresso(tmp_path, monkeypatch):
    output_root = str(tmp_path / "run")
    base_dir = tmp_path / "run.crispresso"
    base_dir.mkdir()
    reads = tmp_path / "reads.fastq"
    reads.write_text("input\n")

    class Pool:
        def __init__(self, *_args):
            pass

        def map_async(self, _function, jobs):
            for job in jobs:
                folder = base_dir / "CRISPResso_on_ampA"
                folder.mkdir()
                (folder / "CRISPResso2_info.json").write_text("{}\n")
                (base_dir / "ampA.finished").write_text("")
            return SimpleNamespace(get=lambda *_args: [
                {"returncode": 0, "error": None, "command": job["command"]}
                for job in jobs
            ])

        def close(self):
            pass

        def join(self):
            pass

    monkeypatch.setattr(crispresso.mp, "Pool", Pool)
    info = {
        "amp_seqs": f"TTTT{MISMATCHED_LOCUS}TTTT",
        "guide_seq": GUIDE,
        "guide_resolution": GuideResolution("flexible", (), (), (GUIDE,), 80),
        "aln_count": "1",
        "reads_r1_file": str(reads),
        "reads_r2_file": str(reads),
    }
    result = crispresso.run_crispresso_commands(
        ["ampA"], {"ampA": info}, output_root, str(base_dir), False, 1,
        alleles=False,
    )
    command = result["ampA"]["crispresso_command"]
    assert " -g " not in command
    assert "--flexiguide_seq" not in command
