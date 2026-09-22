from types import SimpleNamespace

import pytest

from CRISPRSCope import amplicon_assignment, fastq_processing
from CRISPRSCope.cache import (
    CacheManager,
    OutputRequirement,
    large_file_fingerprint,
    small_file_fingerprint,
)
from CRISPRSCope.paths import STAGE_ALIGN, build_stage_filename


def _tool(_command):
    return {"path": "/tools/fake", "version": "1.0"}


def _seed_parse_align_cache(tmp_path, monkeypatch):
    output_root = str(tmp_path / "run")
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    barcodes = tmp_path / "barcodes.txt"
    index = tmp_path / "genome"
    for path, content in ((r1, "r1"), (r2, "r2"), (barcodes, "AAAA\n")):
        path.write_text(content)
    (tmp_path / "genome.1.bt2").write_text("index")
    aligned_bam = build_stage_filename(
        STAGE_ALIGN, "align_merged", ext="bam", output_root=output_root
    )
    with open(aligned_bam, "wb") as handle:
        handle.write(b"bam")
    cell_file = output_root + ".parseReads.cellCount.txt"
    with open(cell_file, "w") as handle:
        handle.write("cellA\t7\n")

    monkeypatch.setattr(fastq_processing, "tool_identity", _tool)
    monkeypatch.setattr(
        "CRISPRSCope.cache.subprocess.run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=0),
    )
    manager = CacheManager(output_root)
    record = manager.new_record(
        "parse_align",
        algorithm_version=1,
        inputs={
            "r1": [large_file_fingerprint(r1)],
            "r2": [large_file_fingerprint(r2)],
            "barcodes": small_file_fingerprint(barcodes),
            "bowtie2_index": [large_file_fingerprint(tmp_path / "genome.1.bt2")],
        },
        parameters={
            "constant1": "AAAA",
            "constant2": "CCCC",
            "allow_barcode_mismatches": False,
            "adapter_DNA": "TTTT",
        },
        tools={"bowtie2": _tool(None), "samtools": _tool(None)},
    )
    requirements = (
        OutputRequirement("aligned_bam", aligned_bam, strategy="stat", validator="bam"),
        OutputRequirement(
            "cell_counts",
            cell_file,
            strategy="sha256",
            allow_empty=True,
            validator="cell_counts",
        ),
    )
    manager.commit(record, requirements)
    return manager, r1, r2, barcodes, index, output_root


def test_parse_align_cache_hit_returns_counts_without_processing(tmp_path, monkeypatch):
    manager, r1, r2, barcodes, index, output_root = _seed_parse_align_cache(
        tmp_path, monkeypatch
    )
    monkeypatch.setattr(
        fastq_processing,
        "get_valid_barcodes",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("cache miss")),
    )

    _bam, counts = fastq_processing.parse_and_align_reads(
        str(r1),
        str(r2),
        "AAAA",
        "CCCC",
        output_root,
        str(barcodes),
        False,
        "TTTT",
        str(index),
        1,
        cache_manager=manager,
    )

    assert counts == {"cellA": 7}
    assert manager.events[-1]["status"] == "hit"


def test_parse_align_changed_parameter_invalidates_cache(tmp_path, monkeypatch):
    manager, r1, r2, barcodes, index, output_root = _seed_parse_align_cache(
        tmp_path, monkeypatch
    )
    monkeypatch.setattr(
        fastq_processing,
        "get_valid_barcodes",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("recompute")),
    )

    with pytest.raises(RuntimeError, match="recompute"):
        fastq_processing.parse_and_align_reads(
            str(r1),
            str(r2),
            "CHANGED",
            "CCCC",
            output_root,
            str(barcodes),
            False,
            "TTTT",
            str(index),
            1,
            cache_manager=manager,
        )
    assert manager.events[-1]["status"] == "invalid"


def _seed_split_cache(tmp_path, monkeypatch):
    output_root = str(tmp_path / "run")
    amp_dir = tmp_path / "run.seq_by_amplicon"
    amp_dir.mkdir()
    aligned_bam = tmp_path / "aligned.bam"
    aligned_bam.write_bytes(b"bam")
    amplicons = tmp_path / "amplicons.tsv"
    amplicons.write_text("ampA\tACGTACGT\tNA\n")
    index = tmp_path / "genome"
    (tmp_path / "genome.1.bt2").write_text("index")
    info_file = tmp_path / "run.splitReads.ampliconInfo.txt"
    r1 = amp_dir / "03_reads_all_cells.ampA.r1.fq.gz"
    r2 = amp_dir / "03_reads_all_cells.ampA.r2.fq.gz"
    import gzip

    for path in (r1, r2):
        with gzip.open(path, "wt") as handle:
            handle.write("")
    header = ["name", "aln_count", "reads_r1_file", "reads_r2_file"]
    info_file.write_text(
        "\t".join(header) + "\n" + f"ampA\t1\t{r1}\t{r2}\n"
    )
    (tmp_path / "run.splitReads.valid_amps.txt").write_text("ampA\t1\n")
    (tmp_path / "run.splitReads.aligned.txt").write_text("Barcode\tAligned Count\ncellA\t1\n")
    (tmp_path / "run.splitReads.unaligned.txt").write_text("Barcode\tUnaligned Count\n")
    (tmp_path / "run.splitReads.amp_classification.txt").write_text(
        "is_valid\tamp1_from_seq\tamp2_from_seq\tamp1_from_align\tamp2_from_align\n"
    )
    monkeypatch.setattr(amplicon_assignment, "tool_identity", _tool)
    manager = CacheManager(output_root)
    record = amplicon_assignment._build_split_cache_record(
        manager,
        str(aligned_bam),
        str(amplicons),
        "",
        str(index),
        18,
        "ADAPTER",
        10,
        False,
        "",
        False,
        "",
        30.0,
    )
    information = {
        "ampA": {
            "name": "ampA",
            "aln_count": "1",
            "reads_r1_file": str(r1),
            "reads_r2_file": str(r2),
        }
    }
    requirements = amplicon_assignment._split_cache_requirements(
        output_root, str(amp_dir), str(info_file), information
    )
    manager.commit(record, requirements)
    return manager, aligned_bam, amplicons, index, amp_dir, output_root


def test_split_cache_hit_uses_validated_amplicon_fastqs(tmp_path, monkeypatch):
    manager, aligned_bam, amplicons, index, amp_dir, output_root = _seed_split_cache(
        tmp_path, monkeypatch
    )
    names, information, _info = amplicon_assignment.split_reads_by_amplicon(
        str(aligned_bam),
        output_root,
        str(amplicons),
        "",
        18,
        str(amp_dir),
        str(index),
        "ADAPTER",
        1,
        False,
        {"cellA": 1},
        10,
        cache_manager=manager,
    )
    assert names == ["ampA"]
    assert information["ampA"]["reads_r1_file"].endswith(".fq.gz")
    assert manager.events[-1]["status"] == "hit"


def test_split_cache_changed_setting_recomputes(tmp_path, monkeypatch):
    manager, aligned_bam, amplicons, index, amp_dir, output_root = _seed_split_cache(
        tmp_path, monkeypatch
    )
    monkeypatch.setattr(
        amplicon_assignment.sb,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=1),
    )
    with pytest.raises(Exception, match="External command failed"):
        amplicon_assignment.split_reads_by_amplicon(
            str(aligned_bam),
            output_root,
            str(amplicons),
            "",
            19,
            str(amp_dir),
            str(index),
            "ADAPTER",
            1,
            False,
            {"cellA": 1},
            10,
            cache_manager=manager,
        )
    assert manager.events[-1]["status"] == "invalid"
