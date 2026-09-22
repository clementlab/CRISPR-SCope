from types import SimpleNamespace

import pytest

from CRISPRSCope import fastq_processing
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
