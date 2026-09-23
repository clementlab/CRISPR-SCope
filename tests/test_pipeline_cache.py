import os
import gzip
from types import SimpleNamespace

import pandas as pd
import pytest

from CRISPRSCope import amplicon_assignment, crispresso, fastq_processing
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


def _seed_crispresso_cache(tmp_path, monkeypatch, *, alleles=False):
    output_root = str(tmp_path / "run")
    base_dir = tmp_path / "run.crispresso"
    run_dir = tmp_path / ("run.crispresso.filtered" if alleles else "run.crispresso")
    run_dir.mkdir(parents=True)
    amp_dir = tmp_path / "run.seq_by_amplicon"
    amp_dir.mkdir()
    if alleles:
        inputs = [amp_dir / "04_alleles_qc_cells.ampA.fq.gz"]
    else:
        inputs = [
            amp_dir / "03_reads_all_cells.ampA.r1.fq.gz",
            amp_dir / "03_reads_all_cells.ampA.r2.fq.gz",
        ]
    for path in inputs:
        path.write_text("fastq\n")
    information = {
        "ampA": {
            "name": "ampA",
            "aln_count": "1",
            "reads_r1_file": str(inputs[0]),
            "reads_r2_file": str(inputs[-1]),
            "amp_seqs": "ACGT",
            "guide_seq": "AC",
        }
    }
    folder = run_dir / "CRISPResso_on_ampA"
    folder.mkdir()
    (folder / "CRISPResso2_info.json").write_text("{}\n")
    (folder / "CRISPResso_output.fastq.gz").write_text(
        "@read\nACGT\n+\nIIII\n"
    )
    (run_dir / "CRISPResso_on_ampA.html").write_text("<html></html>\n")
    finished = run_dir / "ampA.finished"
    finished.write_text("")
    monkeypatch.setattr(crispresso, "tool_identity", _tool)
    manager = CacheManager(output_root)
    record = crispresso._build_crispresso_cache_record(
        manager, "ampA", information["ampA"], False, alleles, [str(path) for path in inputs]
    )
    manager.commit(
        record,
        crispresso._crispresso_cache_requirements(
            str(finished), str(folder), require_report=True
        ),
    )
    return manager, information, output_root, base_dir, folder


@pytest.mark.parametrize("alleles", [False, True])
def test_crispresso_per_amplicon_cache_hits_without_running(tmp_path, monkeypatch, alleles):
    manager, information, output_root, base_dir, folder = _seed_crispresso_cache(
        tmp_path, monkeypatch, alleles=alleles
    )
    monkeypatch.setattr(
        crispresso,
        "run_crispresso_command",
        lambda _job: (_ for _ in ()).throw(AssertionError("cache miss")),
    )
    if alleles:
        monkeypatch.setattr(
            crispresso,
            "_decompressed_fastq_sha256",
            lambda _path: (_ for _ in ()).throw(AssertionError("large file was hashed")),
        )

    result = crispresso.run_crispresso_commands(
        ["ampA"], information, output_root, str(base_dir), False, 1,
        alleles=alleles, cache_manager=manager,
    )
    assert result["ampA"]["status"] == "Completed"
    assert folder.exists()
    assert manager.events[-1]["status"] == "hit"


def test_crispresso_changed_guide_clears_and_reruns_exact_amplicon(tmp_path, monkeypatch):
    manager, information, output_root, base_dir, folder = _seed_crispresso_cache(
        tmp_path, monkeypatch, alleles=False
    )
    stale = folder / "stale.txt"
    stale.write_text("stale\n")
    information["ampA"]["guide_seq"] = "CHANGED"

    class Pool:
        def __init__(self, *_args):
            pass

        def map_async(self, function, jobs):
            values = []
            for job in jobs:
                assert "--no_rerun" not in job["args"]
                assert job["args"][job["args"].index("-o") + 1] == str(base_dir)
                assert job["args"][job["args"].index("-n") + 1] == "ampA"
                assert not stale.exists()
                os.makedirs(job["crispresso_run_folder"], exist_ok=True)
                with open(os.path.join(job["crispresso_run_folder"], "CRISPResso2_info.json"), "w") as handle:
                    handle.write("{}\n")
                with gzip.open(os.path.join(job["crispresso_run_folder"], "CRISPResso_output.fastq.gz"), "wt") as handle:
                    handle.write("fastq\n")
                with open(job["crispresso_run_folder"] + ".html", "w") as handle:
                    handle.write("<html></html>\n")
                with open(job["finished_file"], "w"):
                    pass
                values.append({"returncode": 0, "error": None, "command": job["command"]})
            return SimpleNamespace(get=lambda *_args: values)

        def close(self):
            pass

        def join(self):
            pass

    monkeypatch.setattr(crispresso.mp, "Pool", Pool)
    result = crispresso.run_crispresso_commands(
        ["ampA"], information, output_root, str(base_dir), False, 1,
        alleles=False, cache_manager=manager,
    )
    assert result["ampA"]["status"] == "Completed"
    assert manager.events[-1]["status"] == "invalid"


def test_crispresso_commits_later_success_after_incomplete_amplicon(tmp_path, monkeypatch):
    output_root = str(tmp_path / "run")
    crispresso_dir = tmp_path / "run.crispresso"
    crispresso_dir.mkdir()
    amp_dir = tmp_path / "run.seq_by_amplicon"
    amp_dir.mkdir()
    information = {}
    for amp in ("ampA", "ampB"):
        r1 = amp_dir / f"{amp}.r1.fq.gz"
        r2 = amp_dir / f"{amp}.r2.fq.gz"
        r1.write_text("reads\n")
        r2.write_text("reads\n")
        information[amp] = {
            "name": amp, "aln_count": "1", "reads_r1_file": str(r1),
            "reads_r2_file": str(r2), "amp_seqs": "ACGT", "guide_seq": "",
        }
    monkeypatch.setattr(crispresso, "tool_identity", _tool)

    class Pool:
        def __init__(self, *_args):
            pass

        def map_async(self, _function, jobs):
            values = []
            for job in jobs:
                os.makedirs(job["crispresso_run_folder"], exist_ok=True)
                with open(os.path.join(job["crispresso_run_folder"], "CRISPResso2_info.json"), "w") as handle:
                    handle.write("{}\n")
                if job["amplicon_name"] == "ampB":
                    with gzip.open(os.path.join(job["crispresso_run_folder"], "CRISPResso_output.fastq.gz"), "wt") as handle:
                        handle.write("fastq\n")
                with open(job["finished_file"], "w"):
                    pass
                values.append({"returncode": 0, "error": None, "command": job["command"]})
            return SimpleNamespace(get=lambda *_args: values)

        def close(self):
            pass

        def join(self):
            pass

    monkeypatch.setattr(crispresso.mp, "Pool", Pool)
    manager = CacheManager(output_root)
    with pytest.raises(Exception, match="Cannot commit cache record"):
        crispresso.run_crispresso_commands(
            ["ampA", "ampB"], information, output_root, str(crispresso_dir),
            True, 1, alleles=False, cache_manager=manager,
        )

    assert manager.load("crispresso_reads", "ampA") is None
    assert manager.load("crispresso_reads", "ampB") is not None


def _seed_parse_crispresso_cache(tmp_path):
    output_root = str(tmp_path / "run")
    amp_dir = tmp_path / "run.seq_by_amplicon"
    amp_dir.mkdir()
    folder = tmp_path / "CRISPResso_on_ampA"
    folder.mkdir()
    (folder / "CRISPResso_output.fastq.gz").write_text("fastq\n")
    summary_header = (
        "cell\tall_cell_read_count\tall_cell_mut_pct\tall_cell_allele_string\t"
        "final_cell_read_count\tfinal_cell_mut_allele_pct\tfinal_cell_allele_string\t"
        "final_num_refs_covered\tfinal_cell_allele_mod_string\t"
        "final_cell_allele_mod_types_string\tfinal_cell_allele_readcount_string\t"
        "final_ref_read_count_string\tfinal_ref_mut_allele_fracs_string\n"
    )
    (tmp_path / "CRISPResso_on_ampA.summ").write_text(
        summary_header + "cellA\t10\t0\tNA\t10\t0\tNA\t1\tU\tU\t10\t10\t0\n"
    )
    for suffix in (".summarize_indels.out", ".summarize_alleles.out"):
        (tmp_path / ("CRISPResso_on_ampA" + suffix)).write_text(
            "cell\tread_count\tmod_pct\tallele_0\ncellA\t10\t0\tNA\n"
        )
    (tmp_path / "CRISPResso_on_ampA.summ.finished").write_text(
        "Total reads\t10\nCRISPResso2 aligned reads\t10\nIgnore substitutions\tFalse\n"
    )
    (amp_dir / "03_alleles_all_cells.ampA.fq").write_text("alleles\n")
    amp_info = {
        "ampA": {
            "amp_seqs": "ACGT",
            "input_ref_allele_counts": "1",
        }
    }
    manager = CacheManager(output_root)
    record = crispresso._build_parse_crispresso_cache_record(
        manager, "ampA", amp_info["ampA"], str(folder), False, 5
    )
    manager.commit(
        record,
        crispresso._parse_crispresso_cache_requirements(output_root, "ampA", str(folder)),
    )
    information = {
        "ampA": {"status": "Completed", "crispresso_run_folder": str(folder)}
    }
    return manager, output_root, amp_info, information


def test_parse_crispresso_cache_hit_skips_per_amplicon_parser(tmp_path, monkeypatch):
    manager, output_root, amp_info, information = _seed_parse_crispresso_cache(tmp_path)
    monkeypatch.setattr(
        crispresso,
        "parse_one_crispresso_output",
        lambda _args: (_ for _ in ()).throw(AssertionError("cache miss")),
    )
    result = crispresso.parse_crispresso_outputs(
        ["ampA"], amp_info, str(tmp_path / "unused.tsv"), information,
        output_root, 0, 0, 1, cache_manager=manager,
    )
    assert result.loc["cellA", "totCount.ampA"] == 10
    assert manager.events[-1]["status"] == "hit"


def test_parse_crispresso_setting_change_invalidates(tmp_path, monkeypatch):
    manager, output_root, amp_info, information = _seed_parse_crispresso_cache(tmp_path)
    monkeypatch.setattr(
        crispresso,
        "parse_one_crispresso_output",
        lambda _args: (_ for _ in ()).throw(RuntimeError("reparse")),
    )
    with pytest.raises(RuntimeError, match="reparse"):
        crispresso.parse_crispresso_outputs(
            ["ampA"], amp_info, str(tmp_path / "unused.tsv"), information,
            output_root, 0, 0, 1, ignore_substitutions=True,
            cache_manager=manager,
        )
    assert manager.events[-1]["status"] == "invalid"


def _seed_filter_selected_cache(tmp_path):
    output_root = str(tmp_path / "run")
    amp_dir = tmp_path / "run.seq_by_amplicon"
    amp_dir.mkdir()
    paths = crispresso._filter_selected_paths(output_root, "ampA")
    for key in ("input_r1", "input_r2", "output_r1", "output_r2", "output_alleles"):
        with gzip.open(paths[key], "wt") as handle:
            handle.write("@read:cellA\nACGT\n+\nIIII\n")
    with open(paths["input_alleles"], "w") as handle:
        handle.write("@amp:cellA:1\nACGT\n+\nIIII\n")
    manager = CacheManager(output_root)
    barcode_hash = crispresso._barcode_set_sha256({"cellA"})
    record = crispresso._build_filter_selected_cache_record(
        manager, output_root, "ampA", barcode_hash
    )
    requirements = crispresso._filter_selected_cache_requirements(output_root, "ampA")
    manager.commit(record, requirements)
    return manager, output_root, paths


def test_selected_filter_cache_hit_skips_workers(tmp_path, monkeypatch):
    manager, output_root, _paths = _seed_filter_selected_cache(tmp_path)
    monkeypatch.setattr(
        crispresso.mp,
        "Pool",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("cache miss")),
    )
    parsed = pd.DataFrame({"Color": ["HQ_HI"]}, index=["cellA"])
    result = crispresso.filter_amplicon_reads(
        output_root, parsed, ["ampA"], ["HQ_HI"], 1, cache_manager=manager
    )
    assert not result["read_filter_failures"]
    assert not result["allele_filter_failures"]
    assert manager.events[-1]["status"] == "hit"


def test_removed_amplicon_prunes_only_trusted_stage_owned_outputs(tmp_path):
    output_root = str(tmp_path / "run")
    crispresso_dir = tmp_path / "run.crispresso"
    crispresso_dir.mkdir()
    run_folder = crispresso_dir / "CRISPResso_on_removed"
    run_folder.mkdir()
    (run_folder / "CRISPResso2_info.json").write_text("{}\n")
    with gzip.open(run_folder / "CRISPResso_output.fastq.gz", "wt") as handle:
        handle.write("fastq\n")
    marker = crispresso_dir / "removed.finished"
    marker.write_text("")
    log = crispresso_dir / "removed.log"
    log.write_text("old\n")
    manager = CacheManager(output_root)
    record = manager.new_record(
        "crispresso_reads", "removed", algorithm_version=1
    )
    manager.commit(
        record,
        crispresso._crispresso_cache_requirements(str(marker), str(run_folder)),
    )

    unsafe = manager.new_record(
        "crispresso_reads", "../outside", algorithm_version=1
    )
    manager.commit(unsafe, ())
    outside = tmp_path / "outside.finished"
    outside.write_text("keep\n")

    removed = crispresso.prune_removed_amplicon_caches(
        manager, [], output_root, str(crispresso_dir)
    )

    assert ("crispresso_reads", "removed") in removed
    assert not run_folder.exists()
    assert not marker.exists()
    assert not log.exists()
    assert not os.path.exists(manager.record_path("crispresso_reads", "removed"))
    assert outside.exists()
    assert os.path.exists(manager.record_path("crispresso_reads", "../outside"))


def test_replaced_split_prunes_fastqs_owned_by_prior_record(tmp_path):
    amp_dir = tmp_path / "run.seq_by_amplicon"
    amp_dir.mkdir()
    output_root = str(tmp_path / "run")
    old_r1 = amp_dir / "03_reads_all_cells.old.r1.fq.gz"
    old_r2 = amp_dir / "03_reads_all_cells.old.r2.fq.gz"
    for path in (old_r1, old_r2):
        with gzip.open(path, "wt") as handle:
            handle.write("")
    manager = CacheManager(output_root)
    previous = manager.new_record("split_reads", algorithm_version=1)
    manager.commit(
        previous,
        (
            OutputRequirement("reads:old:r1", str(old_r1), strategy="stat", allow_empty=True, validator="gzip"),
            OutputRequirement("reads:old:r2", str(old_r2), strategy="stat", allow_empty=True, validator="gzip"),
        ),
    )

    removed = amplicon_assignment._prune_obsolete_split_outputs(
        manager, previous, (), str(amp_dir)
    )

    assert set(removed) == {str(old_r1), str(old_r2)}
    assert not old_r1.exists()
    assert not old_r2.exists()


def test_selected_filter_source_change_invalidates_with_same_barcodes(tmp_path):
    manager, output_root, paths = _seed_filter_selected_cache(tmp_path)
    with gzip.open(paths["input_r1"], "at") as handle:
        handle.write("@changed:cellA\nTGCA\n+\nIIII\n")
    record = crispresso._build_filter_selected_cache_record(
        manager,
        output_root,
        "ampA",
        crispresso._barcode_set_sha256({"cellA"}),
    )
    decision = manager.evaluate(
        record, crispresso._filter_selected_cache_requirements(output_root, "ampA")
    )
    assert decision.status == "invalid"
