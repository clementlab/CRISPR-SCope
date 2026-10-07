"""Synthetic read tests for the first-pass genotype policy."""

import gzip
import io
import json

import numpy as np
import pandas as pd
import pytest

from CRISPRSCope.allele_calling import local_quality, parse_allele_support
from CRISPRSCope.crispresso import (
    _parse_cache_matches_ignore_substitutions,
    parse_crispresso_outputs,
    parse_one_crispresso_output,
    write_max_alleles,
    write_filtered_editing_summary_from_filtered_crispresso,
)
from CRISPRSCope.h5ad.builder import CRISPRSCopeAnnDataBuilder
from CRISPRSCope.settings import _parse_allele_calling_config


WT = "A" * 30 + "C" * 10 + "T" * 30


def _mutant(position):
    return WT[:position] + "G" + WT[position + 1:]


def _record(cell, kind, index, *, reference="Reference", poor=False, missing_alignment=False):
    position = {"wt": None, "mut1": 35, "mut2": 36}.get(kind)
    if position is None and kind != "wt":
        position = int(kind[3:])
    sequence = WT if position is None else _mutant(position)
    quality = list("I" * len(sequence))
    if poor and position is not None:
        quality[position] = "#"
    annotation = (
        f"+ ALN={reference} DEL= INS= SUB={'' if position is None else position} "
        f"ALN_REF={WT} ALN_SEQ={sequence}"
    )
    if missing_alignment:
        annotation = f"+ ALN={reference} DEL= INS= SUB={'' if position is None else position}"
    return f"@r{index}:{cell}\n{sequence}\n{annotation}\n{''.join(quality)}\n"


def _run_parser(tmp_path, reads, *, depth=8, support=2, ref_counts="2"):
    run = tmp_path / "CRISPResso_on_ampA"
    run.mkdir()
    with gzip.open(run / "CRISPResso_output.fastq.gz", "wt") as handle:
        handle.write("".join(reads))
    info = tmp_path / "amplicon_info.tsv"
    info.write_text("name\tamp_seqs\nampA\t" + WT + "\n")
    (tmp_path / "run.seq_by_amplicon").mkdir()
    parse_one_crispresso_output({
        "amplicon_name": "ampA",
        "amplicon_info_file": str(info),
        "crispresso_run_folder": str(run),
        "input_ref_allele_counts": ref_counts,
        "min_num_reads_per_cell": 0,
        "min_allele_support": support,
        "min_reads_per_amplicon_for_genotype": depth,
        "ignore_substitutions": False,
        "output_root": str(tmp_path / "run"),
    })
    summary = pd.read_csv(str(run) + ".summ", sep="\t").set_index("cell")
    allele_fastq = tmp_path / "run.seq_by_amplicon" / "03_alleles_all_cells.ampA.fq"
    return summary, run, info, allele_fastq


@pytest.mark.parametrize("read_count,expected_status", [(7, "low_depth"), (8, "called")])
def test_depth_gate_keeps_raw_evidence_and_masks_genotype(tmp_path, read_count, expected_status):
    reads = [_record("cellA", "wt", i) for i in range(read_count - 1)]
    reads.append(_record("cellA", "mut1", read_count))
    summary, run, _, allele_fastq = _run_parser(tmp_path, reads)
    row = summary.loc["cellA"]
    assert row["all_cell_read_count"] == read_count
    assert row["all_cell_mut_pct"] == round(100 / read_count, 2)
    assert row["call_status"] == expected_status
    assert pd.isna(row["final_cell_mut_allele_pct"]) if read_count == 7 else row["final_cell_mut_allele_pct"] == 0
    if read_count == 7:
        assert pd.isna(row["final_cell_allele_string"])
    assert (allele_fastq.read_text() == "") if read_count == 7 else "@ampA:" in allele_fastq.read_text()
    assert "Minimum genotype reads\t8" in (tmp_path / "CRISPResso_on_ampA.summ.finished").read_text()


def test_single_qualifying_allele_fills_expected_copies(tmp_path):
    reads = [_record("cellA", "wt", i) for i in range(7)] + [_record("cellA", "mut1", 7)]
    summary, _, _, _ = _run_parser(tmp_path, reads)
    row = summary.loc["cellA"]
    assert row["call_status"] == "called"
    assert row["final_cell_allele_string"].split(",") == ["Reference:DEL= INS= SUB="] * 2
    assert [item["reads"] for item in json.loads(row["competing_alleles_json"])] == [1]


@pytest.mark.parametrize("support,expected", [(2, 50.0), (0.1, 50.0), (0.26, 0.0), (3, 0.0)])
def test_count_and_fraction_support_are_inclusive(tmp_path, support, expected):
    reads = [_record("cellA", "wt", i) for i in range(6)]
    reads += [_record("cellA", "mut1", i + 6) for i in range(2)]
    summary, _, _, _ = _run_parser(tmp_path, reads, support=support)
    assert summary.loc["cellA", "final_cell_mut_allele_pct"] == expected


@pytest.mark.parametrize("poor,missing,reason", [
    (True, False, "quality"),
    (False, False, "stable_key"),
    (False, True, "stable_key"),
])
def test_boundary_count_tie_uses_local_quality_then_stable_key(tmp_path, poor, missing, reason):
    reads = [_record("cellA", "wt", i) for i in range(4)]
    reads += [_record("cellA", "mut2", i + 4, poor=poor, missing_alignment=missing) for i in range(2)]
    reads += [_record("cellA", "mut1", i + 6, missing_alignment=missing) for i in range(2)]
    summary, run, info, _ = _run_parser(tmp_path, reads)
    row = summary.loc["cellA"]
    assert row["final_cell_mut_allele_pct"] == 50.0
    assert row["count_tie"] == 1
    assert row["tie_resolution"] == reason
    selected = json.loads(row["selected_alleles_json"])
    assert any("SUB=35" in item["allele"] for item in selected)
    assert not any("SUB=36" in item["allele"] for item in selected)
    if poor and not missing:
        parse_crispresso_outputs(
            amplicon_names=["ampA"],
            amplicon_information={"ampA": {"amp_seqs": WT, "input_ref_allele_counts": "2"}},
            amplicon_info_file=str(info),
            crispresso_information={"ampA": {"status": "Completed", "crispresso_run_folder": str(run)}},
            output_root=str(tmp_path / "run"),
            min_total_reads_per_barcode=0,
            min_reads_per_amplicon_per_cell=0,
            min_num_reads_per_cell=0,
            n_processes=1,
        )
        qc = pd.read_csv(tmp_path / "run.alleleCallQC.txt", sep="\t").set_index("cell")
        assert qc.loc["cellA", "count_tie"] == 1
        assert qc.loc["cellA", "tie_resolution"] == "quality"


def test_equal_quality_boundary_tie_prefers_wildtype(tmp_path):
    reads = [_record("cellA", "mut1", i) for i in range(4)]
    reads += [_record("cellA", "wt", i + 4) for i in range(4)]
    summary, _, _, _ = _run_parser(tmp_path, reads, ref_counts="1")
    row = summary.loc["cellA"]
    assert row["count_tie"] == 1
    assert row["tie_resolution"] == "wildtype"
    assert row["final_cell_mut_allele_pct"] == 0.0


def test_no_qualifying_allele_is_no_call(tmp_path):
    reads = [_record("cellA", f"mut{i + 30}", i) for i in range(8)]
    summary, _, _, allele_fastq = _run_parser(tmp_path, reads)
    row = summary.loc["cellA"]
    assert row["call_status"] == "no_supported_allele"
    assert pd.isna(row["final_cell_mut_allele_pct"])
    assert allele_fastq.read_text() == ""


def test_missing_second_reference_withholds_multireference_genotype(tmp_path):
    reads = [_record("cellA", "wt", i) for i in range(8)]
    summary, _, _, _ = _run_parser(tmp_path, reads, ref_counts="1,1")
    assert summary.loc["cellA", "call_status"] == "no_supported_allele"
    assert pd.isna(summary.loc["cellA", "final_cell_mut_allele_pct"])


def test_two_covered_references_are_called_independently(tmp_path):
    reads = [_record("cellA", "wt", i) for i in range(4)]
    reads += [_record("cellA", "mut1", i + 4, reference="Amplicon1") for i in range(4)]
    summary, _, _, _ = _run_parser(tmp_path, reads, ref_counts="1,1")
    row = summary.loc["cellA"]
    assert row["call_status"] == "called"
    assert row["final_cell_mut_allele_pct"] == 50.0
    assert row["final_num_refs_covered"] == 2


def test_aggregation_qc_and_filtered_summary_preserve_low_depth_no_call(tmp_path):
    reads = [_record("cellA", "wt", i) for i in range(6)]
    reads.append(_record("cellA", "mut1", 6))
    _, run, info, _ = _run_parser(tmp_path, reads)
    output_root = str(tmp_path / "run")
    first_pass = {"ampA": {"status": "Completed", "crispresso_run_folder": str(run)}}
    parse_crispresso_outputs(
        amplicon_names=["ampA"],
        amplicon_information={"ampA": {"amp_seqs": WT, "input_ref_allele_counts": "2"}},
        amplicon_info_file=str(info),
        crispresso_information=first_pass,
        output_root=output_root,
        min_total_reads_per_barcode=0,
        min_reads_per_amplicon_per_cell=0,
        min_num_reads_per_cell=0,
        n_processes=1,
    )
    editing = pd.read_csv(output_root + ".editingSummary.txt", sep="\t").set_index("cell")
    pseudobulk = pd.read_csv(output_root + ".editingSummaryPseudobulk.txt", sep="\t").set_index("cell")
    qc = pd.read_csv(output_root + ".alleleCallQC.txt", sep="\t").set_index("cell")
    assert pd.isna(editing.loc["cellA", "modPct.ampA"])
    assert editing.loc["cellA", "totCount.ampA"] == 6
    assert pseudobulk.loc["cellA", "totCount.ampA"] == 7
    assert pseudobulk.loc["cellA", "modPct.ampA"] == 14.29
    assert qc.loc["cellA", "accepted_reads"] == 7
    assert qc.loc["cellA", "genotype_withheld"] == 1
    assert qc.loc["cellA", "call_status"] == "low_depth"

    filtered_run = tmp_path / "CRISPResso_filtered_on_ampA"
    filtered_run.mkdir()
    with gzip.open(filtered_run / "CRISPResso_output.fastq.gz", "wt") as handle:
        handle.write(_record("cellA:1", "wt", 0))
    write_filtered_editing_summary_from_filtered_crispresso(
        ["ampA"], first_pass,
        {"ampA": {"status": "Completed", "crispresso_run_folder": str(filtered_run)}},
        output_root,
    )
    filtered = pd.read_csv(output_root + ".filteredEditingSummary.txt", sep="\t").set_index("cell")
    assert pd.isna(filtered.loc["cellA", "modPct.ampA"])
    assert filtered.loc["cellA", "totCount.ampA"] == 6


def test_settings_migration_and_cache_marker(tmp_path):
    settings = tmp_path / "settings.txt"
    settings.write_text("min_reads_per_amplicon_for_genotype\t6\nmin_allele_support\t0.1\n")
    assert _parse_allele_calling_config(str(settings)) == (6, 0.1)
    settings.write_text("min_allele_count_cutoff\t2\n")
    with pytest.raises(ValueError, match="retired"):
        _parse_allele_calling_config(str(settings))
    marker = tmp_path / "parse.finished"
    marker.write_text("Ignore substitutions\tFalse\nAllele calling version\t2\nMinimum genotype reads\t8\nMinimum allele support\tcount:2\n")
    assert _parse_cache_matches_ignore_substitutions(str(marker), False, 8, 2)
    assert not _parse_cache_matches_ignore_substitutions(str(marker), False, 6, 2)
    assert not _parse_cache_matches_ignore_substitutions(str(marker), False, 8, 0.1)
    assert parse_allele_support("2").mode == "count"
    assert parse_allele_support("0.1").mode == "fraction"
    assert parse_allele_support(2).allows(2, 20)
    assert not parse_allele_support(2).allows(1, 20)
    assert parse_allele_support(0.1).allows(2, 20)
    assert not parse_allele_support(0.1).allows(2, 21)


def test_insertion_and_deletion_qualities_use_edit_flanks():
    insertion = local_quality(
        "ACTGT", "II#II",
        {"ALN_REF": "AC-GT", "ALN_SEQ": "ACTGT", "DEL": "", "INS": "1(1+T)", "SUB": ""},
        False,
    )
    deletion = local_quality(
        "AGT", "I#I",
        {"ALN_REF": "ACGT", "ALN_SEQ": "A-GT", "DEL": "1(1)", "INS": "", "SUB": ""},
        False,
    )
    assert insertion.event_score == 2
    assert deletion.event_score == 2


def test_local_quality_ignores_outside_window_mismatches_and_rejects_unmapped_edits():
    sequence = "G" + WT[1:35] + "G" + WT[36:]
    fields = {"ALN_REF": WT, "ALN_SEQ": sequence, "DEL": "", "INS": "", "SUB": "35"}
    quality = "#" + "I" * (len(sequence) - 1)
    assert local_quality(sequence, quality, fields, False).event_score == 40
    fields["SUB"] = "34"
    assert local_quality(sequence, quality, fields, False) is None


def test_tied_sequence_consensus_is_observed_not_synthetic_wildtype():
    key = "Reference:DEL= INS= SUB=35"
    mutant_g = _mutant(35)
    mutant_t = WT[:35] + "T" + WT[36:]
    output = io.StringIO()
    write_max_alleles(
        {"cellA": {key: {mutant_g: 1, mutant_t: 1}}},
        "cellA", [key], "ampA", "unused", output, WT,
        {(key, mutant_g): 10.0, (key, mutant_t): 30.0},
    )
    assert output.getvalue().splitlines()[1] == mutant_t
    assert WT not in output.getvalue()


def test_h5ad_masks_allele_layers_for_no_call(monkeypatch):
    monkeypatch.setattr(pd, "read_parquet", lambda _: pd.DataFrame({
        "cell_barcode": ["cellA"], "amplicon_name": ["ampA"],
        "count": [7], "allele_sequence": [WT],
    }))
    builder = CRISPRSCopeAnnDataBuilder(
        config={"analysis_parameters": {"zygosity": {
            "wt_max_mod_pct": 20, "het_max_mod_pct": 80,
            "hom_min_mod_pct": 80, "compound_het_min_allele2_pct": 20,
        }}},
        settings={}, amplicons=pd.DataFrame({"sequence": [WT]}, index=["ampA"]),
        editing_summary=pd.DataFrame({"totCount.ampA": [7], "modPct.ampA": [np.nan]}, index=["cellA"]),
        quality_scores=pd.DataFrame({"Color": ["HQ_HI"]}, index=["cellA"]),
        allele_parquet_paths=["synthetic.parquet"],
    )
    adata = builder.build()
    assert adata.layers["counts"][0, 0] == 7
    assert adata.layers["zygosity"][0, 0] == -1
    assert adata.layers["allele_seq_1"][0, 0] == b""
