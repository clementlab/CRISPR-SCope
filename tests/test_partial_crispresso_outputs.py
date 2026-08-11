import gzip

import pandas as pd
import pytest

import CRISPRSCope.cli as cli
from CRISPRSCope.h5ad.loaders import _parse_and_write_parquet
from CRISPRSCope.cli import (
    generate_amplicon_score,
    parse_crispresso_outputs,
    run_crispresso_commands,
    write_filtered_editing_summary_from_filtered_crispresso,
)


def test_parse_crispresso_outputs_preserves_failed_amplicons_as_na(tmp_path):
    output_root = str(tmp_path / "run")
    amplicon_names = ["amp_ok", "amp_failed"]

    completed_run_folder = tmp_path / "CRISPResso_on_amp_ok"
    completed_run_folder.mkdir()
    summ_path = completed_run_folder.with_suffix(completed_run_folder.suffix + ".summ")
    summ_path.write_text(
        "cell	all_cell_read_count	all_cell_mut_pct	all_cell_allele_string	final_cell_read_count	final_cell_mut_allele_pct	final_cell_allele_string	final_num_refs_covered	final_cell_allele_mod_string	final_cell_allele_mod_types_string	final_cell_allele_readcount_string	final_ref_read_count_string	final_ref_mut_allele_fracs_string\n"
        "cellA	12	25.0	NA	10	20.0	NA	1	M	D	10	10	20.0\n"
    )
    finished_marker = completed_run_folder.with_suffix(completed_run_folder.suffix + ".summ.finished")
    finished_marker.write_text("Ignore substitutions\tFalse\n")

    amplicon_information = {
        "amp_ok": {"input_ref_allele_counts": "1"},
        "amp_failed": {"input_ref_allele_counts": "1"},
    }
    crispresso_information = {
        "amp_ok": {
            "status": "Completed",
            "crispresso_run_folder": str(completed_run_folder),
        },
        "amp_failed": {
            "status": "Failed",
            "crispresso_run_folder": str(tmp_path / "CRISPResso_on_amp_failed"),
        },
    }

    summary_df = parse_crispresso_outputs(
        amplicon_names=amplicon_names,
        amplicon_information=amplicon_information,
        amplicon_info_file=str(tmp_path / "amplicons.txt"),
        crispresso_information=crispresso_information,
        output_root=output_root,
        min_total_reads_per_barcode=0,
        min_reads_per_amplicon_per_cell=0,
        n_processes=1,
    )

    assert list(summary_df.index) == ["cellA"]
    assert summary_df.loc["cellA", "totCount.amp_ok"] == 12
    assert summary_df.loc["cellA", "totCount.amp_failed"] == "NA"

    editing_summary = pd.read_csv(f"{output_root}.editingSummary.txt", sep="\t")
    assert list(editing_summary.columns) == [
        "cell",
        "totCount.amp_ok",
        "modPct.amp_ok",
        "totCount.amp_failed",
        "modPct.amp_failed",
    ]
    assert editing_summary.loc[0, "totCount.amp_ok"] == 10
    assert pd.isna(editing_summary.loc[0, "totCount.amp_failed"])
    assert pd.isna(editing_summary.loc[0, "modPct.amp_failed"])

    amp_score = pd.read_csv(f"{output_root}.amplicon_score.txt", sep="\t", index_col=0)
    assert list(amp_score.index) == ["cellA"]
    assert amp_score.loc["cellA", "Read Count"] == 12
    assert not (tmp_path / "run.filteredEditingSummary.txt").exists()
    assert (tmp_path / "run.filteredEditingSummaryPseudobulk.txt").exists()

    filtered_pseudobulk = pd.read_csv(f"{output_root}.filteredEditingSummaryPseudobulk.txt", sep="\t")
    assert list(filtered_pseudobulk["cell"]) == ["cellA"]
    assert filtered_pseudobulk.loc[0, "totCount.amp_ok"] == 12
    assert filtered_pseudobulk.loc[0, "modPct.amp_ok"] == 25.0
    assert pd.isna(filtered_pseudobulk.loc[0, "totCount.amp_failed"])
    assert pd.isna(filtered_pseudobulk.loc[0, "modPct.amp_failed"])


def test_generate_amplicon_score_ignores_all_na_failed_amplicons():
    totals = pd.DataFrame(
        {
            "totCount.amp_ok": [10, 20],
            "totCount.amp_failed": ["NA", "NA"],
        },
        index=["cell_low", "cell_high"],
    )

    amp_score = generate_amplicon_score(
        totals,
        min_reads_per_amplicon_per_cell=0,
        min_total_reads_per_barcode=0,
    )

    assert list(amp_score.index) == ["cell_high", "cell_low"]
    assert amp_score.loc["cell_high", "Read Count"] == 20
    assert amp_score.loc["cell_low", "Read Count"] == 10
    assert amp_score.loc["cell_high", "Amplicon Score"] == 161
    assert amp_score.loc["cell_low", "Amplicon Score"] == 0


def test_generate_amplicon_score_rejects_duplicate_amplicon_columns():
    totals = pd.DataFrame(
        [[10, 12], [20, 22]],
        columns=["totCount.amp_dup", "totCount.amp_dup"],
        index=["cellA", "cellB"],
    )

    with pytest.raises(ValueError, match="Duplicate amplicon names detected"):
        generate_amplicon_score(
            totals,
            min_reads_per_amplicon_per_cell=0,
            min_total_reads_per_barcode=0,
        )


def _write_first_pass_summ(path, rows):
    header = (
        "cell\tall_cell_read_count\tall_cell_mut_pct\tall_cell_allele_string\t"
        "final_cell_read_count\tfinal_cell_mut_allele_pct\tfinal_cell_allele_string\t"
        "final_num_refs_covered\tfinal_cell_allele_mod_string\tfinal_cell_allele_mod_types_string\t"
        "final_cell_allele_readcount_string\tfinal_ref_read_count_string\tfinal_ref_mut_allele_fracs_string\n"
    )
    path.write_text(header + "".join(rows))


def _write_filtered_crispresso_fastq(path, records):
    with gzip.open(path, "wt") as fout:
        for header, sequence, annotation in records:
            fout.write(header + "\n")
            fout.write(sequence + "\n")
            fout.write(annotation + "\n")
            fout.write("a" * len(sequence) + "\n")


def _write_plain_filtered_crispresso_fastq(path, records):
    with open(path, "w") as fout:
        for header, sequence, annotation in records:
            fout.write(header + "\n")
            fout.write(sequence + "\n")
            fout.write(annotation + "\n")
            fout.write("a" * len(sequence) + "\n")


def _make_crispresso_inputs(tmp_path, amplicon_name="ampA"):
    r1 = tmp_path / f"{amplicon_name}.r1.fq.gz"
    r2 = tmp_path / f"{amplicon_name}.r2.fq.gz"
    r1.write_text("")
    r2.write_text("")
    return {
        amplicon_name: {
            "aln_count": "1",
            "reads_r1_file": str(r1),
            "reads_r2_file": str(r2),
            "amp_seqs": "ACGT",
            "guide_seq": "AC",
        }
    }


def _write_completed_crispresso_output(crispresso_dir, amplicon_name="ampA"):
    run_folder = crispresso_dir / f"CRISPResso_on_{amplicon_name}"
    run_folder.mkdir(parents=True, exist_ok=True)
    (run_folder / "CRISPResso2_info.json").write_text("{}\n")
    finished_file = crispresso_dir / f"{amplicon_name}.finished"
    finished_file.write_text("done\n")
    return run_folder, finished_file


def test_run_crispresso_commands_does_not_pass_ignore_substitutions_to_crispresso(tmp_path):
    output_root = str(tmp_path / "run")
    crispresso_dir = tmp_path / "run.crispresso"
    crispresso_dir.mkdir()
    amplicon_information = _make_crispresso_inputs(tmp_path)
    _write_completed_crispresso_output(crispresso_dir)

    result = run_crispresso_commands(
        amplicon_names=["ampA"],
        amplicon_information=amplicon_information,
        output_root=output_root,
        crispresso_dir=str(crispresso_dir),
        suppress_sub_crispresso_plots=False,
        n_processes=1,
        alleles=False,
    )

    assert "--ignore_substitutions" not in result["ampA"]["crispresso_command"]


def test_parse_crispresso_outputs_reparses_when_ignore_substitutions_changes(tmp_path, monkeypatch):
    output_root = str(tmp_path / "run")
    run_folder = tmp_path / "CRISPResso_on_ampA"
    run_folder.mkdir()
    summ_path = run_folder.with_suffix(run_folder.suffix + ".summ")
    _write_first_pass_summ(
        summ_path,
        [
            "cellA\t10\t100.0\tNA\t10\t100.0\tNA\t1\tM\tS\t10\t10\t100.0\n",
        ],
    )
    finished_marker = run_folder.with_suffix(run_folder.suffix + ".summ.finished")
    finished_marker.write_text("Total reads\t10\nCRISPResso2 aligned reads\t10\nIgnore substitutions\tFalse\n")

    calls = []

    def fake_parse_one_crispresso_output(this_args):
        calls.append(this_args)
        _write_first_pass_summ(
            summ_path,
            [
                "cellA\t10\t0.0\tNA\t10\t0.0\tNA\t1\tU\tU\t10\t10\t0.0\n",
            ],
        )
        finished_marker.write_text(
            "Total reads\t10\nCRISPResso2 aligned reads\t10\nIgnore substitutions\tTrue\n"
        )

    monkeypatch.setattr(cli, "parse_one_crispresso_output", fake_parse_one_crispresso_output)

    summary_df = parse_crispresso_outputs(
        amplicon_names=["ampA"],
        amplicon_information={"ampA": {"input_ref_allele_counts": "1"}},
        amplicon_info_file=str(tmp_path / "amplicons.txt"),
        crispresso_information={
            "ampA": {
                "status": "Completed",
                "crispresso_run_folder": str(run_folder),
            },
        },
        output_root=output_root,
        min_total_reads_per_barcode=0,
        min_reads_per_amplicon_per_cell=0,
        n_processes=1,
        ignore_substitutions=True,
    )

    assert len(calls) == 1
    assert calls[0]["ignore_substitutions"] is True
    assert summary_df.loc["cellA", "modPct.ampA"] == 0.0


def test_parse_crispresso_outputs_reparses_legacy_summary_without_ignore_substitutions_metadata(
    tmp_path, monkeypatch
):
    output_root = str(tmp_path / "run")
    run_folder = tmp_path / "CRISPResso_on_ampA"
    run_folder.mkdir()
    summ_path = run_folder.with_suffix(run_folder.suffix + ".summ")
    _write_first_pass_summ(
        summ_path,
        [
            "cellA\t10\t100.0\tNA\t10\t100.0\tNA\t1\tM\tS\t10\t10\t100.0\n",
        ],
    )
    finished_marker = run_folder.with_suffix(run_folder.suffix + ".summ.finished")
    finished_marker.write_text("Total reads\t10\nCRISPResso2 aligned reads\t10\n")

    calls = []

    def fake_parse_one_crispresso_output(this_args):
        calls.append(this_args)
        _write_first_pass_summ(
            summ_path,
            [
                "cellA\t10\t0.0\tNA\t10\t0.0\tNA\t1\tU\tU\t10\t10\t0.0\n",
            ],
        )
        finished_marker.write_text(
            "Total reads\t10\nCRISPResso2 aligned reads\t10\nIgnore substitutions\tFalse\n"
        )

    monkeypatch.setattr(cli, "parse_one_crispresso_output", fake_parse_one_crispresso_output)

    summary_df = parse_crispresso_outputs(
        amplicon_names=["ampA"],
        amplicon_information={"ampA": {"input_ref_allele_counts": "1"}},
        amplicon_info_file=str(tmp_path / "amplicons.txt"),
        crispresso_information={
            "ampA": {
                "status": "Completed",
                "crispresso_run_folder": str(run_folder),
            },
        },
        output_root=output_root,
        min_total_reads_per_barcode=0,
        min_reads_per_amplicon_per_cell=0,
        n_processes=1,
        ignore_substitutions=False,
    )

    assert len(calls) == 1
    assert calls[0]["ignore_substitutions"] is False
    assert summary_df.loc["cellA", "modPct.ampA"] == 0.0


def test_filtered_editing_summary_uses_filtered_crispresso_calls_and_original_support(tmp_path):
    output_root = str(tmp_path / "run")
    pd.DataFrame({"Color": ["HQ_HI", "HQ_HI", "HQ_HI"]}, index=["cell_wt", "cell_het", "cell_mut"]).to_csv(
        f"{output_root}.amplicon_score.txt",
        sep="\t",
    )
    pseudobulk_path = tmp_path / "run.filteredEditingSummaryPseudobulk.txt"
    pseudobulk_path.write_text("sentinel\n")

    first_pass_folder = tmp_path / "CRISPResso_on_ampA"
    first_pass_folder.mkdir()
    _write_first_pass_summ(
        tmp_path / "CRISPResso_on_ampA.summ",
        [
            "cell_wt\t100\t0\tNA\t80\t0\tNA\t1\tU,U\tU,U\t40,40\t80\t0\n",
            "cell_het\t100\t50\tNA\t80\t50\tNA\t1\tU,M\tU,D\t40,40\t80\t50\n",
            "cell_mut\t100\t100\tNA\t80\t100\tNA\t1\tM,M\tD,D\t40,40\t80\t100\n",
        ],
    )

    filtered_folder = tmp_path / "CRISPResso_filtered_on_ampA"
    filtered_folder.mkdir()
    _write_filtered_crispresso_fastq(
        filtered_folder / "CRISPResso_output.fastq.gz",
        [
            ("@ampA:Reference:DEL=_INS=_SUB=:cell_wt:1", "ACGT", "+ ALN=Reference DEL= INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
            ("@ampA:Reference:DEL=_INS=_SUB=:cell_wt:2", "ACGT", "+ ALN=Reference DEL= INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
            ("@ampA:Reference:DEL=_INS=_SUB=:cell_het:1", "ACGT", "+ ALN=Reference DEL= INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
            ("@ampA:Reference:DEL=12_INS=_SUB=:cell_het:2", "ACGT", "+ ALN=Reference DEL=12 INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
            ("@ampA:Reference:DEL=12_INS=_SUB=:cell_mut:1", "ACGT", "+ ALN=Reference DEL=12 INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
            ("@ampA:Reference:DEL=12_INS=_SUB=:cell_mut:2", "ACGT", "+ ALN=Reference DEL=12 INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
        ],
    )

    result = write_filtered_editing_summary_from_filtered_crispresso(
        amplicon_names=["ampA"],
        crispresso_information={"ampA": {"status": "Completed", "crispresso_run_folder": str(first_pass_folder)}},
        crispresso_filtered_information={"ampA": {"status": "Completed", "crispresso_run_folder": str(filtered_folder)}},
        output_root=output_root,
        ignore_substitutions=False,
    )

    filtered_summary = pd.read_csv(f"{output_root}.filteredEditingSummary.txt", sep="\t", index_col=0)
    assert filtered_summary.loc["cell_wt", "totCount.ampA"] == 80
    assert filtered_summary.loc["cell_wt", "modPct.ampA"] == 0
    assert filtered_summary.loc["cell_het", "totCount.ampA"] == 80
    assert filtered_summary.loc["cell_het", "modPct.ampA"] == 50
    assert filtered_summary.loc["cell_mut", "totCount.ampA"] == 80
    assert filtered_summary.loc["cell_mut", "modPct.ampA"] == 100
    assert result.loc["cell_het", "Color"] == "HQ_HI"
    assert pseudobulk_path.read_text() == "sentinel\n"


def test_filtered_editing_summary_ignores_stale_modified_header_when_filtered_call_is_unmodified(tmp_path):
    output_root = str(tmp_path / "run")
    pd.DataFrame({"Color": ["HQ_HI"]}, index=["cellA"]).to_csv(f"{output_root}.amplicon_score.txt", sep="\t")

    first_pass_folder = tmp_path / "CRISPResso_on_ampA"
    first_pass_folder.mkdir()
    _write_first_pass_summ(
        tmp_path / "CRISPResso_on_ampA.summ",
        ["cellA\t100\t50\tNA\t80\t50\tNA\t1\tU,M\tU,S\t40,40\t80\t50\n"],
    )

    filtered_folder = tmp_path / "CRISPResso_filtered_on_ampA"
    filtered_folder.mkdir()
    _write_filtered_crispresso_fastq(
        filtered_folder / "CRISPResso_output.fastq.gz",
        [
            ("@ampA:Reference:DEL=_INS=_SUB=:cellA:1", "ACGT", "+ ALN=Reference DEL= INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
            ("@ampA:Reference:DEL=_INS=_SUB=120:cellA:2", "ACGT", "+ ALN=Reference DEL= INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
        ],
    )

    write_filtered_editing_summary_from_filtered_crispresso(
        amplicon_names=["ampA"],
        crispresso_information={"ampA": {"status": "Completed", "crispresso_run_folder": str(first_pass_folder)}},
        crispresso_filtered_information={"ampA": {"status": "Completed", "crispresso_run_folder": str(filtered_folder)}},
        output_root=output_root,
        ignore_substitutions=False,
    )

    filtered_summary = pd.read_csv(f"{output_root}.filteredEditingSummary.txt", sep="\t", index_col=0)
    assert filtered_summary.loc["cellA", "totCount.ampA"] == 80
    assert filtered_summary.loc["cellA", "modPct.ampA"] == 0


def test_filtered_editing_summary_respects_ignore_substitutions(tmp_path):
    output_root = str(tmp_path / "run")
    pd.DataFrame({"Color": ["HQ_HI"]}, index=["cellA"]).to_csv(f"{output_root}.amplicon_score.txt", sep="\t")

    first_pass_folder = tmp_path / "CRISPResso_on_ampA"
    first_pass_folder.mkdir()
    _write_first_pass_summ(
        tmp_path / "CRISPResso_on_ampA.summ",
        ["cellA\t100\t50\tNA\t80\t50\tNA\t1\tU,M\tU,S\t40,40\t80\t50\n"],
    )

    filtered_folder = tmp_path / "CRISPResso_filtered_on_ampA"
    filtered_folder.mkdir()
    _write_filtered_crispresso_fastq(
        filtered_folder / "CRISPResso_output.fastq.gz",
        [
            ("@ampA:Reference:DEL=_INS=_SUB=:cellA:1", "ACGT", "+ ALN=Reference DEL= INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
            ("@ampA:Reference:DEL=_INS=_SUB=120:cellA:2", "ACGT", "+ ALN=Reference DEL= INS= SUB=120 ALN_REF=ACGT ALN_SEQ=ACGT"),
        ],
    )

    write_filtered_editing_summary_from_filtered_crispresso(
        amplicon_names=["ampA"],
        crispresso_information={"ampA": {"status": "Completed", "crispresso_run_folder": str(first_pass_folder)}},
        crispresso_filtered_information={"ampA": {"status": "Completed", "crispresso_run_folder": str(filtered_folder)}},
        output_root=output_root,
        ignore_substitutions=True,
    )

    filtered_summary = pd.read_csv(f"{output_root}.filteredEditingSummary.txt", sep="\t", index_col=0)
    assert filtered_summary.loc["cellA", "totCount.ampA"] == 80
    assert filtered_summary.loc["cellA", "modPct.ampA"] == 0


def test_filtered_editing_summary_accepts_plain_text_crispresso_output_with_gz_suffix(tmp_path):
    output_root = str(tmp_path / "run")
    pd.DataFrame({"Color": ["HQ_HI"]}, index=["cellA"]).to_csv(f"{output_root}.amplicon_score.txt", sep="\t")

    first_pass_folder = tmp_path / "CRISPResso_on_ampA"
    first_pass_folder.mkdir()
    _write_first_pass_summ(
        tmp_path / "CRISPResso_on_ampA.summ",
        ["cellA\t100\t50\tNA\t80\t50\tNA\t1\tU,M\tU,D\t40,40\t80\t50\n"],
    )

    filtered_folder = tmp_path / "CRISPResso_filtered_on_ampA"
    filtered_folder.mkdir()
    _write_plain_filtered_crispresso_fastq(
        filtered_folder / "CRISPResso_output.fastq.gz",
        [
            ("@ampA:Reference:DEL=_INS=_SUB=:cellA:1", "ACGT", "+ ALN=Reference DEL= INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
            ("@ampA:Reference:DEL=12_INS=_SUB=:cellA:2", "ACGT", "+ ALN=Reference DEL=12 INS= SUB= ALN_REF=ACGT ALN_SEQ=ACGT"),
        ],
    )

    write_filtered_editing_summary_from_filtered_crispresso(
        amplicon_names=["ampA"],
        crispresso_information={"ampA": {"status": "Completed", "crispresso_run_folder": str(first_pass_folder)}},
        crispresso_filtered_information={"ampA": {"status": "Completed", "crispresso_run_folder": str(filtered_folder)}},
        output_root=output_root,
        ignore_substitutions=False,
    )

    filtered_summary = pd.read_csv(f"{output_root}.filteredEditingSummary.txt", sep="\t", index_col=0)
    assert filtered_summary.loc["cellA", "totCount.ampA"] == 80
    assert filtered_summary.loc["cellA", "modPct.ampA"] == 50


def test_h5ad_loader_accepts_plain_text_crispresso_output_with_gz_suffix(tmp_path):
    input_fastq = tmp_path / "CRISPResso_on_ampA" / "CRISPResso_output.fastq.gz"
    input_fastq.parent.mkdir()
    _write_plain_filtered_crispresso_fastq(
        input_fastq,
        [
            ("@ampA:Reference:DEL=_INS=_SUB=:cellA:1", "ACGT", "+ ALN=Reference"),
            ("@ampA:Reference:DEL=_INS=_SUB=:cellA:2", "ACGT", "+ ALN=Reference"),
            ("@ampA:Reference:DEL=_INS=_SUB=:cellB:1", "TGCA", "+ ALN=Reference"),
        ],
    )

    output_parquet = tmp_path / "ampA.parquet"
    result = _parse_and_write_parquet(input_fastq, output_parquet)

    assert result == str(output_parquet)
    parsed = pd.read_parquet(output_parquet).sort_values(
        ["cell_barcode", "allele_sequence"]
    ).reset_index(drop=True)
    assert parsed.to_dict("records") == [
        {
            "cell_barcode": "cellA",
            "amplicon_name": "ampA",
            "allele_sequence": "ACGT",
            "count": 2,
        },
        {
            "cell_barcode": "cellB",
            "amplicon_name": "ampA",
            "allele_sequence": "TGCA",
            "count": 1,
        },
    ]
