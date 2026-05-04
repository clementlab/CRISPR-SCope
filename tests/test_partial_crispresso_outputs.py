import pandas as pd

from CRISPRSCope.cli import parse_crispresso_outputs


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
    finished_marker.write_text("done\n")

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
