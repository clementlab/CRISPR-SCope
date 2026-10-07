import numpy as np
import pandas as pd
import pytest

from CRISPRSCope import crispresso
from CRISPRSCope.fastq_processing import parse_fq_file_pair
from CRISPRSCope.h5ad.builder import CRISPRSCopeAnnDataBuilder
from CRISPRSCope.h5ad.loaders import load_amplicons


def test_load_amplicons_preserves_required_and_optional_columns(tmp_path):
    two_column = tmp_path / "two-column.tsv"
    two_column.write_text("ampA\tACGT\n")
    two_column_frame = load_amplicons(two_column)
    assert two_column_frame.loc["ampA", "sequence"] == "ACGT"
    assert pd.isna(two_column_frame.loc["ampA", "guide"])

    four_column = tmp_path / "four-column.tsv"
    four_column.write_text("ampA\tACGT\tAC\t2\n")
    four_column_frame = load_amplicons(four_column)
    assert four_column_frame.loc["ampA", "sequence"] == "ACGT"
    assert four_column_frame.loc["ampA", "guide"] == "AC"
    assert four_column_frame.loc["ampA", "reference_allele_count"] == "2"


def test_h5ad_zygosity_marks_unobserved_cell_amplicon_as_no_data():
    builder = CRISPRSCopeAnnDataBuilder(
        config={
            "analysis_parameters": {
                "zygosity": {
                    "wt_max_mod_pct": 20.0,
                    "het_max_mod_pct": 80.0,
                    "hom_min_mod_pct": 80.0,
                    "compound_het_min_allele2_pct": 20.0,
                }
            }
        },
        settings={},
        amplicons=pd.DataFrame(
            {"sequence": ["ACGT"]}, index=pd.Index(["ampA"], name="amplicon_name")
        ),
        editing_summary=pd.DataFrame(
            {
                "totCount.ampA": [10, 0],
                "modPct.ampA": [100.0, np.nan],
            },
            index=["called", "uncovered"],
        ),
        quality_scores=pd.DataFrame(
            {"Color": ["HQ_HI", "HQ_HI"]}, index=["called", "uncovered"]
        ),
        allele_parquet_paths=[],
    )

    adata = builder.build()
    assert adata.layers["zygosity"].tolist() == [[2], [-1]]


@pytest.mark.parametrize(
    "r1_records,r2_records,error_message",
    [
        (["read1", "read2"], ["read1"], "different numbers of records"),
        (["read1"], ["other_read"], "different read identifiers"),
    ],
)
def test_raw_fastq_parser_rejects_unsynchronized_mates(
    tmp_path, r1_records, r2_records, error_message
):
    def write_fastq(path, names):
        with open(path, "w") as handle:
            for name in names:
                handle.write(f"@{name}\nACGT\n+\nIIII\n")

    r1_path = tmp_path / "r1.fastq"
    r2_path = tmp_path / "r2.fastq"
    write_fastq(r1_path, r1_records)
    write_fastq(r2_path, r2_records)

    with pytest.raises(ValueError, match=error_message):
        parse_fq_file_pair(
            (
                str(r1_path), str(r2_path), 0, {}, "CCCC", "GGGG", 4, 4,
                "TTTT", "AAAA", str(tmp_path / "run"),
            )
        )


def test_first_pass_parser_discards_reads_without_expected_amplicon_arms(tmp_path, monkeypatch):
    run_folder = tmp_path / "CRISPResso_on_ampA"
    run_folder.mkdir()
    (run_folder / "CRISPResso_output.fastq.gz").write_text(
        "@read:cellA\n"
        + "G" * 70
        + "\n+ ALN=Reference DEL= INS= SUB= ALN_REF=Reference \n"
        + "I" * 70
        + "\n"
    )
    info_file = tmp_path / "amplicon_info.tsv"
    info_file.write_text("name\tamp_seqs\nampA\t" + "A" * 30 + "C" * 10 + "T" * 30 + "\n")
    (tmp_path / "run.seq_by_amplicon").mkdir()
    monkeypatch.setattr(crispresso, "get_wildtype_allele", lambda _: None)

    crispresso.parse_one_crispresso_output(
        {
            "amplicon_name": "ampA",
            "amplicon_info_file": str(info_file),
            "crispresso_run_folder": str(run_folder),
            "input_ref_allele_counts": "1",
            "min_num_reads_per_cell": 0,
            "min_allele_support": 2,
            "min_reads_per_amplicon_for_genotype": 8,
            "ignore_substitutions": False,
            "output_root": str(tmp_path / "run"),
            "min_reads_per_amplicon_per_cell": 0,
        }
    )

    summary_lines = (tmp_path / "CRISPResso_on_ampA.summ").read_text().splitlines()
    assert len(summary_lines) == 1
