import gzip
import importlib.util
from pathlib import Path

import pytest


GENERATOR_PATH = Path(__file__).parent / "data" / "baf3_minimal" / "generate_fixture.py"
SPEC = importlib.util.spec_from_file_location("baf3_fixture_generator", GENERATOR_PATH)
fixture_generator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fixture_generator)


def _write_fastq(path, records):
    with gzip.open(path, "wt") as handle:
        for header, sequence in records:
            handle.write(f"@{header}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def test_select_barcodes_uses_source_rank_cutoff_and_stable_order(monkeypatch):
    monkeypatch.setattr(fixture_generator, "CELLS_PER_CATEGORY", 1)
    monkeypatch.setattr(fixture_generator, "SOURCE_AMP_SCORE_MAX_BARCODE_RANK", 5)
    rows = [
        {"barcode": "b", "Amplicon Score": "1.0", "Barcode Rank": "2"},
        {"barcode": "a", "Amplicon Score": "1.0", "Barcode Rank": "2"},
        {"barcode": "c", "Amplicon Score": "1.0", "Barcode Rank": "6"},
        {"barcode": "d", "Amplicon Score": "0.5", "Barcode Rank": "1"},
        {"barcode": "e", "Amplicon Score": "0.5", "Barcode Rank": "7"},
    ]

    assert fixture_generator.select_barcodes(rows) == {
        "a": "HQ_HI",
        "c": "HQ_LO",
        "d": "LQ_HI",
        "e": "LQ_LO",
    }


def test_select_read_ids_caps_each_cell_and_preserves_amplicon_assignment(tmp_path, monkeypatch):
    monkeypatch.setattr(fixture_generator, "PAIRS_PER_CELL", 2)
    seq_dir = tmp_path / "settings.txt.seq_by_amplicon"
    seq_dir.mkdir()
    barcode_a = "A" * 18
    barcode_c = "C" * 18
    _write_fastq(
        seq_dir / "03_reads_all_cells.ampA.r1.fq.gz",
        [(f"machine:1:{barcode_a}", "ACGT"), (f"machine:2:{barcode_a}", "ACGT"), (f"machine:3:{barcode_c}", "ACGT")],
    )
    _write_fastq(
        seq_dir / "03_reads_all_cells.ampB.r1.fq.gz",
        [(f"machine:4:{barcode_a}", "ACGT"), (f"machine:5:{barcode_c}", "ACGT")],
    )

    selected, counts = fixture_generator.select_read_ids(
        tmp_path,
        ["ampA", "ampB"],
        {barcode_a: "HQ_HI", barcode_c: "LQ_LO"},
    )

    assert len(selected) == 4
    assert sum(counts["ampA"].values()) + sum(counts["ampB"].values()) == 4


def test_extract_raw_pairs_preserves_selected_pairs_and_rejects_mate_mismatch(tmp_path):
    source_r1 = tmp_path / "source-r1.fq.gz"
    source_r2 = tmp_path / "source-r2.fq.gz"
    _write_fastq(source_r1, [("machine:1 1:N:0:1", "AAAA"), ("machine:2 1:N:0:1", "CCCC")])
    _write_fastq(source_r2, [("machine:1 2:N:0:1", "TTTT"), ("machine:2 2:N:0:1", "GGGG")])

    assert fixture_generator.extract_raw_pairs(source_r1, source_r2, {"machine:2"}, tmp_path / "out") == 1
    assert [record[0] for record in fixture_generator.read_fastq_records(tmp_path / "out" / "R1.fastq.gz")] == [
        "@machine:2 1:N:0:1\n"
    ]

    _write_fastq(source_r2, [("machine:wrong 2:N:0:1", "TTTT")])
    with pytest.raises(ValueError, match="different record counts|mate mismatch"):
        fixture_generator.extract_raw_pairs(source_r1, source_r2, {"machine:1"}, tmp_path / "bad")


def test_intermediate_header_requires_appended_barcode():
    assert fixture_generator.intermediate_read_id_and_barcode("@machine:1:" + "A" * 18) == (
        "machine:1",
        "A" * 18,
    )
    with pytest.raises(ValueError, match="Unexpected barcode suffix"):
        fixture_generator.intermediate_read_id_and_barcode("@machine:1:not-a-barcode")


def test_barcode_halves_emits_unique_nine_base_whitelist_entries():
    assert fixture_generator.barcode_halves(["AAAAAAAAACCCCCCCCC", "AAAAAAAAAGGGGGGGGG"]) == [
        "AAAAAAAAA",
        "CCCCCCCCC",
        "GGGGGGGGG",
    ]


def test_tracked_settings_match_fixture_generator(tmp_path):
    generated = tmp_path / "settings.txt"
    fixture_generator.write_settings(generated)
    assert generated.read_text() == (GENERATOR_PATH.parent / "settings.txt").read_text()
