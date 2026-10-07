#!/usr/bin/env python3
"""Build and validate the deterministic minimal BaF3 integration fixture.

The source run is intentionally outside this fixture: regeneration requires the
original BaF3 raw FASTQs and completed run intermediates.  Running the fixture
does not; it uses only the committed files in this directory.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from collections import Counter
from itertools import zip_longest
from pathlib import Path
from typing import Iterable, Iterator


SEED = 20260914
CELLS_PER_CATEGORY = 16
PAIRS_PER_CELL = 450
# Select the balanced cell panel from the completed production run's categories.
# Its rank scale is not comparable to the 64-cell fixture's recomputed ranks.
SOURCE_AMP_SCORE_MAX_BARCODE_RANK = 6000
# This deliberately sits below the 64-cell fixture size so the end-to-end
# pipeline exercises both high- and low-depth score categories.
FIXTURE_AMP_SCORE_MAX_BARCODE_RANK = 32
CATEGORIES = ("HQ_HI", "HQ_LO", "LQ_HI", "LQ_LO")
FIXTURE_DIR = Path(__file__).resolve().parent


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_key(*parts: str) -> str:
    return hashlib.sha256("\0".join((str(SEED), *parts)).encode()).hexdigest()


def read_fastq_records(path: Path) -> Iterator[tuple[str, str, str, str]]:
    """Yield validated four-line FASTQ records from a gzip-compressed file."""
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not sequence or not plus or not quality:
                raise ValueError(f"Truncated FASTQ record in {path}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record in {path}: {header!r}")
            yield header, sequence, plus, quality


def raw_read_id(header: str) -> str:
    """Return the Illumina read identifier shared by R1 and R2."""
    return header.strip().split()[0].lstrip("@")


def intermediate_read_id_and_barcode(header: str) -> tuple[str, str]:
    """Recover original read ID and appended 18-base cell barcode."""
    token = raw_read_id(header)
    try:
        read_id, barcode = token.rsplit(":", 1)
    except ValueError as error:
        raise ValueError(f"Unexpected intermediate header: {header!r}") from error
    if len(barcode) != 18 or any(base not in "ACGT" for base in barcode):
        raise ValueError(f"Unexpected barcode suffix in intermediate header: {header!r}")
    return read_id, barcode


def parse_score_table(path: Path) -> list[dict[str, str]]:
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    if header[0] != "":
        raise ValueError(f"Expected barcode index column in {path}")
    result = []
    for line in lines[1:]:
        values = line.split("\t")
        result.append(dict(zip(header[1:], values[1:])) | {"barcode": values[0]})
    return result


def category_at_rank_cutoff(row: dict[str, str], rank_cutoff: int) -> str:
    """Classify a score-table row at the supplied depth-rank cutoff."""
    high_score = float(row["Amplicon Score"]) >= 2.0 / 3.0
    high_depth = int(row["Barcode Rank"]) <= rank_cutoff
    if high_score:
        return "HQ_HI" if high_depth else "HQ_LO"
    return "LQ_HI" if high_depth else "LQ_LO"


def select_barcodes(score_rows: Iterable[dict[str, str]]) -> dict[str, str]:
    """Select a balanced panel from the completed run's four categories."""
    selected: dict[str, str] = {}
    for category in CATEGORIES:
        candidates = sorted(
            (
                row
                for row in score_rows
                if category_at_rank_cutoff(row, SOURCE_AMP_SCORE_MAX_BARCODE_RANK)
                == category
            ),
            key=lambda row: (int(row["Barcode Rank"]), row["barcode"]),
        )
        if len(candidates) < CELLS_PER_CATEGORY:
            raise ValueError(
                f"Need {CELLS_PER_CATEGORY} {category} barcodes; found {len(candidates)}"
            )
        selected.update({row["barcode"]: category for row in candidates[:CELLS_PER_CATEGORY]})
    return selected


def select_read_ids(
    source_run: Path,
    amplicons: list[str],
    barcode_categories: dict[str, str],
) -> tuple[set[str], dict[str, dict[str, int]]]:
    """Cap each selected cell while preserving its original amplicon breadth."""
    seq_dir = source_run / "settings.txt.seq_by_amplicon"
    candidates = {barcode: [] for barcode in barcode_categories}
    for amplicon in amplicons:
        source = seq_dir / f"03_reads_all_cells.{amplicon}.r1.fq.gz"
        if not source.is_file():
            raise FileNotFoundError(source)
        for header, *_ in read_fastq_records(source):
            read_id, barcode = intermediate_read_id_and_barcode(header)
            if barcode in candidates:
                candidates[barcode].append((stable_key(barcode, read_id), read_id, amplicon))

    selected: set[str] = set()
    counts = {amplicon: {category: 0 for category in CATEGORIES} for amplicon in amplicons}
    for barcode, values in candidates.items():
        for _, read_id, amplicon in sorted(values)[:PAIRS_PER_CELL]:
            if read_id in selected:
                raise ValueError(f"Selected read {read_id} for more than one amplicon")
            selected.add(read_id)
            counts[amplicon][barcode_categories[barcode]] += 1
    if not selected:
        raise ValueError("No read pairs were selected")
    missing_amplicons = [amplicon for amplicon, values in counts.items() if not sum(values.values())]
    if missing_amplicons:
        raise ValueError("No selected reads for amplicons: " + ", ".join(missing_amplicons))
    return selected, counts


def _deterministic_gzip_writer(path: Path):
    raw = path.open("wb")
    gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return raw, gz


def extract_raw_pairs(source_r1: Path, source_r2: Path, selected_ids: set[str], destination: Path) -> int:
    """Write selected original pairs and ensure input mates remain synchronized."""
    destination.mkdir(parents=True, exist_ok=True)
    raw_r1, out_r1 = _deterministic_gzip_writer(destination / "R1.fastq.gz")
    raw_r2, out_r2 = _deterministic_gzip_writer(destination / "R2.fastq.gz")
    found: set[str] = set()
    try:
        for r1, r2 in zip_longest(read_fastq_records(source_r1), read_fastq_records(source_r2)):
            if r1 is None or r2 is None:
                raise ValueError("Input R1 and R2 FASTQs contain different record counts")
            r1_id = raw_read_id(r1[0])
            r2_id = raw_read_id(r2[0])
            if r1_id != r2_id:
                raise ValueError(f"Input mate mismatch: {r1_id} != {r2_id}")
            if r1_id in selected_ids:
                out_r1.write("".join(r1).encode())
                out_r2.write("".join(r2).encode())
                found.add(r1_id)
    finally:
        out_r1.close()
        out_r2.close()
        raw_r1.close()
        raw_r2.close()
    missing = selected_ids - found
    if missing:
        raise ValueError(f"Failed to recover {len(missing)} selected raw read pairs")
    return len(found)


def amplicon_names(path: Path) -> list[str]:
    return [line.split("\t", 1)[0] for line in path.read_text().splitlines() if line]


def barcode_halves(barcodes: Iterable[str]) -> list[str]:
    """Return the 9-base whitelist entries required for full 18-base barcodes."""
    halves = set()
    for barcode in barcodes:
        if len(barcode) != 18 or any(base not in "ACGT" for base in barcode):
            raise ValueError(f"Expected an 18-base barcode, got {barcode!r}")
        halves.update((barcode[:9], barcode[9:]))
    return sorted(halves)


def write_minigenome(amplicons_file: Path, destination: Path) -> None:
    with amplicons_file.open() as source, destination.open("w") as output:
        for line in source:
            name, sequence, *_ = line.rstrip("\n").split("\t")
            output.write(f">{name}\n{sequence.split(',')[0]}\n")


def write_settings(destination: Path) -> None:
    destination.write_text(
        "\n".join(
            [
                "r1\tR1.fastq.gz",
                "r2\tR2.fastq.gz",
                "barcodes\tbarcodes.txt",
                "amplicons\tamplicons.txt",
                "bowtie2_index\tminigenome",
                "constant1\tGTACGTACGAGTC",
                "constant2\tGTACTCGCAGTAGTC",
                "output_root\trun",
                "processes\t1",
                "cache_mode\tauto",
                "allowBarcodeMismatches\tFalse",
                "primer_lookup_len\t18",
                "adapter_DNA\tTGTCTCTTATACACATCTCCGAGCCCACGAG",
                "keep_intermediate_files\tFalse",
                "ignore_substitutions\tFalse",
                "assign_reads_to_all_possible_amplicons\tFalse",
                "suppress_sub_crispresso_plots\tFalse",
                "min_total_reads_per_barcode\t10",
                "min_reads_per_amplicon_per_cell\t0",
                "min_reads_per_amplicon_for_genotype\t8",
                "min_allele_support\t2",
                "amplicon_score_min_reads_per_amplicon\t5",
                "amplicon_score_min_covered_fraction\t0.6666666666666666",
                f"amplicon_score_max_barcode_rank\t{FIXTURE_AMP_SCORE_MAX_BARCODE_RANK}",
                "include_high_score_high_depth\tTrue",
                "include_high_score_low_depth\tFalse",
                "include_low_score_high_depth\tFalse",
                "include_low_score_low_depth\tFalse",
                "write_editing_rate_ci\tTrue",
                "editing_rate_ci_bootstrap_iterations\t100",
                "editing_rate_ci_permutation_iterations\t100",
                "editing_rate_ci_confidence_level\t0.95",
                "editing_rate_ci_seed\t42",
                "editing_rate_ci_coverage_exact_max_reads\t10",
                "editing_rate_ci_coverage_bin_width_reads\t5",
                "write_editing_rate_depth_stability\tFalse",
                "write_h5ad\tFalse",
                "write_output_manifest\tFalse",
                "",
            ]
        )
    )


def source_defaults() -> tuple[Path, Path]:
    workspace = FIXTURE_DIR.parents[4]
    return (
        workspace / "analysis/01_run_on_20200804_BaF3_revision",
        workspace / "data/20200804_BaF3_revision/data",
    )


def build_fixture(destination: Path, source_run: Path, source_data: Path, overwrite: bool) -> dict:
    tracked_outputs = ["R1.fastq.gz", "R2.fastq.gz", "barcodes.txt", "amplicons.txt", "minigenome.fa", "manifest.json"]
    if not overwrite and any((destination / name).exists() for name in tracked_outputs):
        raise FileExistsError(f"Fixture exists at {destination}; use --overwrite to replace it")
    destination.mkdir(parents=True, exist_ok=True)
    source_amplicons = source_data.parent / "amplicons.txt"
    source_scores = source_run / "settings.txt.amplicon_score.txt"
    source_r1 = source_data / "BaF3-NSG_S1_L001_R1_001.fastq.gz"
    source_r2 = source_data / "BaF3-NSG_S1_L001_R2_001.fastq.gz"
    for path in (source_amplicons, source_scores, source_r1, source_r2):
        if not path.is_file():
            raise FileNotFoundError(path)

    amplicons = amplicon_names(source_amplicons)
    if len(amplicons) != 30 or len(set(amplicons)) != 30:
        raise ValueError("Expected exactly 30 unique BaF3 amplicons")
    barcode_categories = select_barcodes(parse_score_table(source_scores))
    selected_ids, per_amplicon_counts = select_read_ids(source_run, amplicons, barcode_categories)
    recovered = extract_raw_pairs(source_r1, source_r2, selected_ids, destination)
    shutil.copyfile(source_amplicons, destination / "amplicons.txt")
    (destination / "barcodes.txt").write_text(
        "\n".join(barcode_halves(barcode_categories)) + "\n"
    )
    write_minigenome(destination / "amplicons.txt", destination / "minigenome.fa")
    subprocess.run(["bowtie2-build", str(destination / "minigenome.fa"), str(destination / "minigenome")], check=True)
    write_settings(destination / "settings.txt")
    manifest = {
        "fixture": "baf3_minimal",
        "selection_seed": SEED,
        "cells_per_category": CELLS_PER_CATEGORY,
        "pairs_per_cell": PAIRS_PER_CELL,
        "amplicon_score_max_barcode_rank": FIXTURE_AMP_SCORE_MAX_BARCODE_RANK,
        "source_sha256": {"r1": sha256_file(source_r1), "r2": sha256_file(source_r2)},
        "barcodes": [{"barcode": barcode, "category": category} for barcode, category in barcode_categories.items()],
        "amplicons": amplicons,
        "per_amplicon_pair_counts": per_amplicon_counts,
        "selected_pair_count": recovered,
    }
    (destination / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def record_golden(destination: Path, output_root: Path, overwrite: bool) -> dict:
    output_files = {
        "valid_amplicons": output_root.with_suffix(".splitReads.valid_amps.txt"),
        "allele_call_qc": output_root.with_suffix(".alleleCallQC.txt"),
        "editing_summary": output_root.with_suffix(".editingSummary.txt"),
        "editing_summary_pseudobulk": output_root.with_suffix(".editingSummaryPseudobulk.txt"),
        "filtered_editing_summary": output_root.with_suffix(".filteredEditingSummary.txt"),
        "filtered_editing_summary_pseudobulk": output_root.with_suffix(".filteredEditingSummaryPseudobulk.txt"),
        "amplicon_score": output_root.with_suffix(".amplicon_score.txt"),
        "editing_rate_ci": output_root.with_suffix(".editingRateConfidenceIntervals.txt"),
        "unconditional_permutation": output_root.with_suffix(".editingRateUnconditionalPermutation.txt"),
    }
    missing = [str(path) for path in output_files.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing golden outputs: " + ", ".join(missing))
    golden_path = destination / "golden.json"
    if golden_path.exists() and not overwrite:
        raise FileExistsError(f"Golden expectations exist at {golden_path}; use --overwrite")
    manifest = json.loads((destination / "manifest.json").read_text())
    valid_amplicons = [line.split("\t", 1)[0] for line in output_files["valid_amplicons"].read_text().splitlines() if line]
    score_lines = output_files["amplicon_score"].read_text().splitlines()
    score_header = score_lines[0].split("\t")
    color_index = score_header.index("Color")
    group_counts = Counter(line.split("\t")[color_index] for line in score_lines[1:] if line)
    golden = {
        "hashes": {name: sha256_file(path) for name, path in output_files.items()},
        "invariants": {
            "amplicon_count": len(valid_amplicons),
            "amplicons": valid_amplicons,
            "selected_pair_count": manifest["selected_pair_count"],
            "scored_cell_count": len(score_lines) - 1,
            "score_color_counts": dict(sorted(group_counts.items())),
            "in_group_scored_cell_count": group_counts["HQ_HI"],
            "out_group_scored_cell_count": sum(
                count for color, count in group_counts.items() if color != "HQ_HI"
            ),
        },
    }
    golden_path.write_text(json.dumps(golden, indent=2, sort_keys=True) + "\n")
    return golden


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=FIXTURE_DIR)
    parser.add_argument("--source-run", type=Path)
    parser.add_argument("--source-data", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--record-golden", type=Path, metavar="OUTPUT_ROOT")
    args = parser.parse_args(argv)
    if args.record_golden:
        golden = record_golden(args.destination, args.record_golden, args.overwrite)
        print(json.dumps(golden["invariants"], indent=2, sort_keys=True))
        return 0
    source_run, source_data = source_defaults()
    manifest = build_fixture(
        args.destination,
        args.source_run or source_run,
        args.source_data or source_data,
        args.overwrite,
    )
    print(json.dumps({"selected_pair_count": manifest["selected_pair_count"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
