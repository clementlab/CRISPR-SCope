"""Golden end-to-end regression test for the compact 30-amplicon BaF3 fixture."""

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from CRISPRSCope.output_artifacts import ARTIFACT_SPECS, INTERMEDIATE_FAMILIES


FIXTURE = Path(__file__).parent / "data" / "baf3_minimal"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@pytest.mark.integration
def test_baf3_minimal_fixture_matches_golden_baseline(tmp_path):
    required_fixture_files = [
        "R1.fastq.gz",
        "R2.fastq.gz",
        "amplicons.txt",
        "barcodes.txt",
        "settings.txt",
        "manifest.json",
        "golden.json",
        "minigenome.1.bt2",
    ]
    missing = [name for name in required_fixture_files if not (FIXTURE / name).is_file()]
    assert not missing, "Generate the fixture before running this integration test: " + ", ".join(missing)

    run_dir = tmp_path / "baf3_minimal"
    shutil.copytree(FIXTURE, run_dir)
    settings_path = run_dir / "settings.txt"
    settings_path.write_text(
        settings_path.read_text().replace(
            "write_output_manifest\tFalse", "write_output_manifest\tTrue"
        )
    )
    subprocess.run(
        [sys.executable, "-m", "CRISPRSCope.cli", "settings.txt"],
        cwd=run_dir,
        check=True,
        text=True,
        capture_output=True,
    )

    manifest = json.loads((run_dir / "manifest.json").read_text())
    golden = json.loads((run_dir / "golden.json").read_text())
    selected_categories = {}
    for row in manifest["barcodes"]:
        category = row["category"]
        selected_categories[category] = selected_categories.get(category, 0) + 1
    assert selected_categories == {"HQ_HI": 16, "HQ_LO": 16, "LQ_HI": 16, "LQ_LO": 16}

    output_root = run_dir / "run"
    output_manifest = json.loads((run_dir / "run.outputManifest.json").read_text())
    assert output_manifest["schema_version"] == 1
    assert output_manifest["status"] == "completed"
    assert [artifact["key"] for artifact in output_manifest["artifacts"]] == [
        spec.key for spec in ARTIFACT_SPECS
    ]
    artifact_statuses = {
        artifact["key"]: artifact["status"]
        for artifact in output_manifest["artifacts"]
    }
    for key in (
        "valid_amplicons",
        "editing_summary",
        "filtered_editing_summary",
        "amplicon_score",
        "editing_rate_ci",
        "editing_rate_unconditional_permutation",
        "editing_rate_unconditional_simulations",
        "report",
        "output_manifest",
    ):
        assert artifact_statuses[key] == "written"
    assert artifact_statuses["h5ad"] == "skipped"
    assert artifact_statuses["editing_rate_depth_stability"] == "skipped"
    family_summaries = {
        family["key"]: family for family in output_manifest["artifact_families"]
    }
    assert list(family_summaries) == [family.key for family in INTERMEDIATE_FAMILIES]
    assert all(family["exists"] and family["file_count"] > 0 for family in family_summaries.values())

    output_files = {
        "valid_amplicons": output_root.with_suffix(".splitReads.valid_amps.txt"),
        "editing_summary": output_root.with_suffix(".editingSummary.txt"),
        "amplicon_score": output_root.with_suffix(".amplicon_score.txt"),
        "editing_rate_ci": output_root.with_suffix(".editingRateConfidenceIntervals.txt"),
        "unconditional_permutation": output_root.with_suffix(".editingRateUnconditionalPermutation.txt"),
    }
    assert all(path.is_file() for path in output_files.values())
    assert {name: _sha256(path) for name, path in output_files.items()} == golden["hashes"]

    valid_amplicons = [line.split("\t", 1)[0] for line in output_files["valid_amplicons"].read_text().splitlines() if line]
    assert valid_amplicons == golden["invariants"]["amplicons"]
    assert len(valid_amplicons) == golden["invariants"]["amplicon_count"] == 30

    score_lines = output_files["amplicon_score"].read_text().splitlines()
    header = score_lines[0].split("\t")
    color_index = header.index("Color")
    color_counts = {}
    for line in score_lines[1:]:
        color = line.split("\t")[color_index]
        color_counts[color] = color_counts.get(color, 0) + 1
    assert len(score_lines) - 1 == golden["invariants"]["scored_cell_count"]
    assert color_counts == golden["invariants"]["score_color_counts"]
    assert color_counts.get("HQ_HI", 0) == golden["invariants"]["in_group_scored_cell_count"]
    assert sum(count for color, count in color_counts.items() if color != "HQ_HI") == golden["invariants"]["out_group_scored_cell_count"]
    assert all(color_counts.get(category, 0) > 0 for category in ["HQ_HI", "HQ_LO", "LQ_HI", "LQ_LO"])

    for suffix in [
        ".01_Log-Log.png",
        ".05_Amplicon_Score.png",
        ".11_EditingRateCoverageAdjustedEffects.png",
        ".12_EditingRateUnconditionalPermutation.png",
        ".13_EditingRateObservedCenteredPermutationSwarm.png",
    ]:
        assert Path(str(output_root) + suffix).is_file()
