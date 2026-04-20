from pathlib import Path

import pytest

from CRISPRSCope.h5ad.paths import default_h5ad_output_path, resolve_h5ad_input_paths


def test_default_h5ad_output_path_uses_output_root(tmp_path):
    output_root = tmp_path / "run_prefix"
    assert default_h5ad_output_path(str(output_root)) == Path(f"{output_root}.h5ad")


def test_resolve_h5ad_input_paths_uses_output_root_when_it_is_settings_file(tmp_path):
    settings_path = tmp_path / "settings.txt"
    settings_path.write_text("amplicons\tamplicons.tsv\n")

    resolved = resolve_h5ad_input_paths(str(settings_path))

    assert resolved.settings_path == settings_path
    assert resolved.editing_summary_path == tmp_path / "settings.txt.filteredEditingSummary.txt"
    assert resolved.quality_scores_path == tmp_path / "settings.txt.amplicon_score.txt"
    assert resolved.crispresso_filtered_dir == tmp_path / "settings.txt.crispresso.filtered"


def test_resolve_h5ad_input_paths_requires_settings_for_custom_output_root(tmp_path):
    output_root = tmp_path / "custom_prefix"

    with pytest.raises(FileNotFoundError):
        resolve_h5ad_input_paths(str(output_root))


def test_resolve_h5ad_input_paths_accepts_explicit_settings_path_for_custom_output_root(tmp_path):
    output_root = tmp_path / "custom_prefix"
    settings_path = tmp_path / "settings.txt"
    settings_path.write_text("amplicons\tamplicons.tsv\n")

    resolved = resolve_h5ad_input_paths(str(output_root), settings_path=str(settings_path))

    assert resolved.settings_path == settings_path
    assert resolved.editing_summary_path == tmp_path / "custom_prefix.filteredEditingSummary.txt"
    assert resolved.quality_scores_path == tmp_path / "custom_prefix.amplicon_score.txt"
    assert resolved.crispresso_filtered_dir == tmp_path / "custom_prefix.crispresso.filtered"
