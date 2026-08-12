import pytest

from CRISPRSCope import cli


def test_editing_rate_ci_is_enabled_with_reproducible_defaults(tmp_path):
    settings = tmp_path / "settings.txt"
    settings.write_text("r1\treads.fastq\n")

    config = cli._parse_editing_rate_ci_config(str(settings))

    assert config.enabled is True
    assert config.bootstrap_iterations == 10_000
    assert config.confidence_level == 0.95
    assert config.seed == 42


def test_editing_rate_ci_settings_are_parsed(tmp_path):
    settings = tmp_path / "settings.txt"
    settings.write_text(
        "write_editing_rate_ci\tTrue\n"
        "editing_rate_ci_bootstrap_iterations\t2500\n"
        "editing_rate_ci_confidence_level\t0.9\n"
        "editing_rate_ci_seed\t123\n"
    )

    config = cli._parse_editing_rate_ci_config(str(settings))

    assert config.enabled is True
    assert config.bootstrap_iterations == 2500
    assert config.confidence_level == 0.9
    assert config.seed == 123


def test_editing_rate_ci_can_be_disabled(tmp_path):
    settings = tmp_path / "settings.txt"
    settings.write_text("write_editing_rate_ci\tFalse\n")

    config = cli._parse_editing_rate_ci_config(str(settings))

    assert config.enabled is False


@pytest.mark.parametrize(
    "line,match",
    [
        ("editing_rate_ci_bootstrap_iterations\t99\n", "must be >= 100"),
        ("editing_rate_ci_confidence_level\t1.0\n", "greater than 0 and less than 1"),
        ("editing_rate_ci_seed\t-1\n", "must be >= 0"),
    ],
)
def test_editing_rate_ci_settings_reject_invalid_values(tmp_path, line, match):
    settings = tmp_path / "settings.txt"
    settings.write_text(line)

    with pytest.raises(ValueError, match=match):
        cli._parse_editing_rate_ci_config(str(settings))


def test_editing_rate_depth_stability_is_enabled_by_default(tmp_path):
    settings = tmp_path / "settings.txt"
    settings.write_text("write_editing_rate_ci\tFalse\n")

    config = cli._parse_editing_rate_depth_stability_config(str(settings))

    assert config.enabled is True
    assert config.iterations == 1_000
    assert config.percentages == (10.0, 25.0, 50.0, 75.0, 90.0)
    assert config.confidence_level == 0.95
    assert config.seed == 42
    assert config.relative_min_hq_edit_pct == 1.0


def test_editing_rate_depth_stability_can_be_disabled(tmp_path):
    settings = tmp_path / "settings.txt"
    settings.write_text("write_editing_rate_depth_stability\tFalse\n")

    config = cli._parse_editing_rate_depth_stability_config(str(settings))

    assert config.enabled is False


def test_editing_rate_depth_stability_settings_are_parsed(tmp_path):
    settings = tmp_path / "settings.txt"
    settings.write_text(
        "write_editing_rate_depth_stability\tTrue\n"
        "editing_rate_depth_stability_iterations\t250\n"
        "editing_rate_depth_stability_percentages\t5,20.5,80\n"
        "editing_rate_depth_stability_relative_min_hq_edit_pct\t2.5\n"
        "editing_rate_ci_confidence_level\t0.9\n"
        "editing_rate_ci_seed\t123\n"
    )

    config = cli._parse_editing_rate_depth_stability_config(str(settings))

    assert config.enabled is True
    assert config.iterations == 250
    assert config.percentages == (5.0, 20.5, 80.0)
    assert config.relative_min_hq_edit_pct == 2.5
    assert config.confidence_level == 0.9
    assert config.seed == 123


@pytest.mark.parametrize(
    "line,match",
    [
        ("editing_rate_depth_stability_iterations\t99\n", "must be >= 100"),
        ("editing_rate_depth_stability_percentages\t10,,25\n", "Invalid value"),
        ("editing_rate_depth_stability_percentages\t0,50\n", "greater than 0"),
        ("editing_rate_depth_stability_percentages\t10,100\n", "less than 100"),
        ("editing_rate_depth_stability_percentages\t25,10\n", "strictly increasing"),
        ("editing_rate_depth_stability_percentages\t10,10\n", "strictly increasing"),
        ("editing_rate_depth_stability_relative_min_hq_edit_pct\t0\n", "greater than 0"),
        ("editing_rate_depth_stability_relative_min_hq_edit_pct\t100.1\n", "no greater than 100"),
        ("editing_rate_depth_stability_relative_min_hq_edit_pct\tnan\n", "greater than 0"),
    ],
)
def test_editing_rate_depth_stability_settings_reject_invalid_values(tmp_path, line, match):
    settings = tmp_path / "settings.txt"
    settings.write_text(line)

    with pytest.raises(ValueError, match=match):
        cli._parse_editing_rate_depth_stability_config(str(settings))
