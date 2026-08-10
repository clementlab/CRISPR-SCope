import pytest

from CRISPRSCope import cli


def test_editing_rate_ci_is_opt_in_with_reproducible_defaults(tmp_path):
    settings = tmp_path / "settings.txt"
    settings.write_text("r1\treads.fastq\n")

    config = cli._parse_editing_rate_ci_config(str(settings))

    assert config.enabled is False
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
