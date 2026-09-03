import numpy as np
import pandas as pd
import pytest

from CRISPRSCope import cli


def _totals(rows, n_amplicons):
    return pd.DataFrame(
        rows,
        columns=[f"totCount.amp{i}" for i in range(n_amplicons)],
        index=[f"cell{i}" for i in range(len(rows))],
    )


def test_supported_breadth_rejects_single_amplicon_jackpot():
    totals = _totals(
        [
            [1000, 0, 0, 0, 0, 0],
            [5, 5, 5, 5, 0, 0],
        ],
        6,
    )

    result = cli.generate_amplicon_score(
        totals,
        min_reads_per_amplicon_per_cell=0,
        min_total_reads_per_barcode=0,
    )

    assert result.loc["cell0", "Read Count"] == 1000
    assert result.loc["cell0", "Amplicon Score"] == pytest.approx(1 / 6)
    assert result.loc["cell0", "Color"] == "LQ_HI"
    assert result.loc["cell1", "Amplicon Score"] == pytest.approx(4 / 6)
    assert result.loc["cell1", "Color"] == "HQ_HI"


def test_default_two_thirds_boundary_uses_required_supported_count():
    totals = _totals(
        [
            [5] * 20 + [0] * 10,
            [5] * 19 + [0] * 11,
        ],
        30,
    )

    result = cli.generate_amplicon_score(totals, 0, 0)

    assert result.loc["cell0", "Supported Amplicons"] == 20
    assert result.loc["cell0", "Color"] == "HQ_HI"
    assert result.loc["cell1", "Supported Amplicons"] == 19
    assert result.loc["cell1", "Color"] == "LQ_HI"


def test_custom_ten_reads_at_three_fourths_requires_23_of_30():
    totals = _totals(
        [
            [10] * 23 + [0] * 7,
            [10] * 22 + [0] * 8,
        ],
        30,
    )
    config = cli.AmpliconScoreConfig(
        min_reads_per_amplicon=10,
        min_covered_fraction=0.75,
        max_barcode_rank=100,
    )

    result = cli.generate_amplicon_score(totals, 0, 0, config=config)

    assert result.loc["cell0", "Color"] == "HQ_HI"
    assert result.loc["cell1", "Color"] == "LQ_HI"


def test_rank_ties_are_resolved_by_barcode():
    totals = pd.DataFrame(
        [[5, 5], [5, 5]],
        columns=["totCount.a", "totCount.b"],
        index=["z_cell", "a_cell"],
    )

    result = cli.generate_amplicon_score(totals, 0, 0)

    assert list(result.index) == ["a_cell", "z_cell"]
    assert result.loc["a_cell", "Barcode Rank"] == 1


def test_existing_all_amplicon_gate_rejects_missing_and_low_counts():
    totals = pd.DataFrame(
        {
            "totCount.a": [5, 5, 5],
            "totCount.b": [5, np.nan, 0],
            "totCount.failed": [np.nan, np.nan, np.nan],
        },
        index=["pass", "missing", "low"],
    )

    result = cli.generate_amplicon_score(
        totals,
        min_reads_per_amplicon_per_cell=1,
        min_total_reads_per_barcode=0,
    )

    assert list(result.index) == ["pass"]
    assert result.loc["pass", "Usable Amplicons"] == 2


def test_no_usable_amplicons_returns_empty_schema():
    totals = pd.DataFrame(
        {"totCount.failed": [np.nan]},
        index=["cell"],
    )

    result = cli.generate_amplicon_score(totals, 0, 0)

    assert result.empty
    assert list(result.columns) == [
        "Amplicon Score",
        "Supported Amplicons",
        "Usable Amplicons",
        "Read Count",
        "Barcode Rank",
        "Color",
    ]


def test_empty_configured_group_fails_before_second_pass():
    classified = pd.DataFrame({"Color": ["LQ_HI", "LQ_LO"]})

    with pytest.raises(ValueError, match="No barcodes match"):
        cli._require_selected_barcodes(classified, ["HQ_HI"])


def test_selected_barcode_fingerprint_is_order_independent():
    assert cli._barcode_set_sha256(["cellB", "cellA"]) == cli._barcode_set_sha256(
        ["cellA", "cellB"]
    )
    assert cli._barcode_set_sha256(["cellA"]) != cli._barcode_set_sha256(
        ["cellA", "cellB"]
    )


def test_amplicon_score_plot_shows_configured_boundaries(tmp_path):
    output_root = str(tmp_path / "run")
    pd.DataFrame(
        {
            "Amplicon Score": [0.8, 0.5],
            "Supported Amplicons": [24, 15],
            "Usable Amplicons": [30, 30],
            "Read Count": [200, 100],
            "Barcode Rank": [1, 2],
            "Color": ["HQ_HI", "LQ_HI"],
        },
        index=["cellA", "cellB"],
    ).to_csv(f"{output_root}.amplicon_score.txt", sep="\t")
    config = cli.AmpliconScoreConfig(
        min_reads_per_amplicon=10,
        min_covered_fraction=0.75,
        max_barcode_rank=42,
    )

    metadata = cli.plot_amp_score(output_root, config=config)

    assert (tmp_path / "run.05_Amplicon_Score.png").exists()
    assert (tmp_path / "run.05_Amplicon_Score.pdf").exists()
    assert metadata.title == "Supported Amplicon Breadth Plot"
    axes = cli.plt.gca()
    assert any(np.allclose(line.get_ydata(), [0.75, 0.75]) for line in axes.lines)
    assert any(np.allclose(line.get_xdata(), [42, 42]) for line in axes.lines)


def test_parse_amplicon_score_settings_defaults_and_custom_values(tmp_path):
    defaults_path = tmp_path / "defaults.txt"
    defaults_path.write_text("# defaults\n")
    defaults = cli._parse_amplicon_score_config(str(defaults_path))
    assert defaults == cli.AmpliconScoreConfig()

    custom_path = tmp_path / "custom.txt"
    custom_path.write_text(
        "amplicon_score_min_reads_per_amplicon\t10\n"
        "amplicon_score_min_covered_fraction\t0.75\n"
        "amplicon_score_max_barcode_rank\t8000\n"
    )
    custom = cli._parse_amplicon_score_config(str(custom_path))
    assert custom.min_reads_per_amplicon == 10
    assert custom.min_covered_fraction == 0.75
    assert custom.max_barcode_rank == 8000


@pytest.mark.parametrize(
    "line,match",
    [
        ("amplicon_score_min_reads_per_amplicon\t0\n", "must be >= 1"),
        ("amplicon_score_min_covered_fraction\t0\n", "must be > 0 and <= 1"),
        ("amplicon_score_min_covered_fraction\t1.1\n", "must be > 0 and <= 1"),
        ("amplicon_score_min_covered_fraction\tnan\n", "must be > 0 and <= 1"),
        ("amplicon_score_max_barcode_rank\t0\n", "must be >= 1"),
    ],
)
def test_invalid_amplicon_score_settings_fail(tmp_path, line, match):
    settings = tmp_path / "settings.txt"
    settings.write_text(line)

    with pytest.raises(ValueError, match=match):
        cli._parse_amplicon_score_config(str(settings))
