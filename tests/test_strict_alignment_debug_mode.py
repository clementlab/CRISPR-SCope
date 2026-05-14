import sys

import pytest

from CRISPRSCope import cli


def _write_settings(tmp_path, extra_lines=None):
	r1 = tmp_path / "r1.fastq"
	r2 = tmp_path / "r2.fastq"
	barcodes = tmp_path / "barcodes.txt"
	amplicons = tmp_path / "amplicons.txt"
	genome = tmp_path / "genome"
	for path in [r1, r2, barcodes, amplicons]:
		path.write_text("\n")
	(tmp_path / "genome.1.bt2").write_text("\n")

	output_root = tmp_path / "run"
	settings = tmp_path / "settings.txt"
	lines = [
		f"r1\t{r1}",
		f"r2\t{r2}",
		"constant1\tACGT",
		"constant2\tTGCA",
		f"barcodes\t{barcodes}",
		f"amplicons\t{amplicons}",
		f"genome\t{genome}",
		f"output_root\t{output_root}",
		"processes\t1",
	]
	if extra_lines:
		lines.extend(extra_lines)
	settings.write_text("\n".join(lines) + "\n")
	return settings


def _parse_settings_for_test(tmp_path, monkeypatch, extra_lines=None):
	settings = _write_settings(tmp_path, extra_lines)
	monkeypatch.setattr(sys, "argv", ["CRISPRSCope", str(settings)])
	monkeypatch.setattr(cli.sb, "check_output", lambda *args, **kwargs: b"ok")
	return cli.parse_settings(sys.argv)


def test_amplicon_assignment_normal_mode_rescues_partial_alignment():
	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		"ampA",
		"NA",
		r1_mean_quality=30.0,
		r2_mean_quality=30.0,
	)

	assert assignment["accepted_amplicon"] == "ampA"
	assert assignment["rescued_by_partial_alignment"] is True
	assert assignment["would_rescue_under_strict"] is False
	assert assignment["reject_reason"] is None


def test_amplicon_assignment_strict_mode_rejects_and_counts_would_rescue():
	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		"ampA",
		"NA",
		require_strict_amplicon_alignment=True,
		r1_mean_quality=30.0,
		r2_mean_quality=30.0,
	)

	assert assignment["accepted_amplicon"] is None
	assert assignment["rescued_by_partial_alignment"] is False
	assert assignment["would_rescue_under_strict"] is True
	assert assignment["reject_reason"] == "strict_alignment_required"


def test_amplicon_assignment_strict_mode_requires_both_alignment_calls():
	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		"ampA",
		"ampA",
		require_strict_amplicon_alignment=True,
	)

	assert assignment["accepted_amplicon"] == "ampA"
	assert assignment["rescued_by_partial_alignment"] is False
	assert assignment["would_rescue_under_strict"] is False
	assert assignment["reject_reason"] is None


@pytest.mark.parametrize(
	"amp1,amp2,amp1_aln,amp2_aln",
	[
		("ampA", "ampB", "ampA", "ampA"),
		("NA", "ampA", "ampA", "ampA"),
		("ampA", "ampA", "ampB", "NA"),
	],
)
def test_amplicon_assignment_rejects_invalid_or_contradictory_calls(amp1, amp2, amp1_aln, amp2_aln):
	assignment = cli._classify_amplicon_assignment(amp1, amp2, amp1_aln, amp2_aln)

	assert assignment["accepted_amplicon"] is None
	assert assignment["rescued_by_partial_alignment"] is False
	assert assignment["would_rescue_under_strict"] is False


def test_inward_boundary_allows_r2_forward_read_at_amplicon_start():
	start_lookup = {"chr1:100": "ampA"}
	end_lookup = {"chr1:199": "ampA"}

	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		"NA",
		cli.inward_alignment_amplicon("chr1", 100, "100M", False, start_lookup, end_lookup),
		r1_mean_quality=30.0,
		r2_mean_quality=30.0,
	)

	assert assignment["accepted_amplicon"] == "ampA"
	assert assignment["rescued_by_partial_alignment"] is True


def test_inward_boundary_allows_r1_reverse_read_at_amplicon_end():
	start_lookup = {"chr1:100": "ampA"}
	end_lookup = {"chr1:199": "ampA"}

	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		cli.inward_alignment_amplicon("chr1", 100, "100M", True, start_lookup, end_lookup),
		"NA",
		r1_mean_quality=30.0,
		r2_mean_quality=30.0,
	)

	assert assignment["accepted_amplicon"] == "ampA"
	assert assignment["rescued_by_partial_alignment"] is True


def test_inward_boundary_rejects_forward_read_at_amplicon_end():
	start_lookup = {"chr1:100": "ampA"}
	end_lookup = {"chr1:199": "ampA"}

	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		cli.inward_alignment_amplicon("chr1", 199, "100M", False, start_lookup, end_lookup),
		"NA",
		r1_mean_quality=30.0,
		r2_mean_quality=30.0,
	)

	assert assignment["accepted_amplicon"] is None
	assert assignment["reject_reason"] == "no_inward_boundary_support"


def test_inward_boundary_rejects_reverse_read_at_amplicon_start():
	start_lookup = {"chr1:100": "ampA"}
	end_lookup = {"chr1:199": "ampA"}

	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		cli.inward_alignment_amplicon("chr1", 1, "100M", True, start_lookup, end_lookup),
		"NA",
		r1_mean_quality=30.0,
		r2_mean_quality=30.0,
	)

	assert assignment["accepted_amplicon"] is None
	assert assignment["reject_reason"] == "no_inward_boundary_support"


def test_partial_rescue_accepts_when_both_mates_q30_or_higher():
	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		"ampA",
		"NA",
		r1_mean_quality=30.0,
		r2_mean_quality=35.0,
	)

	assert assignment["accepted_amplicon"] == "ampA"
	assert assignment["rescued_by_partial_alignment"] is True


@pytest.mark.parametrize("r1_quality,r2_quality", [(29.9, 35.0), (35.0, 29.9)])
def test_partial_rescue_rejects_when_either_mate_below_q30(r1_quality, r2_quality):
	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		"ampA",
		"NA",
		r1_mean_quality=r1_quality,
		r2_mean_quality=r2_quality,
	)

	assert assignment["accepted_amplicon"] is None
	assert assignment["rescued_by_partial_alignment"] is False
	assert assignment["reject_reason"] == "low_mean_quality"


def test_full_alignment_assignment_is_not_quality_gated():
	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		"ampA",
		"ampA",
		r1_mean_quality=10.0,
		r2_mean_quality=10.0,
	)

	assert assignment["accepted_amplicon"] == "ampA"
	assert assignment["rescued_by_partial_alignment"] is False


def test_mean_phred_quality():
	assert cli.mean_phred_quality("??") == 30.0
	assert cli.mean_phred_quality("*") is None


def test_parse_settings_strict_alignment_debug_mode_defaults_false(tmp_path, monkeypatch):
	result = _parse_settings_for_test(tmp_path, monkeypatch)

	assert result[-3] is False
	assert result[-2] == cli.PARTIAL_RESCUE_MIN_MEAN_READ_QUALITY_DEFAULT
	assert result[-4] == ""


def test_parse_settings_strict_alignment_debug_mode_parses_true(tmp_path, monkeypatch):
	result = _parse_settings_for_test(
		tmp_path,
		monkeypatch,
		extra_lines=["debug_require_strict_amplicon_alignment\ttrue"],
	)

	assert result[-3] is True


def test_parse_settings_rescue_quality_threshold_parses_float(tmp_path, monkeypatch):
	result = _parse_settings_for_test(
		tmp_path,
		monkeypatch,
		extra_lines=["partial_rescue_min_mean_read_quality\t28.5"],
	)

	assert result[-2] == 28.5


def test_parse_settings_rejected_rescue_bam_true_uses_output_root(tmp_path, monkeypatch):
	result = _parse_settings_for_test(
		tmp_path,
		monkeypatch,
		extra_lines=["debug_rejected_rescue_reads_bam\ttrue"],
	)

	assert result[-4].endswith("run.splitReads.rejected_rescue_candidates.bam")


def test_parse_settings_rejects_strict_alignment_with_assign_all(tmp_path, monkeypatch):
	with pytest.raises(ValueError, match="debug_require_strict_amplicon_alignment"):
		_parse_settings_for_test(
			tmp_path,
			monkeypatch,
			extra_lines=[
				"debug_require_strict_amplicon_alignment\ttrue",
				"assign_reads_to_all_possible_amplicons\ttrue",
			],
		)
