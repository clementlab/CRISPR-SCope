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
	assignment = cli._classify_amplicon_assignment("ampA", "ampA", "ampA", "NA")

	assert assignment["accepted_amplicon"] == "ampA"
	assert assignment["rescued_by_partial_alignment"] is True
	assert assignment["would_rescue_under_strict"] is False


def test_amplicon_assignment_strict_mode_rejects_and_counts_would_rescue():
	assignment = cli._classify_amplicon_assignment(
		"ampA",
		"ampA",
		"ampA",
		"NA",
		require_strict_amplicon_alignment=True,
	)

	assert assignment["accepted_amplicon"] is None
	assert assignment["rescued_by_partial_alignment"] is False
	assert assignment["would_rescue_under_strict"] is True


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


def test_parse_settings_strict_alignment_debug_mode_defaults_false(tmp_path, monkeypatch):
	result = _parse_settings_for_test(tmp_path, monkeypatch)

	assert result[-2] is False


def test_parse_settings_strict_alignment_debug_mode_parses_true(tmp_path, monkeypatch):
	result = _parse_settings_for_test(
		tmp_path,
		monkeypatch,
		extra_lines=["debug_require_strict_amplicon_alignment\ttrue"],
	)

	assert result[-2] is True


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
