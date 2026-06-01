import os
import sys

import pytest

from CRISPRSCope import cli


def _write_minimal_settings(settings_dir, extra_lines=None):
	(settings_dir / "data").mkdir()
	(settings_dir / "refs").mkdir()
	(settings_dir / "inputs").mkdir()
	(settings_dir / "data" / "r1.fastq").write_text("\n")
	(settings_dir / "data" / "r2.fastq").write_text("\n")
	(settings_dir / "inputs" / "barcodes.txt").write_text("AAAAAAAAAAAAAAAAAA\n")
	(settings_dir / "inputs" / "amplicons.tsv").write_text("ampA\tACGTACGTACGTACGTACGT\tNA\n")
	(settings_dir / "refs" / "genome.1.bt2").write_text("\n")
	lines = [
		"r1\tdata/r1.fastq",
		"r2\tdata/r2.fastq",
		"constant1\tACGT",
		"constant2\tTGCA",
		"barcodes\tinputs/barcodes.txt",
		"amplicons\tinputs/amplicons.tsv",
		"bowtie2_index\trefs/genome",
		"output_root\tresults/run",
		"processes\t1",
	]
	if extra_lines:
		lines.extend(extra_lines)
	settings = settings_dir / "settings.txt"
	settings.write_text("\n".join(lines) + "\n")
	return settings


def _parse_settings(settings, monkeypatch):
	monkeypatch.setattr(sys, "argv", ["CRISPRSCope", str(settings)])
	monkeypatch.setattr(cli.sb, "check_output", lambda *args, **kwargs: b"ok")
	return cli.parse_settings(sys.argv)


def test_parse_settings_resolves_relative_paths_from_settings_file(tmp_path, monkeypatch):
	settings = _write_minimal_settings(
		tmp_path,
		extra_lines=[
			"primer_lookup_len\t12",
			"allowBarcodeMismatches\tFalse",
			"h5ad_output\tresults/custom.h5ad",
			"debug_rescued_reads_bam\tdebug/rescued.bam",
		],
	)

	parsed = _parse_settings(settings, monkeypatch)

	assert parsed[0] == str(tmp_path / "data" / "r1.fastq")
	assert parsed[1] == str(tmp_path / "data" / "r2.fastq")
	assert parsed[4] is False
	assert parsed[5] == str(tmp_path / "inputs" / "barcodes.txt")
	assert parsed[6] == str(tmp_path / "inputs" / "amplicons.tsv")
	assert parsed[7] == 12
	assert parsed[11] == str(tmp_path / "refs" / "genome")
	assert parsed[13] == str(tmp_path / "results" / "run")
	assert parsed[23] == str(tmp_path / "results" / "custom.h5ad")
	assert parsed[25] == str(tmp_path / "debug" / "rescued.bam")


def test_parse_settings_rejects_legacy_primer_lookup_key(tmp_path, monkeypatch):
	settings = _write_minimal_settings(tmp_path, extra_lines=["primerLookupLen\t12"])
	monkeypatch.setattr(sys, "argv", ["CRISPRSCope", str(settings)])

	with pytest.raises(ValueError, match="primerLookupLen is no longer supported"):
		cli.parse_settings(sys.argv)


def test_parse_settings_rejects_non_tab_lines(tmp_path, monkeypatch):
	settings = tmp_path / "settings.txt"
	settings.write_text("r1 data/r1.fastq\n")
	monkeypatch.setattr(sys, "argv", ["CRISPRSCope", str(settings)])

	with pytest.raises(ValueError, match="key<TAB>value"):
		cli.parse_settings(sys.argv)


def test_run_command_raises_on_nonzero_exit():
	with pytest.raises(cli.ExternalCommandError, match="return code 7"):
		cli.run_command(["bash", "-lc", "exit 7"])


def test_multi_reference_names_are_constructed_without_index_error():
	assert cli._build_input_ref_names(4) == ["Reference", "Amplicon1", "Amplicon2", "Amplicon3"]


def test_split_read_cache_validates_referenced_fastqs(tmp_path):
	info_file = tmp_path / "run.splitReads.ampliconInfo.txt"
	r1 = tmp_path / "ampA.r1.fq.gz"
	r2 = tmp_path / "ampA.r2.fq.gz"
	r1.write_text("\n")
	r2.write_text("\n")
	header = "\t".join(["name", "aln_count", "reads_r1_file", "reads_r2_file"])
	info_file.write_text(
		header + "\n"
		+ f"ampA\t1\t{r1}\t{r2}\n"
		+ f"ampB\t1\t{tmp_path / 'missing.r1.fq.gz'}\t{tmp_path / 'missing.r2.fq.gz'}\n"
	)

	cache_is_valid, amplicon_names, amplicon_information = cli._load_split_read_cache(
		str(info_file),
		str(tmp_path),
	)

	assert cache_is_valid is False
	assert amplicon_names == ["ampA", "ampB"]
	assert amplicon_information["ampA"]["reads_r1_file"] == str(r1)
	assert amplicon_information["ampA"]["reads_r2_file"] == str(r2)


@pytest.mark.parametrize("guide", ["", "NA", "none"])
def test_crispresso_command_omits_empty_or_na_guides(tmp_path, guide):
	output_root = str(tmp_path / "run")
	crispresso_dir = str(tmp_path / "run.crispresso")
	seq_dir = tmp_path / "run.seq_by_amplicon"
	run_folder = tmp_path / "run.crispresso" / "CRISPResso_on_ampA"
	seq_dir.mkdir()
	run_folder.mkdir(parents=True)
	r1 = seq_dir / "ampA.r1.fq.gz"
	r2 = seq_dir / "ampA.r2.fq.gz"
	r1.write_text("\n")
	r2.write_text("\n")
	(tmp_path / "run.crispresso" / "ampA.finished").write_text("\n")
	(run_folder / "CRISPResso2_info.json").write_text("{}\n")

	result = cli.run_crispresso_commands(
		["ampA"],
		{
			"ampA": {
				"aln_count": "1",
				"reads_r1_file": str(r1),
				"reads_r2_file": str(r2),
				"amp_seqs": "ACGTACGT",
				"guide_seq": guide,
			}
		},
		output_root,
		crispresso_dir,
		False,
		1,
		alleles=False,
	)

	assert " -g " not in result["ampA"]["crispresso_command"]


def test_example_files_are_tab_delimited():
	repo_root = os.path.dirname(os.path.dirname(__file__))
	for rel_path in ["example/example_settings.txt", "example/amplicon_file.txt"]:
		with open(os.path.join(repo_root, rel_path), "r") as handle:
			for line in handle:
				stripped = line.strip()
				if not stripped or stripped.startswith("#"):
					continue
				assert "\t" in line
