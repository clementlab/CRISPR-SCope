import io

from CRISPRSCope.cli import Metrics
from CRISPRSCope import crispresso, fastq_processing


def test_create_log_str_read_denominators():
    m = Metrics()
    m.tot_reads = 200
    m.has_constant1_count = 180
    m.has_constant2_count = 160
    m.barcodes_valid_count = 150
    m.barcodes_valid_error_correction_count = 5
    m.long_enough_r1_count = 50
    m.no_adapter_read_count = 25
    m.reads_per_cell["AA"] = 10
    m.reads_per_cell["BB"] = 5

    s = m.create_log_str()

    assert "Read 200 reads" in s
    assert "50 (25.00%) have sufficiently-long R1's" in s
    assert "25 (50.00%) did not contain adapter sequences" in s
    assert "assigned to 2 cells" in s


def test_run_alignment_uses_moved_path_helpers(monkeypatch, tmp_path):
    """The multiprocessing alignment worker retains its command-format helpers."""
    calls = []

    class FakeBowtieProcess:
        def __init__(self):
            self.stdout = io.BytesIO()

        def wait(self):
            return 0

    class FakeSamtoolsProcess:
        returncode = 0

        def communicate(self):
            return b"", b""

    def fake_popen(command, **kwargs):
        calls.append((command, kwargs))
        if command[0] == "bowtie2":
            return FakeBowtieProcess()
        return FakeSamtoolsProcess()

    monkeypatch.setattr(fastq_processing.sb, "Popen", fake_popen)
    output_bam = tmp_path / "aligned.bam"

    fastq_processing.run_alignment((
        "reads_R1.fastq.gz",
        "reads_R2.fastq.gz",
        "reference",
        1,
        str(output_bam),
    ))

    assert [command[0][0] for command in calls] == ["bowtie2", "samtools"]
    assert output_bam.is_file()


def test_fastq_processing_retains_bam_command_runner_dependency():
    """The FASTQ stage uses this runner when concatenating and sorting BAM files."""
    assert fastq_processing.run_command is crispresso.run_command
