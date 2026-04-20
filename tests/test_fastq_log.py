from CRISPRSCope.cli import Metrics


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
