import gzip

import pytest

from CRISPRSCope.io_utils import open_text_maybe_gzip


@pytest.mark.parametrize("compressed", [False, True])
def test_open_text_maybe_gzip_uses_file_content_not_suffix(tmp_path, compressed):
    path = tmp_path / "CRISPResso_output.fastq.gz"
    contents = "@read\nACGT\n+\nIIII\n"
    if compressed:
        with gzip.open(path, "wt") as handle:
            handle.write(contents)
    else:
        path.write_text(contents)

    with open_text_maybe_gzip(path) as handle:
        assert handle.read() == contents


def test_open_text_maybe_gzip_rejects_non_text_mode(tmp_path):
    path = tmp_path / "reads.fastq"
    path.write_text("@read\n")

    with pytest.raises(ValueError, match="Unsupported mode"):
        open_text_maybe_gzip(path, "rb")
