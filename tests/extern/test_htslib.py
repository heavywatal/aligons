import logging
from pathlib import Path

import pytest
from aligons.extern import htslib

data_dir = Path(__file__).parent.parent / "_data"


def test_bgzip(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    caplog.set_level(logging.INFO)
    source_gff = data_dir / "genome.gff3"
    sorted_gff = data_dir / "sorted.gff3"
    genome_gff = tmp_path / "genome.gff3.gz"
    with source_gff.open("rb") as fin:
        htslib.bgzip(fin, genome_gff)
    with sorted_gff.open("rb") as fin:
        htslib.bgzip(fin, genome_gff)
        assert "overwriting" in caplog.text
    idx = htslib.try_index(genome_gff)
    assert idx.exists()
    assert idx.suffix == ".csi"
