import gzip
import logging
from pathlib import Path

import pytest
from aligons.util import gff

data_dir = Path(__file__).parent / "data"


def test_split(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    source_gff = data_dir / "genome.gff3"
    genome_gff = tmp_path / "genome.gff3.gz"
    with source_gff.open("rb") as fin, gzip.open(genome_gff, "wb") as fout:
        fout.write(fin.read())
    gff.split_by_seqid(genome_gff)
    assert not caplog.text
    ls = list(tmp_path.glob("genome.chromosome.*"))
    assert len(ls) == 4  # noqa: PLR2004
    chr1_gff_res = tmp_path / "genome.chromosome.1.gff3.gz"
    chr1_gff_exp = data_dir / "genome.chromosome.1.gff3"
    with gzip.open(chr1_gff_res, "rb") as fin:
        res = fin.read()
    with chr1_gff_exp.open("rb") as fin:
        exp = fin.read()
    assert res == exp
    if False:
        with gzip.open(chr1_gff_res, "rb") as fin, (
            data_dir / "genome.chromosome.1.gff3"
        ).open("wb") as fout:
            fout.write(fin.read())


def test_sort(caplog: pytest.LogCaptureFixture):
    source_gff = data_dir / "genome.gff3"
    sorted_gff = data_dir / "sorted.gff3"
    with source_gff.open("rb") as fin:
        content = fin.read()
    res = gff.sort(content)
    assert not caplog.text
    with sorted_gff.open("rb") as fin:
        assert res == fin.read()
    if False:
        import re

        with sorted_gff.open("wb") as fout:
            fout.write(res)
        with (data_dir / "invalid.gff3").open("wb") as fout:
            fout.write(b"# comment to ignore\n")
            fout.write(re.sub(rb"#.+\n", b"", res))


def test_invalid(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    caplog.set_level(logging.DEBUG)
    source_gff = data_dir / "invalid.gff3"
    invalid_gff = tmp_path / "invalid.gff3.gz"
    with source_gff.open("rb") as fin, gzip.open(invalid_gff, "wb") as fout:
        fout.write(fin.read())
    gff.main(["--dry-run", str(invalid_gff)])
    assert "invalid GFF without ##gff-version" in caplog.text
    assert "unfriendly" in caplog.text
    assert "ignoring comments" in caplog.text
    assert "ignoring scaffold" in caplog.text

    with source_gff.open("rb") as fin:
        content = fin.read()
    res = gff.sort(content)
    assert "invalid GFF without ##gff-version" in caplog.text
    assert b"##gff-version 3" in res
