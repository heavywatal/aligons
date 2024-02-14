import gzip
import logging
import re
from pathlib import Path

import pytest
from aligons.extern import htslib

data_dir = Path(__file__).parent.parent / "_data"


def test_split_cat(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    source_gff = data_dir / "genome.gff3"
    genome_gff = tmp_path / "genome.gff3.gz"
    sorted_gff = data_dir / "sorted.gff3"
    fasize = data_dir / "fasize.chrom.sizes"
    genome_gff.symlink_to(source_gff)
    chr_gffs = htslib.split_gff3(genome_gff, fasize)
    assert not caplog.text
    ls = list(tmp_path.glob("genome.chromosome.*"))
    assert len(ls) == 4  # noqa: PLR2004
    chr1_gff_res = tmp_path / "genome.chromosome.1.gff3.gz"
    chr1_gff_exp = data_dir / "genome.chromosome.1.gff3"
    with gzip.open(chr1_gff_res, "rt") as obs, chr1_gff_exp.open("rt") as exp:
        assert obs.read() == exp.read()

    cat_gff = genome_gff.with_name("cat.gff3.gz")
    htslib.concat_bgzip(chr_gffs, cat_gff)
    with gzip.open(cat_gff, "rt") as obs, sorted_gff.open("rt") as exp:
        str_exp = re.sub(r"\bscaffold.+?\n", "", exp.read())
        assert obs.read() == str_exp
    if False:
        with gzip.open(chr1_gff_res, "rb") as fin, (
            data_dir / "genome.chromosome.1.gff3"
        ).open("wb") as fout:
            fout.write(fin.read())


def test_invalid(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    caplog.set_level(logging.DEBUG)
    source_gff = data_dir / "invalid.gff3"
    invalid_gff = tmp_path / "invalid.gff3.gz"
    fasize = data_dir / "fasize.chrom.sizes"
    invalid_gff.symlink_to(source_gff)
    htslib.split_gff3(invalid_gff, fasize)
    assert "invalid GFF without ##gff-version" in caplog.text
    assert "unfriendly" in caplog.text
    assert "ignoring comments" in caplog.text
