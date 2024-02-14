import logging
from pathlib import Path

import pytest
from aligons.util import gff

data_dir = Path(__file__).parent.parent / "_data"


def test_class(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    source_gff = data_dir / "genome.gff3"
    tmp_gff = tmp_path / "genome.gff3"
    obj = gff.GFF(source_gff)
    assert "invalid" not in caplog.text
    with tmp_gff.open("wb") as fout:
        obj.write(fout)
    with tmp_gff.open("rb") as obs, source_gff.open("rb") as exp:
        assert obs.read() == exp.read()


def test_sort(caplog: pytest.LogCaptureFixture):
    source_gff = data_dir / "genome.gff3"
    sorted_gff = data_dir / "sorted.gff3"
    obj = gff.GFF(source_gff)
    obj.sanitize()
    assert not caplog.text
    with sorted_gff.open("rt") as fin:
        assert obj.to_string() == fin.read()
    if False:
        import re

        with sorted_gff.open("wb") as fout:
            obj.write(fout)
        with (data_dir / "invalid.gff3").open("wt") as fout:
            fout.write("# comment to ignore\n")
            fout.write(re.sub(r"#.+\n", "", obj.to_string()))


def test_invalid(caplog: pytest.LogCaptureFixture):
    caplog.set_level(logging.DEBUG)
    source_gff = data_dir / "invalid.gff3"
    obj = gff.GFF(source_gff)
    assert "invalid GFF without ##gff-version" in caplog.text
    assert obj.header[0] == b"##gff-version 3\n"

    caplog.set_level(logging.DEBUG)
    source_gff = data_dir / "invalid.gff3"
    obj = gff.GFF(source_gff)
    obj.sanitize()
    assert "invalid GFF without ##gff-version" in caplog.text
    assert obj.to_string().startswith("##gff-version 3")
