from pathlib import Path

from aligons.util import fs, maf

data_dir = Path(__file__).parent.parent / "_data"


def test_maf_block_ranges(tmp_path: Path):
    source_maf = data_dir / "ucsc.maf"
    expected_bed = data_dir / "ucsc.bed"
    tmp_maf = tmp_path / "ucsc.maf"
    fs.symlink(source_maf, tmp_maf)
    tmp_bed = maf.maf_block_ranges(tmp_maf)
    with tmp_bed.open() as obs, expected_bed.open() as exp:
        assert obs.read() == exp.read()
