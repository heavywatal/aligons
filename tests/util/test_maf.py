from pathlib import Path

from polars.testing import assert_series_equal

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


def test_read_bed_and_to_one():
    expected_bed = data_dir / "ucsc.bed"
    lf0 = maf.read_bed(expected_bed)
    lf1 = maf.to_one_based_inclusive(lf0)
    df0 = lf0.collect()
    df1 = lf1.collect()
    # depends on strand, but currently it is "+" for all records
    assert_series_equal(df0.get_column("chrom"), df1.get_column("chrom"))
    assert_series_equal(df0.get_column("start") + 1, df1.get_column("start"))
    assert_series_equal(df0.get_column("end"), df1.get_column("end"))
    assert_series_equal(df0.get_column("name"), df1.get_column("name"))
    assert_series_equal(df0.get_column("score"), df1.get_column("score"))
    assert_series_equal(df0.get_column("strand"), df1.get_column("strand"))
