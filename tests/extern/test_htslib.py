import logging
from pathlib import Path

import pytest
from aligons.extern import htslib

data_dir = Path(__file__).parent.parent / "_data"


def test_fasta(tmp_path: Path):
    genome_fa = data_dir / "genome.fa"
    genome_fa_gz = tmp_path / "genome.fa.gz"
    assert htslib.to_be_bgzipped(genome_fa.name)
    assert htslib.to_be_faidxed(genome_fa.name)
    assert not htslib.to_be_tabixed(genome_fa.name)
    with genome_fa.open("rb") as fin:
        htslib.bgzip(fin, genome_fa_gz)
    assert genome_fa_gz.exists()
    idx = htslib.try_index(genome_fa_gz)
    assert idx.exists()
    assert idx.suffix == ".fai"
    chromosomes: list[Path] = []
    for seqid in ["1", "2"]:
        chr_fa = genome_fa_gz.with_name(f"chromosome.{seqid}.fa")
        chromosomes.append(chr_fa)
        htslib.faidx_query(genome_fa_gz, seqid, chr_fa)
    for seqid in ["3", "4"]:
        chr_fa = genome_fa_gz.with_name(f"chromosome.{seqid}.fa.gz")
        chromosomes.append(chr_fa)
        with (
            htslib.popen_bgzip(chr_fa) as bgzip,
            htslib.popen_faidx_query(genome_fa_gz, seqid, stdout=bgzip.stdin) as p,
        ):
            p.communicate()
    cat_fa_gz = genome_fa_gz.with_name("cat.fa.gz")
    htslib.concat_bgzip(chromosomes, cat_fa_gz)
    with cat_fa_gz.open("rb") as obs, genome_fa_gz.open("rb") as exp:
        assert obs.read() == exp.read()


def test_gff(tmp_path: Path, caplog: pytest.LogCaptureFixture):
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
