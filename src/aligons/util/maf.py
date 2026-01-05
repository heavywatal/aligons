"""Utilities for MAF files."""

import io
import logging
import re
from pathlib import Path
from typing import Any

import polars as pl

from aligons.db import api, phylo
from aligons.extern import htslib

from . import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    """CLI for manual execution and testing."""
    parser = cli.ArgumentParser()
    parser.add_argument("-F", "--flank", type=int, default=0)
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    bed = maf_block_ranges(args.infile)
    fa = subseqs_from_bed(bed, args.flank)
    fs.print_if_exists(fa)


def subseqs_from_bed(infile: Path, flank: int = 0) -> Path:
    """Extract subsequences from genome FASTA based on BED file.

    :param infile: Input BED file.
    :param flank: Number of flanking bases to include on each side.
    :returns: Generated FASTA file.
    """
    outfile = infile.with_suffix(".fa")
    with subp.open_(outfile, "wb", if_=fs.is_outdated(outfile, infile)) as fout:
        bed = read_bed(infile)
        regions = to_one_based_inclusive(bed).with_columns(
            start=pl.col("start") - flank, end=pl.col("end") + flank
        )
        for row in regions.collect().iter_rows(named=True):
            fout.write(subseq(**row))
    return outfile


def subseq(chrom: str, start: int, end: int, name: str, strand: str, **_: Any) -> bytes:
    """Extract a subsequence from the genome FASTA.

    :param chrom: Chromosome name.
    :param start: Start position (1-based).
    :param end: End position (inclusive).
    :param name: Species name.
    :param strand: Strand, "+" or "-".
    :returns: Extracted sequence.
    """
    bgz = api.genome_fa(name)
    region = f"{chrom}:{start}-{end}"
    p = htslib.popen_faidx_query(bgz, region, strand=strand)
    return p.communicate()[0]


def maf_block_ranges(infile: Path) -> Path:
    """Extract alignment block ranges from a MAF file and save as BED.

    :param infile: Input MAF file.
    :returns: Generated BED file.
    """
    outfile = infile.with_suffix(".bed")
    if fs.is_outdated(outfile, infile):
        maf_s = read_s(infile)
        _log.debug(maf_s.collect().write_csv(separator="\t", include_header=False))
        bed = to_bed(maf_s)
        bed.collect().write_csv(outfile, separator="\t", include_header=False)
    return outfile


def to_one_based_inclusive(bed: pl.LazyFrame) -> pl.LazyFrame:
    """Convert BED to one-based inclusive coordinates.

    :param bed: LazyFrame from `read_bed()`, zero-based, half-closed-half-open.
    :returns: LazyFrame in one-based inclusive coordinates.
    """
    return bed.with_columns(
        start=pl.when(pl.col("strand") == "-")
        .then(pl.col("start"))
        .otherwise(pl.col("start") + 1),
        end=pl.when(pl.col("strand") == "-")
        .then(pl.col("end") - 1)
        .otherwise(pl.col("end")),
    )


def read_bed(file: Path) -> pl.LazyFrame:
    """Read a BED file.

    Note: zero-based, half-closed-half-open.
    E.g., the first 100 bases: start=0, end=100.

    :param file: Input BED file.
    :returns: LazyFrame with columns: chrom, start, end, name, score, strand.
    """
    return pl.read_csv(
        file,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end", "name", "score", "strand"],
    ).lazy()


def to_bed(maf_s: pl.LazyFrame) -> pl.LazyFrame:
    """Convert LazyFrame of MAF "s" lines to BED format.

    :param maf_s: LazyFrame from `read_s()`.
    :returns: LazyFrame in BED format with columns:
        chrom, start, end, name, score, strand.
    """
    return (
        maf_s.with_columns(pl.col("seqid").str.split_exact(".", 1))
        .unnest("seqid")
        .rename({"field_0": "species", "field_1": "chrom"})
        .with_columns(
            start=pl.when(pl.col("strand") == "-")
            .then(pl.col("fasize") - pl.col("start") - pl.col("size") + 1)
            .otherwise(pl.col("start"))
        )
        .with_columns(end=(pl.col("start") + pl.col("size")))
        .with_columns(
            name=pl.col("species").map_elements(phylo.lengthen, return_dtype=pl.String)
        )
        .with_columns(score=pl.lit("."))
        .select(["chrom", "start", "end", "name", "score", "strand"])
    )


def read_s(file: Path) -> pl.LazyFrame:
    """Read lines starting with "s " from a MAF file.

    Note: "This is a zero-based number".

    <https://genome.ucsc.edu/FAQ/FAQformat.html#format5>

    :param file: Input MAF file.
    :returns: LazyFrame with columns: seqid, start, size, strand, fasize.
    """
    source = io.BytesIO()
    with file.open("rb") as fin:
        for line in fin:
            if line.startswith(b"s "):
                source.write(re.sub(rb"[ \t]+", b"\t", line))
    return pl.read_csv(
        source,
        separator="\t",
        has_header=False,
        columns=range(1, 6),
        new_columns=["seqid", "start", "size", "strand", "fasize"],
    ).lazy()


def is_empty(file: Path) -> bool:
    """Check if a MAF file has no alignment blocks.

    :param file: Input MAF file.
    :returns: True if the file has no alignment blocks.
    """
    with file.open("rb") as fin:
        for line in fin:
            if line.startswith(b"a "):
                return False
    return True


if __name__ == "__main__":
    main()
