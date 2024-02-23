import io
import logging
import re
from pathlib import Path
from typing import Any

import polars as pl

from aligons.db import api, phylo

from . import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-F", "--flank", type=int, default=0)
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    bed = maf_block_ranges(args.infile)
    fa = subseqs_from_bed(bed, args.flank)
    print(fa)


def subseqs_from_bed(infile: Path, flank: int = 0) -> Path:
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
    """Note: 1-based, inclusive coordinates.

    E.g., the first 100 bases: 1-100.
    """
    bgz = api.genome_fa(name)
    region = f"{chrom}:{start}-{end}"
    args: subp.Args = ["samtools", "faidx", bgz, region]
    if strand == "-":
        args.append("-i")
    p = subp.run(args, stdout=subp.PIPE)
    return p.stdout


def maf_block_ranges(infile: Path) -> Path:
    outfile = infile.with_suffix(".bed")
    if fs.is_outdated(outfile, infile):
        maf_s = read_s(infile)
        _log.debug(maf_s.collect().write_csv(separator="\t", include_header=False))
        bed = to_bed(maf_s)
        bed.collect().write_csv(outfile, separator="\t", include_header=False)
    return outfile


def to_one_based_inclusive(bed: pl.LazyFrame) -> pl.LazyFrame:
    return bed.with_columns(
        start=pl.when(pl.col("strand") == "-")
        .then(pl.col("start"))
        .otherwise(pl.col("start") + 1),
        end=pl.when(pl.col("strand") == "-")
        .then(pl.col("end") - 1)
        .otherwise(pl.col("end")),
    )


def read_bed(file: Path) -> pl.LazyFrame:
    """Note: zero-based, half-closed-half-open.

    E.g., the first 100 bases: start=0, end=100.
    """
    return pl.read_csv(
        file,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end", "name", "score", "strand"],
    ).lazy()


def to_bed(maf_s: pl.LazyFrame) -> pl.LazyFrame:
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
        .with_columns(name=pl.col("species").map_elements(phylo.lengthen))
        .with_columns(score=pl.lit("."))
        .select(["chrom", "start", "end", "name", "score", "strand"])
    )


def read_s(file: Path) -> pl.LazyFrame:
    """Note: "This is a zero-based number".

    <https://genome.ucsc.edu/FAQ/FAQformat.html#format5>
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


if __name__ == "__main__":
    main()
