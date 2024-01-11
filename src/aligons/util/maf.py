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
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    print(maf_to_fa(args.infile))


def maf_to_fa(infile: Path) -> Path:
    outfile = infile.with_suffix(".fa")
    maf_s = read_s(infile)
    _log.debug(maf_s.write_csv(separator="\t", include_header=False))
    bed = to_bed(maf_s)
    _log.debug(bed.write_csv(separator="\t", include_header=False))
    with subp.open_(outfile, "wb", if_=fs.is_outdated(outfile, infile)) as fout:
        for row in bed.iter_rows(named=True):
            fout.write(subseq(**row))
    return outfile


def subseq(chrom: str, start: int, end: int, name: str, strand: str, **_: Any) -> bytes:
    bgz = api.genome_fa(name)
    region = f"{chrom}:{start}-{end}"
    args: subp.Args = ["samtools", "faidx", bgz, region]
    if strand == "-":
        args.append("-i")
    p = subp.run(args, stdout=subp.PIPE)
    return p.stdout


def to_bed(maf_s: pl.DataFrame) -> pl.DataFrame:
    return (
        maf_s.with_columns(pl.col("seqid").str.split_exact(".", 1))
        .unnest("seqid")
        .rename({"field_0": "species", "field_1": "chrom"})
        .with_columns(
            pl.when(pl.col("strand") == "-")
            .then(pl.col("fasize") - pl.col("start") - pl.col("size") + 1)
            .otherwise(pl.col("start"))
            .alias("start")
        )
        .with_columns((pl.col("start") + pl.col("size")).alias("end"))
        .with_columns(pl.col("species").map_elements(phylo.lengthen).alias("name"))
        .with_columns(score=pl.lit("."))
        .select(["chrom", "start", "end", "name", "score", "strand"])
    )


def read_s(file: Path) -> pl.DataFrame:
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
        new_columns=[
            "seqid",
            "start",
            "size",
            "strand",
            "fasize",
        ],
    )


if __name__ == "__main__":
    main()
