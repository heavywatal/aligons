"""TRF: Tandem Repeats Finder.

https://tandem.bu.edu/trf/home

src: {basename}.fa
dat: {basename}.fa.2.5.7.80.10.40.500.dat
out: {basename}.fa.trf.bed.gz
"""
import logging
from pathlib import Path

import polars as pl

from aligons.extern import htslib
from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path, nargs="+")
    args = parser.parse_args(argv or None)
    for infile in args.infile:
        cli.thread_submit(run, infile)


def run(infile: Path) -> Path:
    dat = trf(infile)
    bed = dat_to_bed(dat)
    outfile = infile.with_suffix(infile.suffix + ".trf.bed.gz")
    with outfile.open("wb") as fout:
        fout.write(htslib.bgzip_compress(bed.encode()))
    return outfile


def trf(infile: Path):
    assert infile.suffix != ".gz"
    param_defaults = {
        "match": 2,
        "mismatch": 5,
        "delta": 7,
        "pm": 80,
        "pi": 10,
        "minscore": 40,
        "maxperiod": 500,
    }
    params = [str(x) for x in param_defaults.values()]
    dat = infile.parent / ".".join([infile.name, *params, "dat"])
    args: subp.Args = ["trf", infile]
    args.extend(params)
    args.extend(["-d", "-h", "-l", "10"])
    subp.run_if(fs.is_outdated(dat, infile), args)
    _log.info(f"{dat}")
    return dat


def dat_to_bed(infile: Path) -> str:
    seqid = read_dat_seqid(infile)
    bed_df = (
        read_dat_body(infile)
        .with_columns(
            pl.lit(seqid).alias("seqid"),
            pl.lit(".").alias("name"),
            pl.lit(".").alias("strand"),
        )
        .select(["seqid", "start", "end", "name", "score", "strand"])
        .sort("start")
    )
    return bed_df.write_csv(separator="\t", has_header=False)


def read_dat_body(infile: Path):
    return pl.read_csv(
        infile,
        has_header=False,
        separator=" ",
        skip_rows=15,
        new_columns=[
            "start",
            "end",
            "width",
            "copies",
            "consensus",
            "matches",
            "indels",
            "score",
            "A",
            "C",
            "G",
            "T",
            "entropy",
            "seq1",
            "seq2",
        ],
    )


def read_dat_seqid(infile: Path):
    with infile.open("rt") as fin:
        for line in fin:
            if line.startswith("Sequence:"):
                return line.removeprefix("Sequence: ").split(" ", 1)[0]
    _log.warning("Sequence ID not found")
    return ""


if __name__ == "__main__":
    main()
