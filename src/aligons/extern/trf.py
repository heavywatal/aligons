"""TRF: Tandem Repeats Finder.

https://tandem.bu.edu/trf/home

src: {basename}.fa
dat: {basename}.fa.2.5.7.80.10.40.500.dat.gz
out: {basename}.fa.trf.bed.gz
"""
import logging
import re
from pathlib import Path

import polars as pl

from aligons.extern import htslib
from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path, nargs="+")
    args = parser.parse_args(argv or None)
    fts = [cli.thread_submit(run, f) for f in args.infile]
    cli.wait_raise(fts)


def run(infile: Path) -> Path:
    dat = trf(infile)
    outfile = infile.with_suffix(infile.suffix + ".trf.bed.gz")
    if fs.is_outdated(outfile, [infile, dat]) and not cli.dry_run:
        bed = dat_to_bed(dat)
        htslib.bgzip(bed.encode(), outfile)
    _log.info(f"{outfile}")
    return outfile


def trf(infile: Path) -> Path:
    """Be careful of messy and dirty output from trf.

    - without `-ngs`
      - returns non-zero even when it is "Done" successsfully.
      - writes ./{infile}.{params}.dat including verbose header.
    - with `-ngs`
      - returns 0 if successful.
      - writes to stdout in a different format.
        - no header; concise sequence info with leading @.
        - has two extra columns with DNAs, which bloats the file ~50%.

    https://github.com/Benson-Genomics-Lab/TRF/tree/master/src
    """
    fs.expect_suffix(infile, ".gz", negate=True)
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
    dat = infile.with_name(".".join([infile.name, *params, "dat.gz"]))
    args: subp.Args = ["trf", infile]
    args.extend(params)
    args.extend(["-d", "-h", "-l", "10"])
    is_to_run = fs.is_outdated(dat, infile) and not cli.dry_run
    p = subp.run(args, if_=is_to_run, check=False)
    if is_to_run:
        pwd_dat = Path(dat.name.removesuffix(".gz"))
        _log.info(f"trf returned {p.returncode}; wrote {pwd_dat}")
        with pwd_dat.open("rb") as fin:
            subp.gzip(fin, dat)
        pwd_dat.unlink()
    _log.info(f"{dat}")
    return dat


def dat_to_bed(infile: Path) -> str:
    with infile.open("rb") as fin:
        content = fs.gzip_decompress(fin.read())
    blocks = content.split(b"\n\nSequence: ")[1:]
    return "".join([_block_to_bed(x) for x in blocks])


def _block_to_bed(block: bytes) -> str:
    assert (mobj := re.search(rb"\S+", block))
    seqid = mobj.group(0).decode()
    _log.debug(seqid)
    return (
        _read_dat_body(block)
        .with_columns(
            pl.lit(seqid).alias("seqid"),
            pl.lit(".").alias("name"),
            pl.lit(".").alias("strand"),
        )
        .select(["seqid", "start", "end", "name", "score", "strand"])
        .sort("start")
        .write_csv(separator="\t", has_header=False)
    )


def _read_dat_body(source: Path | str | bytes) -> pl.DataFrame:
    return pl.read_csv(
        source,
        has_header=False,
        separator=" ",
        skip_rows=7,
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


if __name__ == "__main__":
    main()
