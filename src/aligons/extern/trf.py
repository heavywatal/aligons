"""TRF: Tandem Repeats Finder.

<https://tandem.bu.edu/trf/home>

src: {basename}.fa
dat: {basename}.fa.2.5.7.80.10.40.500.dat.gz
out: {basename}.fa.trf.bed.gz
"""

import logging
import re
from collections.abc import Iterator
from pathlib import Path

import polars as pl

from aligons.util import cli, fs, subp

from . import htslib

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
        dat_to_bed(dat, outfile)
    return fs.print_if_exists(outfile)


def trf(infile: Path) -> Path:
    """Be careful of messy and dirty output from trf.

    - without `-ngs`
      - returns non-zero even when it is "Done" successfully.
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
        "min_score": 40,
        "max_period": 500,
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
    return fs.print_if_exists(dat)


def dat_to_bed(dat: Path, bed: Path) -> Path:
    content = b"".join(_dat_to_bed_iter(dat))
    return htslib.bgzip(content, bed)


def _dat_to_bed_iter(infile: Path) -> Iterator[bytes]:
    content = subp.run_zcat(infile).stdout
    for block in content.split(b"\n\nSequence: ")[1:]:
        yield _block_to_bed(block).encode()


def _block_to_bed(block: bytes) -> str:
    mobj = re.search(rb"\S+", block)
    assert mobj, block
    seqid = mobj.group(0).decode()
    _log.debug(seqid)
    return (
        _read_dat_body(block)
        .with_columns(seqid=pl.lit(seqid), name=pl.lit("."), strand=pl.lit("."))
        .select(["seqid", "start", "end", "name", "score", "strand"])
        .sort("start")
        .collect()
        .write_csv(separator="\t", include_header=False)
    )


def _read_dat_body(source: Path | str | bytes) -> pl.LazyFrame:
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
    ).lazy()


if __name__ == "__main__":
    main()
