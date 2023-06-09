"""RepeatMasker finds interspersed repeats and low complexity sequences.

https://www.repeatmasker.org/
https://github.com/rmhubley/RepeatMasker

src: {basename}.fa
dst: {basename}.fa.out.gff
"""
import io
import logging
import re
from pathlib import Path

import polars as pl

from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-S", "--species")
    parser.add_argument("infile", type=Path, nargs="+")
    args = parser.parse_args(argv or None)
    for infile in args.infile:
        cli.thread_submit(repeatmasker, infile, args.species)


def repeatmasker(infile: Path, species: str = "", *, soft: bool = True):
    assert infile.suffix != ".gz"
    outfile = infile.parent / (infile.name + ".out.gff")
    parallel: int = 2
    args: subp.Args = ["RepeatMasker", "-e", "rmblast", "-gff"]
    if soft:
        args.append("-xsmall")
    if parallel > 1:
        args.extend(["-pa", f"{parallel}"])  # RMBlast uses 4 cores per "-pa"
    if species:
        args.extend(["-species", species])
    args.append(infile)
    subp.run_if(fs.is_outdated(outfile, infile), args)
    _log.info(f"{outfile}")
    return outfile


def read_out(infile: Path):
    assert infile.suffix == ".out"
    with infile.open("rb") as fin:
        content = re.sub(rb" *\n *", rb"\n", fin.read())
        content = re.sub(rb" +", rb"\t", content)
    return pl.read_csv(
        io.BytesIO(content),
        has_header=False,
        separator="\t",
        skip_rows=3,
        new_columns=[
            "sw_score",
            "perc_div",
            "perc_del",
            "perc_ins",
            "seqid",
            "begin",
            "end",
            "left",
            "strand",
            "repeat",
            "class",
            "rep_begin",
            "rep_end",
            "rep_left",
            "id",
        ],
    )


if __name__ == "__main__":
    main()
