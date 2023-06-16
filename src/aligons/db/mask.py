"""Softmask DNA sequences.

src: {stem}.dna.chromosome.{seqid}.fa
dst: {stem}.dna_sm.chromosome.{seqid}.fa.gz
"""
import logging
import re
from pathlib import Path

from aligons.extern import bedtools, repeatmasker, sdust, trf
from aligons.util import cli, fs

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-S", "--species")
    parser.add_argument("infile", type=Path, nargs="+")
    args = parser.parse_args(argv or None)
    cli.wait_raise(run(x, args.species) for x in args.infile)


def run(
    infile: Path, species: str | None = None, outfile: Path | None = None
) -> cli.FuturePath:
    assert infile.suffix != ".gz"
    if not species:
        species = "angiosperms"
    if outfile is None:
        patt = re.compile(r"\.dna\.")
        assert patt.search(infile.name)
        outfile = infile.parent / (patt.sub(".dna_sm.", infile.name) + ".gz")
    assert outfile.suffix == ".gz"
    fts: list[cli.FuturePath] = []
    fts.append(cli.thread_submit(repeatmasker.repeatmasker, infile, species))
    fts.append(cli.thread_submit(sdust.run, infile))
    fts.append(cli.thread_submit(trf.run, infile))
    if fs.is_outdated(outfile, infile) and not cli.dry_run:
        with infile.open("rb") as fin:
            fi = fin.read()
    else:
        fi = b""
    return cli.thread_submit(bedtools.wait_maskfasta, fi, fts, outfile)


if __name__ == "__main__":
    main()
