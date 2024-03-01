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


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-S", "--species")
    parser.add_argument("infile", type=Path, nargs="+")
    args = parser.parse_args(argv or None)
    cli.wait_raise(submit(x, args.species) for x in args.infile)


def submit(
    infile: cli.Future[Path] | Path,
    species: str | None = None,
    outfile: Path | None = None,
) -> cli.Future[Path]:
    infile = cli.result(infile)
    if ".dna_sm." in infile.name:
        fs.expect_suffix(infile, ".gz")
        return cli.thread_submit(cli.result, infile)
    fs.expect_suffix(infile, ".gz", negate=True)
    if not species:
        species = "angiosperms"
    if outfile is None:
        (outname, count) = re.subn(r"\.dna\.", ".dna_sm.", infile.name)
        assert count == 1, infile
        outfile = infile.with_name(outname + ".gz")
    fs.expect_suffix(outfile, ".gz")
    fts: list[cli.Future[Path]] = []
    fts.append(cli.thread_submit(repeatmasker.repeatmasker, infile, species))
    fts.append(cli.thread_submit(sdust.run, infile))
    fts.append(cli.thread_submit(trf.run, infile))
    return cli.thread_submit(bedtools.wait_maskfasta, infile, fts, outfile)


if __name__ == "__main__":
    main()
