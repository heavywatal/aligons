"""Softmask DNA sequences.

src: {stem}.dna.chromosome.{seqid}.fa.gz
dst: {stem}.dna_sm.chromosome.{seqid}.fa.gz
"""
import gzip
import logging
import re
from pathlib import Path

from aligons.extern import bedtools, repeatmasker, sdust, trf
from aligons.util import cli, fs

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-S", "--species", default="rice")
    parser.add_argument("infile", type=Path, nargs="+")
    args = parser.parse_args(argv or None)
    cli.wait_raise(run(x, args.species) for x in args.infile)


def run(infile: Path, species: str, outfile: Path | None = None) -> cli.FuturePath:
    if outfile is None:
        patt = re.compile(r"\.dna\.")
        assert patt.search(infile.name)
        outfile = infile.parent / patt.sub(".dna_sm.", infile.name)
    assert outfile.suffix == ".gz"
    fa = prepare(infile)
    fts: list[cli.FuturePath] = []
    fts.append(cli.thread_submit(repeatmasker.repeatmasker, fa, species))
    fts.append(cli.thread_submit(sdust.run, fa))
    fts.append(cli.thread_submit(trf.run, fa))
    if fs.is_outdated(outfile, infile) and not cli.dry_run:
        with fa.open("rb") as fin:
            fi = fin.read()
    else:
        fi = b""
    return cli.thread_submit(bedtools.wait_maskfasta, fi, fts, outfile)


def prepare(infile: Path):
    outdir: Path = infile.parent / "_mask"
    outfile = outdir / infile.name.removesuffix(".gz")
    if fs.is_outdated(outfile, infile) and not cli.dry_run:
        outdir.mkdir(0o755, exist_ok=True)
        if infile.suffix == ".gz":
            with gzip.open(infile, "rb") as fin, outfile.open("wb") as fout:
                fout.write(fin.read())
        else:
            outfile.symlink_to(infile)
    return outfile


if __name__ == "__main__":
    main()
