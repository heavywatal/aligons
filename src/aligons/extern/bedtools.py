"""bedtools: a powerful toolset for genome arithmetic.

https://bedtools.readthedocs.io/
"""
import logging
from pathlib import Path

from aligons.util import cli, fs, subp

from . import htslib

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    parser.add_argument("bed", type=Path)
    args = parser.parse_args(argv or None)
    with args.infile.open("rb") as fin:
        print(maskfasta(fin.read(), args.bed).decode())


def wait_maskfasta(
    fi: bytes, fts: list[cli.FuturePath], fo: Path, *, soft: bool = True
) -> Path:
    fs.expect_suffix(fo, ".gz")
    beds = [f.result() for f in fts]
    if fs.is_outdated(fo, beds) or cli.dry_run:
        for bed in beds:
            fi = maskfasta(fi, bed, soft=soft)
        htslib.bgzip(fi, fo)
    _log.info(f"{fo}")
    return fo


def maskfasta(fi: bytes, bed: Path, *, soft: bool = True) -> bytes:
    """https://bedtools.readthedocs.io/en/latest/content/tools/maskfasta.html."""
    args: subp.Args = ["bedtools", "maskfasta", "-fullHeader"]
    if soft:
        args.append("-soft")
    args.extend(["-fi", "/dev/stdin"])
    args.extend(["-bed", bed])  # <BED/GFF/VCF>(.gz)
    args.extend(["-fo", "/dev/stdout"])
    p = subp.run(args, input=fs.gzip_decompress(fi), stdout=subp.PIPE)
    return p.stdout


if __name__ == "__main__":
    main()
