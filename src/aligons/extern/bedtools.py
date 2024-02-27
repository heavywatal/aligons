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
    fi = subp.run_zcat(args.infile).stdout
    _log.info(maskfasta(fi, args.bed).decode())


def wait_maskfasta(
    fi: Path, fts: list[cli.Future[Path]], fo: Path, *, soft: bool = True
) -> Path:
    fs.expect_suffix(fo, ".gz")
    beds = [f.result() for f in fts]
    if fs.is_outdated(fo, [fi, *beds]) or cli.dry_run:
        content = subp.run_zcat(fi).stdout
        for bed in beds:
            content = maskfasta(content, bed, soft=soft)
        htslib.bgzip(content, fo)
    return fs.print_if_exists(fo)


def maskfasta(fi: bytes, bed: Path, *, soft: bool = True) -> bytes:
    """https://bedtools.readthedocs.io/en/latest/content/tools/maskfasta.html."""
    args: subp.Args = ["bedtools", "maskfasta", "-fullHeader"]
    if soft:
        args.append("-soft")
    args.extend(["-fi", "/dev/stdin"])
    args.extend(["-bed", bed])  # <BED/GFF/VCF>(.gz)
    args.extend(["-fo", "/dev/stdout"])
    p = subp.run(args, input=fi, stdout=subp.PIPE)
    return p.stdout


def subtract(a: Path, b: Path) -> bytes:
    sub: subp.Args = ["bedtools", "subtract", "-a", a, "-b", b]
    return subp.run(sub, stdout=subp.PIPE).stdout


def remove_short(b: bytes, min_width: int = 0) -> bytes:
    awk = ["awk", f"($3-$2) >= {min_width}"]
    return subp.run(awk, input=b, stdout=subp.PIPE).stdout


if __name__ == "__main__":
    main()
