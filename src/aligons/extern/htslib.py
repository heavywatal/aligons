import logging
from pathlib import Path
from typing import IO

from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    try_index(args.infile)


def faidx_query(bgz: Path, region: str, *, if_: bool = True) -> subp.Popen[bytes]:
    args: subp.Args = ["samtools", "faidx", bgz, region]
    return subp.popen(args, stdout=subp.PIPE, if_=if_)


def concat_bgzip(infiles: list[Path], outfile: Path) -> Path:
    if fs.is_outdated(outfile, infiles) and not cli.dry_run:
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        with popen_bgzip(outfile) as bgz:
            assert bgz.stdin is not None
            for infile in infiles:
                subp.run_zcat(infile, bgz.stdin)
    _log.info(f"{outfile}")
    return outfile


def bgzip(data: bytes | IO[bytes] | None, outfile: Path, *, if_: bool = True) -> Path:
    """https://www.htslib.org/doc/bgzip.html."""
    fs.expect_suffix(outfile, ".gz")
    if outfile.exists():
        _log.info(f"overwriting {outfile}")
    if_ = if_ and bool(data) and not cli.dry_run
    with subp.open_(outfile, "wb", if_=if_) as fout:
        if isinstance(data, bytes):
            subp.run(["bgzip", "-@2"], input=data, stdout=fout, if_=if_)
        else:
            subp.run(["bgzip", "-@2"], stdin=data, stdout=fout, if_=if_)
    return outfile


def popen_bgzip(
    outfile: Path, *, stdin: subp.FILE = subp.PIPE, if_: bool = True
) -> subp.Popen[bytes]:
    fs.expect_suffix(outfile, ".gz")
    if outfile.exists():
        _log.info(f"overwriting {outfile}")
    if_ = if_ and not cli.dry_run
    with subp.open_(outfile, "wb", if_=if_) as fout:
        return subp.popen(["bgzip", "-@2"], stdin=stdin, stdout=fout, if_=if_)


def try_index(bgz: Path | cli.Future[Path]) -> Path:
    bgz = cli.result(bgz)
    if to_be_tabixed(bgz.name):
        return tabix(bgz)
    if to_be_faidxed(bgz.name):
        return faidx(bgz)
    return bgz


def faidx(bgz: Path | cli.Future[Path]) -> Path:
    """https://www.htslib.org/doc/samtools-faidx.html."""
    bgz = cli.result(bgz)
    outfile = bgz.with_suffix(bgz.suffix + ".fai")
    subp.run(["samtools", "faidx", bgz], if_=fs.is_outdated(outfile, bgz))
    _log.info(f"{outfile}")
    return outfile


def tabix(bgz: Path | cli.Future[Path]) -> Path:
    """https://www.htslib.org/doc/tabix.html.

    Use .csi instead of .tbi for chromosomes >512 Mbp e.g., atau, hvul.
    """
    bgz = cli.result(bgz)
    outfile = bgz.with_suffix(bgz.suffix + ".csi")
    subp.run(["tabix", "--csi", bgz], if_=fs.is_outdated(outfile, bgz))
    _log.info(f"{outfile}")
    return outfile


def index(cram: Path | cli.Future[Path]) -> Path:
    """https://www.htslib.org/doc/samtools-index.html."""
    cram = cli.result(cram)
    outfile = cram.with_suffix(cram.suffix + ".crai")
    subp.run(["samtools", "index", cram], if_=fs.is_outdated(outfile, cram))
    _log.info(f"{outfile}")
    return outfile


def to_be_bgzipped(filename: str) -> bool:
    return to_be_faidxed(filename) or to_be_tabixed(filename)


def to_be_faidxed(filename: str) -> bool:
    ext = (".fa", ".fas", ".fasta", ".fna")
    return filename.removesuffix(".gz").removesuffix(".zip").endswith(ext)


def to_be_tabixed(filename: str) -> bool:
    ext = (".gff", ".gff3", ".gtf", ".bed")
    return filename.removesuffix(".gz").removesuffix(".zip").endswith(ext)


if __name__ == "__main__":
    main()
