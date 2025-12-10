"""HTSlib: C library for reading/writing high-throughput sequencing data.

<https://www.htslib.org/>
"""

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


def concat_bgzip(infiles: list[Path], outfile: Path) -> Path:
    """Concatenate files into a bgzipped file.

    :param infiles: A list of files to concatenate, can be gzipped.
    :param outfile: A bgzipped output file.
    :returns: The same path as `outfile`.
    """
    if fs.is_outdated(outfile, infiles) and not cli.dry_run:
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        with popen_bgzip(outfile) as bgz:
            assert bgz.stdin is not None
            subp.run_zcat(infiles, stdout=bgz.stdin)
    return fs.print_if_exists(outfile)


def bgzip(data: bytes | IO[bytes] | None, outfile: Path, *, if_: bool = True) -> Path:
    """Compress data in BGZF format.

    <https://www.htslib.org/doc/bgzip.html>

    :param data: Bytes or a binary stream to compress.
    :param outfile: Compressed BGZF file.
    """
    fs.expect_suffix(outfile, ".gz")
    if outfile.exists():
        _log.info(f"overwriting {outfile}")
    cmd = ["bgzip", "-@2", "-o", outfile]
    if isinstance(data, bytes):
        subp.run(cmd, input=data, if_=if_)
    else:
        subp.run(cmd, stdin=data, if_=if_)
    return outfile


def popen_bgzip(
    outfile: Path, *, stdin: subp.FILE = subp.PIPE, if_: bool = True
) -> subp.Popen[bytes]:
    """Open a BGZF compression process for writing.

    :param outfile: Compressed BGZF file.
    :param stdin: Passed to `subprocess.Popen`.
    :returns: A `subprocess.Popen` object.
    """
    fs.expect_suffix(outfile, ".gz")
    if outfile.exists():
        _log.info(f"overwriting {outfile}")
    cmd = ["bgzip", "-@2", "-o", outfile]
    return subp.popen(cmd, stdin=stdin, if_=if_)


def try_index(bgz: Path | cli.Future[Path]) -> Path:
    """Call `faidx` or `tabix` on the given BGZF file.

    :param bgz: BGZF-compressed file or a future of it.
    :returns: The index file if created, else the original file.
    """
    bgz = cli.result(bgz)
    if to_be_tabixed(bgz.name):
        return tabix(bgz)
    if to_be_faidxed(bgz.name):
        return faidx(bgz)
    return bgz


def read_fai(genome: Path) -> dict[str, int]:
    """Read a FASTA index (`.fai`) file.

    :param genome: A bgzipped genome FASTA. Use `api.genome_fa(species)`.
    :returns: A mapping of sequence IDs to their lengths.
    """
    fai = faidx(genome)
    res: dict[str, int] = {}
    with fai.open("rt") as fin:
        for line in fin:
            (seqid, size, _rest) = line.split(maxsplit=2)
            res[seqid] = int(size)
    return res


def faidx(bgz: Path | cli.Future[Path]) -> Path:
    """Index a bgzipped FASTA file.

    <https://www.htslib.org/doc/samtools-faidx.html>

    :param bgz: A bgzipped FASTA file or a future of it.
    :returns: The FASTA index file (`.fai`).
    """
    bgz = cli.result(bgz)
    fs.expect_suffix(bgz, ".gz")
    outfile = bgz.with_suffix(bgz.suffix + ".fai")
    subp.run(["samtools", "faidx", bgz], if_=fs.is_outdated(outfile, bgz))
    return fs.print_if_exists(outfile)


def faidx_query(bgz: Path | cli.Future[Path], region: str, outfile: Path) -> Path:
    """Query a region from an indexed bgzipped FASTA file.

    :param bgz: A bgzipped FASTA file or a future of it.
    :param region: The region to query like `seqid:start-end`.
    :param outfile: The output FASTA for the queried region.
    """
    bgz = cli.result(bgz)
    fs.expect_suffix(bgz, ".gz")
    args = ["samtools", "faidx", bgz, region, "-o", outfile]
    subp.run(args, if_=fs.is_outdated(outfile, bgz))
    return outfile


def popen_faidx_query(
    bgz: Path | cli.Future[Path],
    region: str,
    *,
    strand: str = "+",
    stdout: subp.FILE = subp.PIPE,
    if_: bool = True,
) -> subp.Popen[bytes]:
    """Query a region from an indexed bgzipped FASTA file.

    Note: 1-based, inclusive coordinates, e.g., the first 100 bases: 1-100.

    :param bgz: A bgzipped FASTA file or a future of it.
    :param region: The region to query like `seqid:start-end`.
    :param strand: "+" or "-".
    :param stdout: Passed to `subprocess.Popen`.
    :returns: A `subprocess.Popen` object.
    """
    bgz = cli.result(bgz)
    fs.expect_suffix(bgz, ".gz")
    args: subp.Args = ["samtools", "faidx", bgz, region]
    if strand == "-":
        args.append("-i")
    return subp.popen(args, stdout=stdout, if_=if_)


def tabix(bgz: Path | cli.Future[Path]) -> Path:
    """Index a bgzipped TAB-delimited file.

    Use `.csi` instead of `.tbi` for chromosomes >512 Mbp.
    <https://www.htslib.org/doc/tabix.html>

    :param bgz: A bgzipped TAB-delimited file or a future of it.
    :returns: The index file (`.csi`).
    """
    bgz = cli.result(bgz)
    outfile = bgz.with_suffix(bgz.suffix + ".csi")
    subp.run(["tabix", "--csi", bgz], if_=fs.is_outdated(outfile, bgz))
    return fs.print_if_exists(outfile)


def index(cram: Path | cli.Future[Path]) -> Path:
    """Index a CRAM file.

    <https://www.htslib.org/doc/samtools-index.html>

    :param cram: A CRAM file or a future of it.
    :returns: The CRAM index file (`.crai`).
    """
    cram = cli.result(cram)
    outfile = cram.with_suffix(cram.suffix + ".crai")
    subp.run(["samtools", "index", cram], if_=fs.is_outdated(outfile, cram))
    return fs.print_if_exists(outfile)


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
