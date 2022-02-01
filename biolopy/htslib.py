import gzip
import logging
import re
import shutil
import subprocess
from collections.abc import Iterable
from pathlib import Path
from typing import Any, IO

from . import cli, fs

_log = logging.getLogger(__name__)
_dry_run = False


def main(argv: list[str] | None = None):
    import argparse

    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument(
        "-f", "--fasta", action="store_const", const="fasta", dest="format"
    )
    parser.add_argument(
        "-g", "--gff3", action="store_const", const="gff3", dest="format"
    )
    parser.add_argument("path", type=Path)
    args = parser.parse_args(argv or None)
    global _dry_run
    _dry_run = args.dry_run
    cli.logging_config(args.loglevel)
    if args.format == "fasta":
        index_fasta(args.path)
    elif args.format == "gff3":
        index_gff3(args.path)


def index_fasta(path: Path):
    outfile = create_genome_bgzip(path)
    for chromosome in fs.sorted_naturally(path.glob(r"*.chromosome.*.fa.gz")):
        print(faToTwoBit(chromosome))
    print(outfile)
    print(faidx(outfile))
    print(faToTwoBit(outfile))
    print(faSize(outfile))
    return outfile


def index_gff3(path: Path):
    outfile = create_genome_bgzip(path, "gff3")
    print(outfile)
    print(tabix(outfile))
    return outfile


def create_genome_bgzip(path: Path, ext: str = "fa"):
    files = fs.sorted_naturally(path.glob(rf"*.chromosome.*.{ext}.gz"))
    assert files
    _log.debug(str(files))
    name = files[0].name
    (outname, count) = re.subn(rf"\.chromosome\..+\.{ext}", rf".genome.{ext}", name)
    assert count == 1
    outfile = path / outname
    return bgzip(files, outfile, ext)


def bgzip(infiles: Iterable[Path], outfile: Path, format: str = "fasta"):
    if fs.is_outdated(outfile) and not _dry_run:
        with open(outfile, "wb") as fout:
            bgzip = popen("bgzip -@2", stdin=subprocess.PIPE, stdout=fout)
            assert bgzip.stdin
            for file in infiles:
                if format.startswith("gff"):
                    p = sort_clean_chromosome_gff3(file)
                    assert p.stdout
                    bgzip.stdin.write(p.stdout.read())
                else:
                    with gzip.open(file, "rb") as fin:
                        bgzip.stdin.write(fin.read())
            bgzip.communicate()
    return outfile


def faidx(bgz: Path):
    """http://www.htslib.org/doc/samtools-faidx.html"""
    outfile = bgz.with_suffix(bgz.suffix + ".fai")
    if fs.is_outdated(outfile, bgz):
        run(["samtools", "faidx", str(bgz)])
    return outfile


def tabix(bgz: Path):
    """http://www.htslib.org/doc/tabix.html"""
    outfile = bgz.with_suffix(bgz.suffix + ".tbi")
    if fs.is_outdated(outfile, bgz):
        run(["tabix", str(bgz)])
    return outfile


def faToTwoBit(fa_gz: Path):
    # TODO: separate to kent.py
    """https://github.com/ENCODE-DCC/kentUtils/"""
    outfile = fa_gz.with_suffix("").with_suffix(".2bit")
    if fs.is_outdated(outfile, fa_gz) and not _dry_run:
        with gzip.open(fa_gz, "rb") as fin:
            args = ["faToTwoBit", "stdin", str(outfile)]
            p = popen(args, stdin=subprocess.PIPE)
            assert p.stdin
            shutil.copyfileobj(fin, p.stdin)
            p.stdin.close()
        p.communicate()
    return outfile


def faSize(genome_fa_gz: Path):
    # TODO: separate to kent.py
    """https://github.com/ENCODE-DCC/kentUtils/"""
    if not genome_fa_gz.name.endswith("genome.fa.gz"):
        _log.warning(f"expecting *.genome.fa.gz: {genome_fa_gz}")
    outfile = genome_fa_gz.parent / "fasize.chrom.sizes"
    if fs.is_outdated(outfile, genome_fa_gz) and not _dry_run:
        with open(outfile, "wb") as fout:
            run(["faSize", "-detailed", str(genome_fa_gz)], stdout=fout)
    return outfile


def sort_clean_chromosome_gff3(infile: Path):
    # TODO: jbrowse2 still needs billzt/gff3sort precision?
    p1 = popen(f"zgrep -v '^#' {str(infile)}", stdout=subprocess.PIPE)
    p2 = popen("grep -v '\tchromosome\t'", stdin=p1.stdout, stdout=subprocess.PIPE)
    return popen("sort -k4,4n", stdin=p2.stdout, stdout=subprocess.PIPE)


def popen(
    args: list[str] | str,
    stdin: IO[Any] | int | None = None,
    stdout: IO[Any] | int | None = None,
):  # kwargs hinders type inference to Popen[bytes]
    (args, cmd) = cli.prepare_args(args, _dry_run)
    _log.info(cmd)
    return subprocess.Popen(args, stdin=stdin, stdout=stdout)


def run(
    args: list[str] | str,
    stdin: IO[Any] | int | None = None,
    stdout: IO[Any] | int | None = None,
):
    (args, cmd) = cli.prepare_args(args, _dry_run)
    _log.info(cmd)
    return subprocess.run(args, stdin=stdin, stdout=stdout)


if __name__ == "__main__":
    main()
