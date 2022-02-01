#!/usr/bin/env python3
"""Integrate chromosome wigs into a genome-wide bigwig

src: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz
dst: ./multiple/{target}/{clade}/phastcons.bw

https://github.com/ucscGenomeBrowser/kent
"""
import fileinput
import gzip
import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, IO

from . import cli, fs
from .db import ensemblgenomes

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    import argparse

    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("clade", type=Path)
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    outfile = integrate_wigs(args.clade)
    print(outfile)
    _log.info(bigWigInfo(outfile).stdout)


def integrate_wigs(clade: Path):
    species = clade.parent.name
    chrom_sizes = ensemblgenomes.get_file("*chrom.sizes", species)
    name = "phastcons.wig.gz"
    wigs = [p / name for p in fs.sorted_naturally(clade.glob("chromosome.*"))]
    _log.info(f"{[str(x) for x in wigs]}")
    outfile = clade / "phastcons.bw"
    with fileinput.input(wigs, mode="rb", openhook=fileinput.hook_compressed) as reader:
        wigToBigWig(reader, chrom_sizes, outfile)
    return outfile


def wigToBigWig(
    infile: IO[Any] | fileinput.FileInput[Any], chrom_sizes: Path, outfile: Path
):
    args = ["wigToBigWig", "stdin", str(chrom_sizes), str(outfile)]
    (args, cmd) = cli.prepare_args(args, not outfile.exists())
    _log.info(cmd)
    p = subprocess.Popen(args, stdin=subprocess.PIPE)
    if not outfile.exists() and not cli.dry_run:
        assert p.stdin
        for line in infile:
            p.stdin.write(line)
    return p.communicate()


def bigWigInfo(path: Path):
    args = ["bigWigInfo", str(path)]
    return subprocess.run(args, stdout=subprocess.PIPE, text=True)


def faToTwoBit(fa_gz: Path):
    outfile = fa_gz.with_suffix("").with_suffix(".2bit")
    if fs.is_outdated(outfile, fa_gz) and not cli.dry_run:
        with gzip.open(fa_gz, "rb") as fin:
            args = ["faToTwoBit", "stdin", str(outfile)]
            p = popen(args, stdin=subprocess.PIPE)
            assert p.stdin
            shutil.copyfileobj(fin, p.stdin)
            p.stdin.close()
        p.communicate()
    return outfile


def faSize(genome_fa_gz: Path):
    if not genome_fa_gz.name.endswith("genome.fa.gz"):
        _log.warning(f"expecting *.genome.fa.gz: {genome_fa_gz}")
    outfile = genome_fa_gz.parent / "fasize.chrom.sizes"
    if fs.is_outdated(outfile, genome_fa_gz) and not cli.dry_run:
        with open(outfile, "wb") as fout:
            run(["faSize", "-detailed", str(genome_fa_gz)], stdout=fout)
    return outfile


def popen(
    args: list[str] | str,
    stdin: IO[Any] | int | None = None,
    stdout: IO[Any] | int | None = None,
):  # kwargs hinders type inference to Popen[bytes]
    (args, cmd) = cli.prepare_args(args)
    _log.info(cmd)
    return subprocess.Popen(args, stdin=stdin, stdout=stdout)


def run(
    args: list[str] | str,
    stdin: IO[Any] | int | None = None,
    stdout: IO[Any] | int | None = None,
):
    (args, cmd) = cli.prepare_args(args)
    _log.info(cmd)
    return subprocess.run(args, stdin=stdin, stdout=stdout)


if __name__ == "__main__":
    main()
