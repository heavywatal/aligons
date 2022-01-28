#!/usr/bin/env python3
"""Integrate chromosome wigs into a genome-wide bigwig

src: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz
dst: ./multiple/{target}/{clade}/phastcons.bw

https://github.com/ENCODE-DCC/kentUtils/blob/master/src/utils/wigToBigWig/wigToBigWig.c
"""
import fileinput
import logging
import os
import subprocess
from pathlib import Path
from typing import Any, IO

from . import cli
from .db import ensemblgenomes


_log = logging.getLogger(__name__)
_dry_run = False


def main(argv: list[str] = []):
    import argparse

    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("clade", type=Path)
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    global _dry_run
    _dry_run = args.dry_run
    outfile = integrate_wigs(args.clade)
    print(outfile)
    _log.info(bigWigInfo(outfile).stdout)


def integrate_wigs(clade: Path):
    species = clade.parent.name
    chrom_sizes = ensemblgenomes.get_file("*chrom.sizes", species)
    wigs = list(sorted_wigs(clade, chrom_sizes))
    _log.info(f"{[str(x) for x in wigs]}")
    outfile = clade / "phastcons.bw"
    with fileinput.input(wigs, mode="rb", openhook=fileinput.hook_compressed) as reader:
        wigToBigWig(reader, chrom_sizes, outfile)
    return outfile


def wigToBigWig(
    infile: IO[Any] | fileinput.FileInput[Any], chrom_sizes: Path, outfile: Path
):
    args = ["wigToBigWig", "stdin", str(chrom_sizes), str(outfile)]
    (args, cmd) = cli.prepare_args(args, outfile.exists() or _dry_run)
    _log.info(cmd)
    p = subprocess.Popen(args, stdin=subprocess.PIPE)
    if not outfile.exists() and not _dry_run:
        assert p.stdin
        for line in infile:
            p.stdin.write(line)
    return p.communicate()


def bigWigInfo(path: Path):
    args = ["bigWigInfo", str(path)]
    return subprocess.run(args, stdout=subprocess.PIPE, text=True)


def sorted_wigs(path: Path, chrom_sizes: Path):
    with open(chrom_sizes, "r") as fin:
        for line in fin:
            seqname = line.split("\t")[0]
            yield path / f"chromosome.{seqname}" / "phastcons.wig.gz"


if __name__ == "__main__":
    main()
