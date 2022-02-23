#!/usr/bin/env python3
"""Integrate chromosome wigs into a genome-wide bigwig

src: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz
dst: ./multiple/{target}/{clade}/phastcons.bw

https://github.com/ucscGenomeBrowser/kent
"""
import argparse
import gzip
import logging
import os
import shutil
from pathlib import Path

from .. import cli, fs
from ..db import ensemblgenomes

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("clade", type=Path)
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    run(args.clade)


def run(clade: Path):
    if (bigwig := integrate_wigs(clade)).exists():
        print(bigwig)
        _log.info(bigWigInfo(bigwig).rstrip())


def integrate_wigs(clade: Path):
    species = clade.parent.name
    chrom_sizes = ensemblgenomes.get_file("*chrom.sizes", species)
    name = "phastcons.wig.gz"
    wigs = [p / name for p in fs.sorted_naturally(clade.glob("chromosome.*"))]
    _log.debug(f"{[str(x) for x in wigs]}")
    outfile = clade / "phastcons.bw"
    is_to_run = not cli.dry_run and fs.is_outdated(outfile, wigs)
    args = ["wigToBigWig", "stdin", chrom_sizes, outfile]
    p = cli.popen_if(is_to_run, args, stdin=cli.PIPE)
    if is_to_run:
        assert p.stdin
        for wig in wigs:
            with gzip.open(wig, "rb") as fin:
                p.stdin.write(fin.read())
                p.stdin.flush()
    p.communicate()
    return outfile


def bigWigInfo(path: Path):
    args = ["bigWigInfo", path]
    return cli.run(args, stdout=cli.PIPE, text=True).stdout


def faToTwoBit(fa_gz: Path):
    outfile = fa_gz.with_suffix("").with_suffix(".2bit")
    if fs.is_outdated(outfile, fa_gz) and not cli.dry_run:
        with gzip.open(fa_gz, "rb") as fin:
            args = ["faToTwoBit", "stdin", outfile]
            p = cli.popen(args, stdin=cli.PIPE)
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
            cli.run(["faSize", "-detailed", genome_fa_gz], stdout=fout)
    return outfile


if __name__ == "__main__":
    main()
