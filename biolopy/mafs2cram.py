"""
Convert MAF to SAM/BAM/CRAM for visualzation.
"""
import concurrent.futures as confu
import logging
import os
import re
import subprocess
from pathlib import Path
from typing import IO

from . import cli
from .db import name
from .db import ensemblgenomes

_log = logging.getLogger(__name__)


def mafs2cram(path: Path, jobs: int = 1):
    (target, _query) = path.name.split("_")
    target_species = ensemblgenomes.expand_shortname(target)
    reference = ensemblgenomes.get_file("*.genome.fa.gz", target_species)
    crams: list[str] = []
    with confu.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures: list[confu.Future[Path]] = []
        for dir in list(path.glob("chromosome.*")):
            futures.append(executor.submit(maf2cram, dir, reference))
        for future in futures:
            sing_cram = str(future.result())
            _log.info(sing_cram)
            crams.append(sing_cram)
    outfile = path / f"pairwise-{path.name}.cram"
    samtools_merge(["-f", "-o", str(outfile)] + crams)
    return outfile


def maf2cram(chromosome_dir: Path, reference: Path):
    infile = chromosome_dir / "sing.maf"
    outfile = chromosome_dir / "sing.cram"
    assert infile.exists()
    mafconv = maf_convert(infile)
    samview = sanitize_cram(reference, mafconv.stdout)
    samtools_sort(["-o", str(outfile)], stdin=samview.stdout)
    return outfile


def sanitize_cram(reference: Path, stdin: IO[str] | None):
    shortname = name.shorten(reference.parent.parent.name)
    patt_refseq = re.compile(fr"(?<=\t){shortname}\.")
    patt_cigar = re.compile(r"\d+H")
    viewargs = ["--reference", str(reference)]
    samview = samtools_view(viewargs, stdin=subprocess.PIPE)
    assert samview.stdin
    for line in stdin or []:
        line = patt_refseq.sub("", line, 1)
        line = patt_cigar.sub("", line, 2)
        samview.stdin.write(line)
    samview.stdin.close()
    return samview


def sanitize_cram_sed(reference: Path, stdin: IO[str] | None):
    shortname = name.shorten(reference.parent.parent.name)
    sed_refseq = sed(f"s/{shortname}.//", stdin=stdin)
    sed_cigar = sed(r"s/[0-9]\+H//g", stdin=sed_refseq.stdout)
    viewargs = ["--reference", str(reference)]
    return samtools_view(viewargs, stdin=sed_cigar.stdout)


def maf_convert(maf: Path, format: str = "sam"):
    """[LAST](https://gitlab.com/mcfrith/last)"""
    args = ["maf-convert", format, str(maf)]
    _log.info(" ".join(args))
    return subprocess.Popen(args, stdout=subprocess.PIPE, text=True)


def samtools_view(args: list[str], stdin: IO[str] | int | None):
    args = ["samtools", "view", "--no-PG", "-h", "-C", "-@", "2"] + args
    _log.info(" ".join(args))
    return subprocess.Popen(args, stdin=stdin, stdout=subprocess.PIPE, text=True)


def samtools_sort(args: list[str], stdin: IO[str] | int | None):
    args = ["samtools", "sort", "--no-PG", "-O", "CRAM", "-@", "2"] + args
    _log.info(" ".join(args))
    return subprocess.run(args, stdin=stdin, text=True)


def samtools_merge(args: list[str]):
    args = ["samtools", "merge", "--no-PG", "-O", "CRAM", "-@", "2"] + args
    _log.info(" ".join(args))
    return subprocess.run(args, text=True)


def sed(expr: str, stdin: IO[str] | int | None):
    args = ["sed", "-e", expr]
    _log.info(" ".join(args))
    return subprocess.Popen(args, stdin=stdin, stdout=subprocess.PIPE, text=True)


def main():
    import argparse

    parser = argparse.ArgumentParser(parents=[cli.logging_argparser("v")])
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("paths", nargs="+", type=Path)
    args = parser.parse_args()
    cli.logging_config(args.loglevel)

    for path in args.paths:
        outfile = mafs2cram(path, args.jobs)
        print(outfile)


if __name__ == "__main__":
    main()
