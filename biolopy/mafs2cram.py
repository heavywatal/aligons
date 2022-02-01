"""Convert MAF to SAM/BAM/CRAM for visualzation.

src: ./pairwise/{target}/{query}/{chromosome}/sing.maf
dst: ./pairwise/{target}/{query}/cram/genome.cram
"""
import concurrent.futures as confu
import logging
import os
import re
import subprocess
from subprocess import PIPE
from pathlib import Path
from typing import IO

from . import cli, fs
from .db import name
from .db import ensemblgenomes

_log = logging.getLogger(__name__)
_dry_run = False


def main():
    import argparse

    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("query", nargs="+", type=Path)
    args = parser.parse_args()
    cli.logging_config(args.loglevel)
    global _dry_run
    _dry_run = args.dry_run

    for path in args.query:
        outfile = mafs2cram(path, args.jobs)
        print(outfile)


def mafs2cram(path: Path, jobs: int = 1):
    target_species = path.parent.name
    reference = ensemblgenomes.get_file("*.genome.fa.gz", target_species)
    outdir = path / "cram"
    outdir.mkdir(0o755, exist_ok=True)
    outfile = outdir / "genome.cram"
    crams: list[str] = []
    with confu.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures: list[confu.Future[Path]] = []
        for dir in list(path.glob("chromosome.*")):
            maf = dir / "sing.maf"
            cram = outdir / (dir.name + ".cram")
            futures.append(executor.submit(maf2cram, maf, cram, reference))
        for future in futures:
            cram = str(future.result())
            _log.info(cram)
            crams.append(cram)
    cmd = f"samtools merge --no-PG -O CRAM -@ 2 -f -o {str(outfile)} "
    outfile_is_outdated = fs.is_outdated(outfile, Path(crams[0]))
    popen_if(outfile_is_outdated, cmd + " ".join(crams)).communicate()
    popen_if(outfile_is_outdated, f"samtools index {str(outfile)}").communicate()
    return outfile


def maf2cram(infile: Path, outfile: Path, reference: Path):
    assert infile.exists() and outfile.parent.exists()
    cond = fs.is_outdated(outfile)
    mafconv = popen_if(cond, f"maf-convert sam {str(infile)}", stdout=PIPE)
    samview = sanitize_cram(cond, reference, mafconv.stdout)
    cmd = f"samtools sort --no-PG -O CRAM -@ 2 -o {str(outfile)}"
    popen_if(cond, cmd, stdin=samview.stdout).communicate()
    return outfile


def sanitize_cram(cond: bool, reference: Path, stdin: IO[bytes] | None):
    shortname = name.shorten(reference.parent.parent.name)
    patt_refseq = re.compile(rf"(?<=\t){shortname}\.".encode())
    patt_cigar = re.compile(rb"\d+H")
    cmd = f"samtools view --no-PG -h -C -@ 2 -T {str(reference)}"
    samview = popen_if(cond, cmd, stdin=PIPE, stdout=PIPE)
    assert samview.stdin
    for line in stdin or []:
        line = patt_refseq.sub(b"", line, 1)
        line = patt_cigar.sub(b"", line, 2)
        samview.stdin.write(line)
    samview.stdin.close()
    return samview


def sanitize_cram_sed(cond: bool, reference: Path, stdin: IO[bytes] | None):
    shortname = name.shorten(reference.parent.parent.name)
    sed_refseq = popen_if(cond, f"sed -e s/{shortname}.//", stdin=stdin, stdout=PIPE)
    sed_cigar = popen_if(
        cond, r"sed -e s/[0-9]\+H//g", stdin=sed_refseq.stdout, stdout=PIPE
    )
    cmd = f"samtools view --no-PG -h -C -@ 2 -T {str(reference)}"
    return popen_if(cond, cmd, stdin=sed_cigar.stdout, stdout=PIPE)


def popen_if(
    cond: bool,
    args: list[str] | str,
    stdin: IO[bytes] | int | None = None,
    stdout: IO[bytes] | int | None = None,
):
    (args, cmd) = cli.prepare_args(args, (not cond) or _dry_run)
    _log.info(cmd)
    return subprocess.Popen(args, stdin=stdin, stdout=stdout)


if __name__ == "__main__":
    main()
