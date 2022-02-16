"""Convert MAF to SAM/BAM/CRAM for visualzation.

src: ./pairwise/{target}/{query}/{chromosome}/sing.maf
dst: ./pairwise/{target}/{query}/cram/genome.cram
"""
import argparse
import concurrent.futures as confu
import logging
import os
import re
from pathlib import Path
from subprocess import PIPE

from . import cli, fs
from .db import ensemblgenomes

_log = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-t", "--test", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("query", nargs="+", type=Path)
    args = parser.parse_args()
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run

    if args.test:
        for path in args.query:
            target_species = "oryza_sativa"
            reference = ensemblgenomes.get_file("*.genome.fa.gz", target_species)
            maf2cram(path, Path(path.parent.name + ".cram"), reference)
        return

    for path in args.query:
        outfile = mafs2cram(path, args.jobs)
        print(outfile)


def mafs2cram(path: Path, jobs: int = 1):
    target_species = path.parent.name
    reference = ensemblgenomes.get_file("*.genome.fa.gz", target_species)
    outdir = path / "cram"
    outdir.mkdir(0o755, exist_ok=True)
    outfile = outdir / "genome.cram"
    crams: list[Path] = []
    with confu.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures: list[confu.Future[Path]] = []
        for dir in list(path.glob("chromosome.*")):
            maf = dir / "sing.maf"
            cram = outdir / (dir.name + ".cram")
            futures.append(executor.submit(maf2cram, maf, cram, reference))
        for future in futures:
            cram = future.result()
            _log.info(cram)
            crams.append(cram)
    cmd = f"samtools merge --no-PG -O CRAM -@ 2 -f -o {str(outfile)} "
    is_to_run = fs.is_outdated(outfile, crams)
    cli.run_if(is_to_run, cmd + " ".join([str(x) for x in crams]))
    cli.run_if(is_to_run, ["samtools", "index", outfile])
    return outfile


def maf2cram(infile: Path, outfile: Path, reference: Path):
    assert infile.exists() and outfile.parent.exists()
    is_to_run = fs.is_outdated(outfile, infile)
    mafconv = cli.popen_if(is_to_run, ["maf-convert", "sam", infile], stdout=PIPE)
    (stdout, _stderr) = mafconv.communicate()
    content = sanitize_cram(is_to_run, reference, stdout)
    cmd = f"samtools sort --no-PG -O CRAM -@ 2 -o {str(outfile)}"
    cli.popen_if(is_to_run, cmd, stdin=PIPE).communicate(content)
    return outfile


def sanitize_cram(cond: bool, reference: Path, sam: bytes):
    def repl(mobj: re.Match[bytes]):
        if int(mobj["flag"]) & 16:  # reverse strand
            qstart = int(mobj["tail_cigar"]) + 1
        else:
            qstart = int(mobj["head_cigar"]) + 1
        qend = qstart + len(mobj["seq"]) - 1
        cells = [
            mobj["qname"] + f":{qstart}-{qend}".encode(),
            mobj["flag"],
            mobj["rname"],
            mobj["pos"],
            mobj["mapq"],
            mobj["cigar"],
            mobj["rnext"],
            mobj["pnext"],
            mobj["tlen"],
            mobj["seq"],
            mobj["misc"],
        ]
        return b"\t".join(cells)

    patt = re.compile(
        rb"^(?P<qname>\S+)\t(?P<flag>\d+)\t"
        rb"\w+\.(?P<rname>\S+)\t(?P<pos>\d+)\t(?P<mapq>\d+)\t"
        rb"(?P<head_cigar>\d+)H(?P<cigar>\S+?)(?P<tail_cigar>\d+)H\t"
        rb"(?P<rnext>\S+)\t(?P<pnext>\w+)\t(?P<tlen>\d+)\t"
        rb"(?P<seq>\w+)\t(?P<misc>.+$)"
    )
    lines = [patt.sub(repl, line) for line in sam.splitlines(keepends=True)]
    cmd = f"samtools view --no-PG -h -C -@ 2 -T {str(reference)}"
    samview = cli.popen_if(cond, cmd, stdin=PIPE, stdout=PIPE)
    (stdout, _stderr) = samview.communicate(b"".join(lines))
    return stdout


if __name__ == "__main__":
    main()
