"""Convert MAF to SAM/BAM/CRAM for visualzation.

src: ./pairwise/{target}/{query}/{chromosome}/sing.maf
dst: ./pairwise/{target}/{query}/cram/genome.cram
"""
import concurrent.futures as confu
import logging
import os
import re
from pathlib import Path
from typing import Iterable

from ..db import ensemblgenomes, phylo
from ..util import cli, config, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = cli.ArgumentParser()
    parser.add_argument("-t", "--test", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("query", nargs="*", type=Path)  # pairwise/{target}/{query}
    args = parser.parse_args(argv or None)
    subdir = "kmer" if config["db"]["kmer"] else ""
    if args.test:
        for path in args.query:
            target = path.parent.parent.parent.name
            reference = ensemblgenomes.get_file("*.genome.fa.gz", target, subdir)
            stem = str(path.parent).replace("/", "_")
            maf2cram(path, Path(stem + ".cram"), reference)
        return
    _run(args.query, args.jobs)


def run(target: Path, clade: str, jobs: int):
    tree = phylo.newicks[clade]
    tips = phylo.extract_names(tree)
    query_names = ensemblgenomes.sanitize_queries(target.name, tips)
    _run([target / q for q in query_names], jobs)


def _run(queries: Iterable[Path], jobs: int):
    for path in queries:
        if (cram := mafs2cram(path, jobs)).exists():
            print(cram)


def mafs2cram(path: Path, jobs: int = 1):
    target_species = path.parent.name
    subdir = "kmer" if config["db"]["kmer"] else ""
    reference = ensemblgenomes.get_file("*.genome.fa.gz", target_species, subdir)
    outdir = path / "cram"
    if not cli.dry_run:
        outdir.mkdir(0o755, exist_ok=True)
    crams: list[Path] = []
    with confu.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures: list[confu.Future[Path]] = []
        for dir in fs.sorted_naturally(path.glob("chromosome.*")):
            maf = dir / "sing.maf"
            if not maf.exists():
                _log.warning(f"not found {maf}")
                continue
            cram = outdir / (dir.name + ".cram")
            futures.append(executor.submit(maf2cram, maf, cram, reference))
        for future in futures:
            cram = future.result()
            _log.info(f"{cram}")
            crams.append(cram)
    outfile = outdir / "genome.cram"
    is_to_run = bool(crams) and fs.is_outdated(outfile, crams)
    cmd = f"samtools merge --no-PG -O CRAM -@ 2 -f -o {str(outfile)} "
    cmd += " ".join([str(x) for x in crams])
    subp.run_if(is_to_run, cmd)
    subp.run_if(is_to_run, ["samtools", "index", outfile])
    return outfile


def maf2cram(infile: Path, outfile: Path, reference: Path):
    is_to_run = fs.is_outdated(outfile, infile)
    mafconv = subp.popen_if(is_to_run, ["maf-convert", "sam", infile], stdout=subp.PIPE)
    (stdout, _stderr) = mafconv.communicate()
    content = sanitize_cram(is_to_run, reference, stdout)
    cmd = f"samtools sort --no-PG -O CRAM -@ 2 -o {str(outfile)}"
    subp.popen_if(is_to_run, cmd, stdin=subp.PIPE).communicate(content)
    return outfile


def sanitize_cram(cond: bool, reference: Path, sam: bytes):
    def repl(mobj: re.Match[bytes]):
        qstart = 0
        if int(mobj["flag"]) & 16:  # reverse strand
            if tail := mobj["tail_cigar"]:
                qstart = int(tail.rstrip(b"H")) + 1
        else:
            if head := mobj["head_cigar"]:
                qstart = int(head.rstrip(b"H")) + 1
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
        rb"(?P<head_cigar>\d+H)?(?P<cigar>\S+?)(?P<tail_cigar>\d+H)?\t"
        rb"(?P<rnext>\S+)\t(?P<pnext>\w+)\t(?P<tlen>\d+)\t"
        rb"(?P<seq>\w+)\t(?P<misc>.+$)"
    )
    lines = [patt.sub(repl, line) for line in sam.splitlines(keepends=True)]
    cmd = f"samtools view --no-PG -h -C -@ 2 -T {str(reference)}"
    samview = subp.popen_if(cond, cmd, stdin=subp.PIPE, stdout=subp.PIPE)
    (stdout, _stderr) = samview.communicate(b"".join(lines))
    return stdout


if __name__ == "__main__":
    main()
