"""Convert MAF to SAM/BAM/CRAM for visualzation.

src: ./pairwise/{target}/{query}/{chromosome}/sing.maf
dst: ./pairwise/{target}/{query}/cram/genome.cram
"""
import concurrent.futures as confu
import logging
import re
from pathlib import Path

from aligons.db import ensemblgenomes, phylo
from aligons.util import cli, config, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-t", "--test", action="store_true")
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
    cli.wait_raise([mafs2cram(path) for path in args.query])


def run(target: Path, clade: str):
    tree = phylo.newicks[clade]
    tips = phylo.extract_names(tree)
    query_names = ensemblgenomes.sanitize_queries(target.name, tips)
    return [mafs2cram(target / q) for q in query_names]


def mafs2cram(path: Path):
    target_species = path.parent.name
    subdir = "kmer" if config["db"]["kmer"] else ""
    reference = ensemblgenomes.get_file("*.genome.fa.gz", target_species, subdir)
    outdir = path / "cram"
    if not cli.dry_run:
        outdir.mkdir(0o755, exist_ok=True)
    pool = cli.ThreadPool()
    futures: list[confu.Future[Path]] = []
    for chr_dir in fs.sorted_naturally(path.glob("chromosome.*")):
        maf = chr_dir / "sing.maf"
        if not maf.exists():
            _log.warning(f"not found {maf}")
            continue
        cram = outdir / (chr_dir.name + ".cram")
        futures.append(pool.submit(maf2cram, maf, cram, reference))
    return pool.submit(merge_crams, futures)


def merge_crams(futures: list[confu.Future[Path]]):
    crams: list[Path] = [f.result() for f in futures]
    outdir = crams[0].parent
    outfile = outdir / "genome.cram"
    is_to_run = bool(crams) and fs.is_outdated(outfile, crams)
    cmd = f"samtools merge --no-PG -O CRAM -@ 2 -f -o {outfile!s} "
    cmd += " ".join([str(x) for x in crams])
    subp.run_if(is_to_run, cmd)
    subp.run_if(is_to_run, ["samtools", "index", outfile])
    if outfile.exists():
        print(outfile)
    return outfile


def maf2cram(infile: Path, outfile: Path, reference: Path):
    is_to_run = fs.is_outdated(outfile, infile)
    mafconv = subp.popen_if(is_to_run, ["maf-convert", "sam", infile], stdout=subp.PIPE)
    (stdout, _stderr) = mafconv.communicate()
    content = sanitize_cram(is_to_run, reference, stdout)
    cmd = f"samtools sort --no-PG -O CRAM -@ 2 -o {outfile!s}"
    subp.popen_if(is_to_run, cmd, stdin=subp.PIPE).communicate(content)
    _log.info(f"{outfile}")
    return outfile


def sanitize_cram(cond: bool, reference: Path, sam: bytes):  # noqa: FBT001
    def repl(mobj: re.Match[bytes]):
        qstart = 0
        if int(mobj["flag"]) & 16:  # reverse strand
            if tail := mobj["tail_cigar"]:
                qstart = int(tail.rstrip(b"H")) + 1
        elif head := mobj["head_cigar"]:
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
    cmd = f"samtools view --no-PG -h -C -@ 2 -T {reference!s}"
    samview = subp.popen_if(cond, cmd, stdin=subp.PIPE, stdout=subp.PIPE)
    (stdout, _stderr) = samview.communicate(b"".join(lines))
    return stdout


if __name__ == "__main__":
    main()
