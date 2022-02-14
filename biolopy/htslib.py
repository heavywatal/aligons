import gzip
import itertools
import logging
import re
from collections.abc import Iterable
from pathlib import Path
from subprocess import PIPE

from . import cli, fs

_log = logging.getLogger(__name__)


def create_genome_bgzip(path: Path):
    """Combine chromosome files and bgzip it."""
    if (ext := path.parent.name) != "gff3":
        ext = "fa"
    files = fs.sorted_naturally(path.glob(rf"*.chromosome.*.{ext}.gz"))
    assert files
    _log.debug(str(files))
    name = files[0].name
    (outname, count) = re.subn(rf"\.chromosome\..+\.{ext}", rf".genome.{ext}", name)
    assert count == 1
    outfile = path / outname
    return bgzip(files, outfile)


def bgzip(infiles: Iterable[Path], outfile: Path):
    if fs.is_outdated(outfile) and not cli.dry_run:
        with open(outfile, "wb") as fout:
            bgzip = cli.popen("bgzip -@2", stdin=PIPE, stdout=fout)
            assert bgzip.stdin
            if ".gff" in outfile.name:
                infiles, it2 = itertools.tee(infiles)
                header = collect_gff3_header(it2)
                bgzip.stdin.write(header)
                _log.debug(header.decode())
            for file in infiles:
                if ".gff" in outfile.name:
                    p = sort_clean_chromosome_gff3(file)
                    assert p.stdout
                    bgzip.stdin.write(p.stdout.read())
                else:
                    with gzip.open(file, "rb") as fin:
                        bgzip.stdin.write(fin.read())
            bgzip.communicate()
    return outfile


def collect_gff3_header(infiles: Iterable[Path]):
    header = b"##gff-version 3\n"
    for file in infiles:
        with gzip.open(file, "rt") as fin:
            for line in fin:
                if line.startswith("##sequence-region"):
                    header += line.encode()
                    break
                if not line.startswith("#"):
                    break
    return header


def faidx(bgz: Path):
    """http://www.htslib.org/doc/samtools-faidx.html"""
    outfile = bgz.with_suffix(bgz.suffix + ".fai")
    if fs.is_outdated(outfile, bgz):
        cli.run(["samtools", "faidx", str(bgz)])
    return outfile


def tabix(bgz: Path):
    """http://www.htslib.org/doc/tabix.html"""
    outfile = bgz.with_suffix(bgz.suffix + ".tbi")
    if fs.is_outdated(outfile, bgz):
        cli.run(["tabix", str(bgz)])
    return outfile


def sort_clean_chromosome_gff3(infile: Path):
    # TODO: jbrowse2 still needs billzt/gff3sort precision?
    p1 = cli.popen(f"zgrep -v '^#' {str(infile)}", stdout=PIPE)
    p2 = cli.popen("grep -v '\tchromosome\t'", stdin=p1.stdout, stdout=PIPE)
    return cli.popen("sort -k4,4n", stdin=p2.stdout, stdout=PIPE)
