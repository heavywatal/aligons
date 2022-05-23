import gzip
import logging
import re
from collections.abc import Iterable
from pathlib import Path

from ..util import cli, fs, subp

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


def bgzip(infiles: list[Path], outfile: Path):
    if fs.is_outdated(outfile, infiles) and not cli.dry_run:
        with open(outfile, "wb") as fout:
            bgzip = subp.popen("bgzip -@2", stdin=subp.PIPE, stdout=fout)
            assert bgzip.stdin
            if ".gff" in outfile.name:
                header = collect_gff3_header(infiles)
                bgzip.stdin.write(header)
                bgzip.stdin.flush()
                _log.debug(header.decode())
            for file in infiles:
                if ".gff" in outfile.name:
                    p = sort_clean_chromosome_gff3(file)
                    (stdout, _stderr) = p.communicate()
                    bgzip.stdin.write(stdout)
                    bgzip.stdin.flush()
                else:
                    with gzip.open(file, "rb") as fin:
                        bgzip.stdin.write(fin.read())
                        bgzip.stdin.flush()
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
        subp.run(["samtools", "faidx", bgz])
    return outfile


def tabix(bgz: Path):
    """http://www.htslib.org/doc/tabix.html
    Use .csi instead of .tbi for chromosomes >512 Mbp e.g., atau, hvul
    """
    outfile = bgz.with_suffix(bgz.suffix + ".csi")
    if fs.is_outdated(outfile, bgz):
        subp.run(["tabix", "--csi", bgz])
    return outfile


def sort_clean_chromosome_gff3(infile: Path):
    # TODO: jbrowse2 still needs billzt/gff3sort precision?
    p1 = subp.popen(f"zgrep -v '^#' {str(infile)}", stdout=subp.PIPE, quiet=True)
    p2 = subp.popen(
        "grep -v '\tchromosome\t'", stdin=p1.stdout, stdout=subp.PIPE, quiet=True
    )
    if p1.stdout:
        p1.stdout.close()
    p3 = subp.popen("sort -k4,4n", stdin=p2.stdout, stdout=subp.PIPE, quiet=True)
    if p2.stdout:
        p2.stdout.close()
    return p3
