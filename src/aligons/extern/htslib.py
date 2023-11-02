import concurrent.futures as confu
import gzip
import logging
import re
from collections.abc import Iterable
from pathlib import Path
from typing import IO

from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-O", "--outdir", type=Path)
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    fts = split_fa_gz(args.infile, outdir=args.outdir)
    cli.wait_raise(fts)


def split_fa_gz(
    bgz: Path,
    fmt: str = "{stem}.part_{seqid}.fa.gz",
    sub: tuple[str, str] = ("", ""),
    outdir: Path | None = None,
) -> list[cli.FuturePath]:
    if outdir is None:
        outdir = bgz.parent
    elif not cli.dry_run:
        outdir.mkdir(0o755, parents=True, exist_ok=True)
    stem = bgz.with_suffix("").stem
    if sub[0]:
        stem = re.sub(sub[0], sub[1], stem)
    fai = faidx(bgz)
    if cli.dry_run and not fai.exists():
        return []
    with fai.open("rt") as fin:
        seqids = [line.split()[0] for line in fin]
    fts: list[cli.FuturePath] = []
    for seqid in seqids:
        if re.search(r"scaffold|contig", seqid):
            _log.debug("ignoring {seqid} in {bgz}")
            continue
        outfile = outdir / fmt.format(stem=stem, seqid=seqid)
        fts.append(cli.thread_submit(faidx_query, bgz, seqid, outfile))
    return fts


def faidx_query(bgz: Path, region: str, outfile: Path) -> Path:
    args: subp.Args = ["samtools", "faidx", bgz, region]
    is_to_run = fs.is_outdated(outfile, bgz)
    if outfile.suffix == ".gz":
        with subp.popen(args, stdout=subp.PIPE, if_=is_to_run) as p:
            if is_to_run:
                bgzip(p.stdout, outfile)
    else:
        args.extend(["-o", outfile])
        subp.run(args, if_=is_to_run)
    _log.info(f"{outfile}")
    return outfile


def concat_bgzip(infiles: list[Path], outfile: Path) -> Path:
    if fs.is_outdated(outfile, infiles) and not cli.dry_run:
        with outfile.open("wb") as fout:
            bgzip = subp.popen("bgzip -@2", stdin=subp.PIPE, stdout=fout)
            assert bgzip.stdin is not None
            if ".gff" in outfile.name:
                header = collect_gff3_header(infiles)
                bgzip.stdin.write(header)
                bgzip.stdin.flush()
                _log.debug(header.decode())
            for file in infiles:
                if ".gff" in outfile.name:
                    sort_clean_chromosome_gff3(file, bgzip.stdin)
                else:
                    subp.popen_zcat(file, bgzip.stdin).communicate()
            bgzip.communicate()
    _log.info(f"{outfile}")
    return outfile


def collect_gff3_header(infiles: Iterable[Path]) -> bytes:
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


def bgzip(data: bytes | IO[bytes] | None, outfile: Path) -> Path:
    """https://www.htslib.org/doc/bgzip.html."""
    assert outfile.suffix == ".gz", outfile
    if outfile.exists():
        _log.info("overwriting {outfile}")
    if data and not cli.dry_run:
        with outfile.open("wb") as fout:
            if isinstance(data, bytes):
                subp.run(["bgzip", "-@2"], input=data, stdout=fout)
            else:
                subp.run(["bgzip", "-@2"], stdin=data, stdout=fout)
    return outfile


def try_index(bgz: Path | cli.FuturePath) -> Path:
    if isinstance(bgz, confu.Future):
        bgz = bgz.result()
    if to_be_tabixed(bgz.name):
        return tabix(bgz)
    if to_be_faidxed(bgz.name):
        return faidx(bgz)
    return bgz


def faidx(bgz: Path | cli.FuturePath) -> Path:
    """https://www.htslib.org/doc/samtools-faidx.html."""
    if isinstance(bgz, confu.Future):
        bgz = bgz.result()
    outfile = bgz.with_suffix(bgz.suffix + ".fai")
    subp.run(["samtools", "faidx", bgz], if_=fs.is_outdated(outfile, bgz))
    _log.info(f"{outfile}")
    return outfile


def tabix(bgz: Path | cli.FuturePath) -> Path:
    """https://www.htslib.org/doc/tabix.html.

    Use .csi instead of .tbi for chromosomes >512 Mbp e.g., atau, hvul.
    """
    if isinstance(bgz, confu.Future):
        bgz = bgz.result()
    outfile = bgz.with_suffix(bgz.suffix + ".csi")
    subp.run(["tabix", "--csi", bgz], if_=fs.is_outdated(outfile, bgz))
    _log.info(f"{outfile}")
    return outfile


def index(cram: Path | cli.FuturePath) -> Path:
    """https://www.htslib.org/doc/samtools-index.html."""
    if isinstance(cram, confu.Future):
        cram = cram.result()
    outfile = cram.with_suffix(cram.suffix + ".crai")
    subp.run(["samtools", "index", cram], if_=fs.is_outdated(outfile, cram))
    _log.info(f"{outfile}")
    return outfile


def to_be_bgzipped(filename: str) -> bool:
    return to_be_faidxed(filename) or to_be_tabixed(filename)


def to_be_faidxed(filename: str) -> bool:
    ext = (".fa", ".fas", ".fasta", ".fna")
    return filename.removesuffix(".gz").removesuffix(".zip").endswith(ext)


def to_be_tabixed(filename: str) -> bool:
    ext = (".gff", ".gff3", ".gtf", ".bed")
    return filename.removesuffix(".gz").removesuffix(".zip").endswith(ext)


def sort_clean_chromosome_gff3(infile: Path, stdout: IO[bytes]) -> None:
    # TODO: jbrowse2 still needs billzt/gff3sort precision?
    p1 = subp.popen(f"zgrep -v '^#' {infile!s}", stdout=subp.PIPE, quiet=True)
    p2 = subp.popen(
        "grep -v '\tchromosome\t'", stdin=p1.stdout, stdout=subp.PIPE, quiet=True
    )
    if p1.stdout:
        p1.stdout.close()
    p3 = subp.popen("sort -k4,4n", stdin=p2.stdout, stdout=stdout, quiet=True)
    if p2.stdout:
        p2.stdout.close()
    p3.communicate()


if __name__ == "__main__":
    main()
