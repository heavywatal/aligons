import logging
import re
from collections.abc import Iterable
from pathlib import Path

from aligons.extern import htslib, kent
from aligons.util import cli, dl, fs, gff, subp

from . import _rsrc, jgi, mask

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("url")
    args = parser.parse_args(argv or None)
    dl_mirror_db(args.url)


def fetch_and_bgzip(
    entry: _rsrc.DataSet, prefix: Path
) -> tuple[cli.Future[Path], cli.Future[Path]]:
    url_prefix = entry["url_prefix"]
    species = entry["species"]
    annotation = entry["annotation"]
    sequences = entry["sequences"]
    softmasked = re.search(r"dna_sm|softmasked", sequences[0])
    dna_sm = "dna_sm" if softmasked else "dna"
    stem = f"{species}_{ver}" if (ver := entry.get("version", None)) else species
    out_fa = prefix / species / f"{stem}.{dna_sm}.toplevel.fa.gz"
    out_gff = prefix / species / f"{stem}.genome.gff3.gz"
    responses_fa = [dl_mirror_db(url_prefix + s) for s in sequences]
    ft_fa = cli.thread_submit(_index_compress_concat, responses_fa, out_fa)
    response_gff = dl_mirror_db(url_prefix + annotation)
    ft_gff = cli.thread_submit(index_bgzip, response_gff, out_gff)
    return ft_fa, ft_gff


def dl_mirror_db(url: str) -> dl.Response:
    if "jgi.doe.gov" in url:
        return jgi.session.mirror(url, _rsrc.db_root())
    return dl.mirror(url, _rsrc.db_root())


def process_genome(
    futures: tuple[cli.Future[Path], cli.Future[Path]],
) -> list[cli.Future[Path]]:
    ft_fa, ft_gff = futures
    ft_masked = mask.submit(ft_fa)
    return [ft_gff, cli.thread_submit(from_genome, ft_masked)]


def _index_compress_concat(responses: Iterable[dl.Response], outfile: Path) -> Path:
    return index_compress_concat([res.path for res in responses], outfile)


def index_compress_concat(infiles: list[Path], outfile: Path) -> Path:
    infiles = fs.sorted_naturally(infiles)
    htslib.concat_bgzip(infiles, outfile)
    htslib.try_index(outfile)
    return outfile


def index_bgzip(infile: Path | dl.Response, outfile: Path) -> Path:
    if isinstance(infile, dl.Response):
        infile = infile.path
    if fs.is_outdated(outfile, infile):
        if any(s in (".gff", ".gff3") for s in outfile.suffixes):
            gff3 = gff.GFF(infile)
            with htslib.popen_bgzip(outfile) as bgzip:
                assert bgzip.stdin
                gff3.sanitize().write(bgzip.stdin)
        else:
            with subp.popen_zcat(infile) as zcat:
                htslib.bgzip(zcat.stdout, outfile)
    htslib.try_index(outfile)
    return outfile


def from_genome(genome: Path | cli.Future[Path]) -> list[Path]:
    genome = cli.result(genome)
    fasize = read_fasize(genome)
    min_size = 1000000
    chr_2bits: list[Path] = []
    for seqid, size in fasize.items():
        if re.search(r"scaffold|contig", seqid):
            _log.info(f"ignoring {seqid} in {genome}")
            continue
        if size < min_size:
            _log.warning(f"{genome}:{seqid} {size} < {min_size}")
            continue
        name = genome.name.replace(".genome.fa.gz", f".chromosome.{seqid}.2bit", 1)
        chr2bit = genome.with_name(name)
        if fs.is_outdated(chr2bit, genome):
            chr_2bits.append(faidx_twobit(genome, seqid, chr2bit))
    return chr_2bits


def faidx_twobit(fa_gz: Path, seqid: str, twobit: Path) -> Path:
    with htslib.faidx_query(fa_gz, seqid) as p:
        return kent.faToTwoBit(p.stdout, twobit)


def read_fasize(genome: Path) -> dict[str, int]:
    fasize = kent.faSize(genome)
    res: dict[str, int] = {}
    with fasize.open("rt") as fin:
        for line in fin:
            (seqid, size) = line.split(maxsplit=1)
            res[seqid] = int(size)
    return res


def bgzip_or_symlink(infile: Path | dl.Response, outfile: Path) -> Path:
    """Decompress/compress/symlink depending on file names."""
    if isinstance(infile, dl.Response):
        infile = infile.path
    if fs.is_outdated(outfile, infile) and not cli.dry_run:
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        if htslib.to_be_bgzipped(outfile.name):
            index_bgzip(infile, outfile)
        elif outfile.suffix == infile.suffix:
            fs.symlink(infile, outfile, relative=True)
        elif outfile.suffix == ".gz":
            _log.debug(f"rare case: gzip {infile = } to {outfile = }")
            with subp.popen_zcat(infile) as zcat:
                subp.gzip(zcat.stdout, outfile)
        else:
            msg = f"unexpected formats: {infile = }, {outfile = }"
            raise ValueError(msg)
    _log.info(f"{outfile}")
    return outfile


if __name__ == "__main__":
    main()
