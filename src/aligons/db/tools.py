"""Preprocessing tools for genome databases."""

import logging
import re
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable

from aligons.extern import htslib, kent
from aligons.util import cli, dl, fs, gff, subp

from . import _rsrc, jgi, mask

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    """CLI for manual execution and testing."""
    parser = cli.ArgumentParser()
    parser.add_argument("url")
    args = parser.parse_args(argv or None)
    dl_mirror_db(args.url)


def fetch_and_bgzip(
    entry: _rsrc.DataSet, prefix: Path
) -> tuple[cli.Future[Path], cli.Future[Path]]:
    """Fetch genome and annotation files, bgzip and index them.

    :param entry: Dataset entry from `_rsrc`.
    :param prefix: Output directory prefix.
    :returns: A tuple of futures:
        (`{stem}.{dna_sm}.genome.fa.gz`, `{stem}.genome.gff3.gz`).
    """
    url_prefix = entry["url_prefix"]
    species = entry["species"]
    annotation = entry["annotation"]
    sequences = entry["sequences"]
    softmasked = re.search(r"dna_sm|softmasked", sequences[0])
    dna_sm = "dna_sm" if softmasked else "dna"
    stem = f"{species}_{ver}" if (ver := entry.get("version", None)) else species
    out_fa = prefix / species / f"{stem}.{dna_sm}.genome.fa.gz"
    out_gff = prefix / species / f"{stem}.genome.gff3.gz"
    responses_fa = [dl_mirror_db(url_prefix + s) for s in sequences]
    ft_fa = cli.thread_submit(_index_compress_concat, responses_fa, out_fa)
    response_gff = dl_mirror_db(url_prefix + annotation)
    ft_gff = cli.thread_submit(index_bgzip, response_gff, out_gff)
    return ft_fa, ft_gff


def dl_mirror_db(url: str) -> dl.Response:
    """Call `dl.mirror` or `jgi.session.mirror` based on URL.

    :param url: URL to download.
    :returns: Download response.
    """
    if "jgi.doe.gov" in url:
        return jgi.session.mirror(url, _rsrc.db_root())
    return dl.mirror(url, _rsrc.db_root())


def _index_compress_concat(responses: Iterable[dl.Response], outfile: Path) -> Path:
    return index_compress_concat([res.path for res in responses], outfile)


def index_compress_concat(
    infiles: list[Path] | list[cli.Future[Path]], outfile: Path
) -> Path:
    """Concatenate FASTA or GFF files into a bgzipped file and index if applicable.

    :param infiles: A list of FASTA, GFF files, or futures of them.
    :param outfile: The output path of the concatenated bgzip.
    :returns: The same path as `outfile`.
    """
    infiles = [cli.result(f) for f in infiles]
    infiles = fs.sorted_naturally(infiles)
    if len(infiles) == 1:
        if not cli.dry_run:
            outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
            fs.symlink(infiles[0], outfile)
        try:
            htslib.try_index(outfile)
        except subprocess.CalledProcessError as e:
            _log.debug(f"indexing failed for {outfile}: {e}")
            outfile.unlink()
        else:
            return outfile
    htslib.concat_bgzip(infiles, outfile)
    htslib.try_index(outfile)
    return outfile


def index_bgzip(infile: Path | dl.Response, outfile: Path) -> Path:
    """Bgzip FASTA or GFF file and index if applicable.

    :param infile: A FASTA or GFF file, possibly a download response.
    :param outfile: The output bgzipped file.
    :returns: The same path as `outfile`, not index file.
    """
    if isinstance(infile, dl.Response):
        infile = infile.path
    if fs.is_outdated(outfile, infile):
        if not cli.dry_run:
            outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
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


def unzip_index_bgzip(infile_zip: Path, name: str, outfile: Path) -> Path:
    """Make bgzip and index of a specific file inside a zip archive.

    :param infile_zip: Input zip file like NCBI genome package.
    :param name: File inside the zip to be extracted.
    :param outfile: Output bgzipped file.
    :returns: The same path as `outfile`, not index file.
    """
    if fs.is_outdated(outfile, infile_zip):
        if not cli.dry_run:
            outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        if any(s in (".gff", ".gff3") for s in outfile.suffixes):
            with subp.popen_zcat([infile_zip, Path(name)]) as zcat:
                assert zcat.stdout
                gff3 = gff.GFF(zcat.stdout)
            with htslib.popen_bgzip(outfile) as bgzip:
                assert bgzip.stdin
                gff3.sanitize().write(bgzip.stdin)
        else:
            with subp.popen_zcat([infile_zip, Path(name)]) as zcat:
                htslib.bgzip(zcat.stdout, outfile)
    htslib.try_index(outfile)
    return outfile


def softmask(genome: Path, species: str = "") -> cli.Future[Path]:
    """Softmask a genome FASTA.

    :param genome: A bgzipped genome FASTA file.
    :param species: A taxonomic group name for RepeatMasker/Dfam.
        If empty, the parent directory name of the genome file is used.
    :returns: A future of the softmasked genome FASTA with ".dna_sm."
        in the file name.
    """
    if ".dna_sm." in genome.name:
        fs.expect_suffix(genome, ".gz")
        return cli.thread_submit(cli.result, genome)
    (outname, count) = re.subn(r"\.dna\.", ".dna_sm.", genome.name)
    assert count == 1, genome
    species = species or genome.parent.name
    masked = genome.parent / outname
    ft_chromosomes = _split_genome_fa(genome, "_work")
    fts = [mask.submit(ft, species) for ft in ft_chromosomes]
    return cli.thread_submit(index_compress_concat, fts, masked)


def _split_genome_fa(genome: Path, subdir: str) -> list[cli.Future[Path]]:
    fai = htslib.read_fai(genome)
    _log.debug(f"{fai = }")
    min_size = 1000000
    workdir = genome.parent / subdir
    workdir.mkdir(0o755, exist_ok=True)
    fts: list[cli.Future[Path]] = []
    for seqid, size in fai.items():
        if re.search(r"scaffold|contig|Egra_v1_0", seqid):
            _log.info(f"ignoring {seqid} in {genome}")
            continue
        if size < min_size:
            _log.info(f"{genome}:{seqid} {size} < {min_size}")
            continue
        name = genome.name.replace(".genome.fa.gz", f".chromosome.{seqid}.fa", 1)
        chr_fa = workdir / name
        fts.append(cli.thread_submit(htslib.faidx_query, genome, seqid, chr_fa))
    return fts


def genome_to_twobits(genome: Path | cli.Future[Path]) -> list[cli.Future[Path]]:
    """Split a genome FASTA into chromosomes in 2bit format.

    :param genome: A bgzipped genome FASTA file.
    :returns: A list of futures of 2bit files for each chromosome.
    """
    genome = cli.result(genome)
    fai = htslib.read_fai(genome)
    _log.debug(f"{fai = }")
    min_size = 1000000
    fts: list[cli.Future[Path]] = []
    for seqid, size in fai.items():
        if re.search(r"scaffold|contig", seqid):
            _log.info(f"ignoring {seqid} in {genome}")
            continue
        if size < min_size:
            _log.info(f"{genome}:{seqid} {size} < {min_size}")
            continue
        name = genome.name.replace(".genome.fa.gz", f".chromosome.{seqid}.2bit", 1)
        chr2bit = genome.with_name(name)
        fts.append(cli.thread_submit(faidx_twobit, genome, seqid, chr2bit))
    return fts


def faidx_twobit(fa_gz: Path, seqid: str, twobit: Path) -> Path:
    """Make 2bit file of a specific sequence from a bgzipped FASTA.

    :param fa_gz: A bgzipped FASTA file.
    :param seqid: A sequence name contained in `fa_gz`.
    :param twobit: The output 2bit file.
    :returns: The same path as `twobit`.
    """
    if_ = fs.is_outdated(twobit, fa_gz)
    with htslib.popen_faidx_query(fa_gz, seqid, if_=if_) as p:
        return kent.faToTwoBit(p.stdout, twobit, if_=if_)


def bgzip_or_symlink(infile: Path | dl.Response, outfile: Path) -> Path:
    """Decompress/compress/symlink depending on file names."""
    if isinstance(infile, dl.Response):
        infile = infile.path
    if fs.is_outdated(outfile, infile) and not cli.dry_run:
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        if htslib.to_be_bgzipped(outfile.name):
            index_bgzip(infile, outfile)
        elif outfile.suffix == infile.suffix or (
            outfile.suffix == ".bw" and infile.suffix.lower() == ".bigwig"
        ):
            fs.symlink(infile, outfile, relative=True)
        elif outfile.suffix == ".gz":
            _log.debug(f"rare case: gzip {infile = } to {outfile = }")
            with subp.popen_zcat(infile) as zcat:
                subp.gzip(zcat.stdout, outfile)
        else:
            msg = f"unexpected formats: {infile = }, {outfile = }"
            raise ValueError(msg)
    return fs.print_if_exists(outfile)


if __name__ == "__main__":
    main()
