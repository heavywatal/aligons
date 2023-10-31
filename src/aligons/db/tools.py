import concurrent.futures as confu
import logging
import os
import re
from pathlib import Path

from aligons import db
from aligons.db import DataSet, jgi, mask
from aligons.extern import htslib, jellyfish, kent
from aligons.util import cli, dl, fs, gff

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    gff.split_by_seqid(args.infile)


def fetch_and_bgzip(entry: DataSet, prefix: Path) -> list[cli.FuturePath]:
    url_prefix = entry["url_prefix"]
    species = entry["species"]
    annotation = entry["annotation"]
    sequences = entry["sequences"]
    softmasked = re.search(r"dna_sm|softmasked", sequences[0])
    dna_sm = "dna_sm" if softmasked else "dna"
    stem = f"{species}_{ver}" if (ver := entry.get("version", None)) else species
    out_fa = prefix / f"fasta/{species}/{stem}.{dna_sm}.toplevel.fa.gz"
    out_gff = prefix / f"gff3/{species}/{stem}.gff3.gz"
    fts: list[cli.FuturePath] = []
    if len(sequences) == 1:
        response_fa = dl_mirror_db(url_prefix + sequences[0])
        fts.append(cli.thread_submit(index_compress, response_fa, out_fa))
    elif fs.is_outdated(out_fa):
        chr_fa = [dl_mirror_db(url_prefix + s).content for s in sequences]
        content_fa = b"".join(chr_fa)
        fts.append(cli.thread_submit(bgzip_index, content_fa, out_fa))
    response_gff = dl_mirror_db(url_prefix + annotation)
    fts.append(cli.thread_submit(index_compress, response_gff, out_gff))
    return fts


def dl_mirror_db(url: str) -> dl.Response:
    if "jgi.doe.gov" in url:
        return jgi.session.mirror(url, db.path_mirror())
    return dl.mirror(url, db.path_mirror())


def index_as_completed(futures: list[cli.FuturePath]) -> list[cli.FuturePath]:
    fts: list[cli.FuturePath] = []
    for ft in confu.as_completed(futures):
        file = ft.result()
        if file.name.endswith("dna.toplevel.fa.gz"):
            fts.append(split_mask_index(file))
        elif file.name.endswith("dna_sm.toplevel.fa.gz"):
            fts.append(index_fasta([file]))
        elif file.name.endswith("gff3.gz"):
            fts.append(index_gff3([file]))
    return fts


def split_mask_index(toplevel_fa_gz: Path) -> cli.FuturePath:
    future_chromosomes = _split_toplevel_fa_work(toplevel_fa_gz)
    future_masked = [mask.submit(f) for f in future_chromosomes]
    links = [_symlink_masked(f) for f in future_masked]
    return index_fasta(links)


def _symlink_masked(ft: cli.FuturePath) -> Path:
    masked = ft.result()
    link = masked.parent.parent / masked.name
    return fs.symlink(masked, link, relative=True)


def index_compress(response: dl.Response, outfile: Path) -> Path:
    htslib.try_index(compress_lazy(response, outfile))
    return outfile


def bgzip_index(content: bytes, outfile: Path) -> Path:
    htslib.try_index(compress(content, outfile))
    return outfile


def index_fasta(paths: list[Path]) -> cli.FuturePath:
    """Create bgzipped and indexed genome.fa."""
    if len(paths) == 1:
        if "chromosome" in paths[0].name:
            _log.warning(f"splitting chromosome? {paths[0]}")
        paths = [f.result() for f in _split_toplevel_fa(paths[0])]
    return cli.thread_submit(_create_genome_bgzip, paths)


def index_gff3(paths: list[Path]) -> cli.FuturePath:  # gff3/{species}
    """Create bgzipped and indexed genome.gff3."""
    if len(paths) == 1:
        if "chromosome" in paths[0].name:
            _log.warning(f"splitting chromosome? {paths[0]}")
        paths = gff.split_by_seqid(paths[0])
    return cli.thread_submit(_create_genome_bgzip, paths)


def _create_genome_bgzip(files: list[Path]) -> Path:
    """Combine chromosome files and bgzip it."""
    files = fs.sorted_naturally(files)
    _log.debug(str(files))
    if cli.dry_run and not files:
        return Path(os.devnull)
    name = files[0].name
    ext = files[0].with_suffix("").suffix
    (outname, count) = re.subn(rf"\.chromosome\..+{ext}", rf".genome{ext}", name)
    assert count == 1, name
    outfile = files[0].with_name(outname)
    htslib.concat_bgzip(files, outfile)
    htslib.try_index(outfile)
    if outfile.name.endswith(".fa.gz"):
        kent.faSize(outfile)
    return outfile


def _split_toplevel_fa_work(fa_gz: Path) -> list[cli.FuturePath]:
    fmt = "{stem}.{seqid}.fa"
    outdir = fa_gz.with_name("_work")
    return htslib.split_fa_gz(fa_gz, fmt, (r"toplevel", "chromosome"), outdir)


def _split_toplevel_fa(fa_gz: Path) -> list[cli.FuturePath]:
    if "toplevel" not in fa_gz.name:
        _log.warning(f"'toplevel' not in {fa_gz.name}")
    fmt = "{stem}.{seqid}.fa.gz"
    return htslib.split_fa_gz(fa_gz, fmt, (r"toplevel", "chromosome"))


def softmask(species: str) -> cli.FuturePath:
    masked = jellyfish.run(species)
    return index_fasta(masked)


def recompress(infile: Path, outfile: Path) -> Path:
    if fs.is_outdated(outfile, infile) and not cli.dry_run:
        with infile.open("rb") as fin:
            content = fin.read()
        compress(content, outfile)
    return outfile


def compress_lazy(response: dl.Response, outfile: Path) -> Path:
    if fs.is_outdated(outfile, response.path) and not cli.dry_run:
        compress(response.content, outfile)
    return outfile


def compress(content: bytes, outfile: Path) -> Path:
    """Decompress/compress depending on file names.

    - .gff |> decompress |> sort |> bgzip
    - .bed |> decompress |> bgzip
    - .fa |> decompress |> bgzip
    - .zip |> decompress
    - gzip if outfile has new .gz
    """
    if not cli.dry_run and fs.is_outdated(outfile):
        if fs.is_zip(content):
            fs.expect_suffix(outfile, ".zip", negate=True)
            content = fs.zip_decompress(content)
        if htslib.to_be_bgzipped(outfile.name):
            fs.expect_suffix(outfile, ".gz")
            content = fs.gzip_decompress(content)
            if ".gff" in outfile.name:
                content = gff.sort(content)
            if outfile.name.endswith(".fa.gz") and not looks_like_fasta(content):
                msg = f"invalid fasta: {outfile} not written"
                raise ValueError(msg)
            content = htslib.bgzip_compress(content)
        elif outfile.suffix == ".gz":
            content = fs.gzip_compress(content)
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        with outfile.open("wb") as fout:
            fout.write(content)
    _log.info(f"{outfile}")
    return outfile


def looks_like_fasta(content: bytes) -> bool:
    if fs.is_gz(content):
        return True  # TODO: test first few bytes
    return content.startswith(b">")


if __name__ == "__main__":
    main()
