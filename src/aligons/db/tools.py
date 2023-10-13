import logging
import re
from collections.abc import Iterator
from pathlib import Path

from aligons import db
from aligons.db import DataSet, api, mask
from aligons.extern import htslib, jellyfish, kent
from aligons.util import cli, dl, fs, gff, resources_data, tomllib

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    gff.split_by_seqid(args.infile)


def iter_dataset(filename: str) -> Iterator[DataSet]:
    with resources_data(filename).open("rb") as fin:
        meta = tomllib.load(fin)
    for dic in meta["dataset"]:
        if dic.get("draft", False):
            continue
        yield DataSet(dic)


def retrieve(entry: DataSet, prefix: Path) -> list[cli.FuturePath]:
    url_prefix = entry["url_prefix"]
    species = entry["species"]
    annotation = entry["annotation"]
    sequences = entry["sequences"]
    stem = f"{species}_{ver}" if (ver := entry.get("version", None)) else species
    out_fa = prefix / f"fasta/{species}/{stem}.dna.toplevel.fa.gz"
    out_gff = prefix / f"gff3/{species}/{stem}.gff3.gz"
    fts: list[cli.FuturePath] = []
    if not fs.is_outdated(out_fa):
        content_fa = b""
    elif len(sequences) > 1:
        chr_fa = [dl_mirror_db(url_prefix + s).content for s in sequences]
        content_fa = b"".join(chr_fa)
    else:
        content_fa = dl_mirror_db(url_prefix + sequences[0]).content
    fts.append(cli.thread_submit(bgzip_index, content_fa, out_fa))
    if not fs.is_outdated(out_gff):
        content_gff = b""
    else:
        content_gff = dl_mirror_db(url_prefix + annotation).content
    fts.append(cli.thread_submit(bgzip_index, content_gff, out_gff))
    return fts


def dl_mirror_db(url: str) -> dl.Response:
    return dl.mirror(url, db.path_mirror())


def prepare_fasta(species: str) -> cli.FuturePath:
    toplevel_fa_gz = api.get_file("*.dna.toplevel.fa.gz", species)
    future_chromosomes = _split_toplevel_fa_work(toplevel_fa_gz)
    future_masked = [mask.submit(f) for f in future_chromosomes]
    links = [_symlink_masked(f) for f in future_masked]
    return cli.thread_submit(index_fasta, links)


def _symlink_masked(ft: cli.FuturePath) -> Path:
    masked = ft.result()
    link = masked.parent.parent / masked.name
    return fs.symlink(masked, link)


def bgzip_index(content: bytes, outfile: Path):
    htslib.try_index(compress(content, outfile))
    return outfile


def index_fasta(paths: list[Path]):
    """Create bgzipped and indexed genome.fa."""
    if len(paths) == 1:
        paths = [f.result() for f in _split_toplevel_fa(paths[0])]
    fts = [cli.thread_submit(kent.faToTwoBit, x) for x in paths]
    genome = _create_genome_bgzip(paths)
    htslib.faidx(genome)
    kent.faToTwoBit(genome)
    kent.faSize(genome)
    cli.wait_raise(fts)
    return genome


def index_gff3(paths: list[Path]):  # gff3/{species}
    """Create bgzipped and indexed genome.gff3."""
    if len(paths) == 1:
        assert "chromosome" not in paths[0].name, paths[0]
        paths = gff.split_by_seqid(paths[0])
    genome = _create_genome_bgzip(paths)
    htslib.tabix(genome)
    return genome


def _create_genome_bgzip(files: list[Path]):
    """Combine chromosome files and bgzip it."""
    files = fs.sorted_naturally(files)
    _log.debug(str(files))
    if cli.dry_run and not files:
        return Path("/dev/null")
    name = files[0].name
    ext = files[0].with_suffix("").suffix
    (outname, count) = re.subn(rf"\.chromosome\..+{ext}", rf".genome{ext}", name)
    assert count == 1, name
    outfile = files[0].with_name(outname)
    return htslib.concat_bgzip(files, outfile)


def _split_toplevel_fa_work(fa_gz: Path) -> list[cli.FuturePath]:
    fmt = "{stem}.{seqid}.fa"
    outdir = fa_gz.with_name("_work")
    return htslib.split_fa_gz(fa_gz, fmt, (r"toplevel", "chromosome"), outdir)


def _split_toplevel_fa(fa_gz: Path) -> list[cli.FuturePath]:
    assert "toplevel" in fa_gz.name, fa_gz
    fmt = "{stem}.{seqid}.fa.gz"
    return htslib.split_fa_gz(fa_gz, fmt, (r"toplevel", "chromosome"))


def softmask(species: str):
    masked = jellyfish.run(species)
    return index_fasta(masked)


def compress(content: bytes, outfile: Path) -> Path:
    """Uncompress/compress depending on file names.

    - .gff |> uncompress |> sort |> bgzip
    - .bed |> uncompress |> bgzip
    - .fa |> uncompress |> bgzip
    - .zip |> uncompress
    - gzip if outfile has new .gz
    """
    if not cli.dry_run and fs.is_outdated(outfile):
        if fs.is_zip(content):
            assert outfile.suffix != ".zip", outfile
            content = fs.zip_decompress(content)
        if htslib.to_be_bgzipped(outfile.name):
            assert outfile.suffix == ".gz", outfile
            content = fs.gzip_decompress(content)
            if ".gff" in outfile.name:
                content = gff.sort(content)
            content = htslib.bgzip_compress(content)
        elif outfile.suffix == ".gz":
            content = fs.gzip_compress(content)
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        with outfile.open("wb") as fout:
            fout.write(content)
    _log.info(f"{outfile}")
    return outfile


if __name__ == "__main__":
    main()
