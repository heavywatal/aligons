"""https://solgenomics.net/ftp/genomes/."""
import logging
from collections.abc import Generator
from pathlib import Path
from typing import TypedDict

from aligons import db
from aligons.db import api, mask
from aligons.extern import htslib
from aligons.util import cli, fs, resources_data, tomllib

from . import tools

_log = logging.getLogger(__name__)


class DataSet(TypedDict):
    species: str
    version: str
    draft: bool
    label: str
    clade: str
    sequences: list[str]
    annotation: str


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-M", "--mask", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    args = parser.parse_args(argv or None)
    if args.download:
        fts: list[cli.FuturePath] = []
        for entry in iter_dataset():
            fts.extend(retrieve(entry))
        cli.wait_raise(fts)
    if args.mask:
        fts: list[cli.FuturePath] = []
        for entry in iter_dataset():
            species = entry["species"]
            fts.append(prepare_fasta(species))
            gff3_gz = api.get_file("*.gff3.gz", species)
            fts.append(cli.thread_submit(tools.index_gff3, [gff3_gz]))
        cli.wait_raise(fts)


def prepare_fasta(species: str) -> cli.FuturePath:
    toplevel_fa_gz = api.get_file("*.dna.toplevel.fa.gz", species)
    future_chromosomes = _split_toplevel_fa(toplevel_fa_gz)
    future_masked = [mask.submit(f) for f in future_chromosomes]
    links = [_symlink_masked(f) for f in future_masked]
    return cli.thread_submit(tools.index_fasta, links)


def _symlink_masked(ft: cli.FuturePath) -> Path:
    masked = ft.result()
    link = masked.parent.parent / masked.name
    return fs.symlink(masked, link)


def _split_toplevel_fa(fa_gz: Path) -> list[cli.FuturePath]:
    fmt = "{stem}.{seqid}.fa"
    outdir = fa_gz.parent / "_work"
    return htslib.split_fa_gz(fa_gz, fmt, (r"toplevel", "chromosome"), outdir)


def retrieve(entry: DataSet) -> list[cli.FuturePath]:
    remote_prefix = "https://solgenomics.net/ftp/genomes/"
    species = entry["species"]
    annotation = entry["annotation"]
    sequences = entry["sequences"]
    stem = f"{species}_{ver}" if (ver := entry.get("version", None)) else species
    out_fa = db_prefix() / f"fasta/{species}/{stem}.dna.toplevel.fa.gz"
    out_gff = db_prefix() / f"gff3/{species}/{stem}.gff3.gz"
    fts: list[cli.FuturePath] = []
    if not fs.is_outdated(out_fa):
        content_fa = b""
    elif len(sequences) > 1:
        chr_fa = [tools.retrieve_content(remote_prefix + s) for s in sequences]
        content_fa = b"".join(chr_fa)
    else:
        content_fa = tools.retrieve_content(remote_prefix + sequences[0])
    fts.append(cli.thread_submit(bgzip_index, content_fa, out_fa))
    if not fs.is_outdated(out_gff):
        content_gff = b""
    else:
        content_gff = tools.retrieve_content(remote_prefix + annotation)
    fts.append(cli.thread_submit(bgzip_index, content_gff, out_gff))
    return fts


def bgzip_index(content: bytes, outfile: Path):
    htslib.try_index(tools.compress(content, outfile))
    return outfile


def iter_dataset() -> Generator[DataSet, None, None]:
    with resources_data("solgenomics.toml").open("rb") as fin:
        meta = tomllib.load(fin)
    for dic in meta["dataset"]:
        if dic.get("draft", False):
            continue
        yield DataSet(dic)


def db_prefix():
    return db.path("solgenomics")


if __name__ == "__main__":
    main()
