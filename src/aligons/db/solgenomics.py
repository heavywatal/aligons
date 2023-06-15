"""https://solgenomics.net/ftp/genomes/.

{local_db_root}/fasta/{species}/dna/{stem}.chromosome.{chr}.fa.gz
{local_db_root}/gff3/{species}/{stem}.gff3.gz
"""
import concurrent.futures as confu
import logging
import re
from collections.abc import Generator
from pathlib import Path
from typing import TypedDict

from aligons import db
from aligons.db import mask
from aligons.extern import htslib
from aligons.util import cli, resources_data, tomllib

from . import tools

_log = logging.getLogger(__name__)
_futures: list[cli.FuturePath] = []


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
    for entry in iter_dataset():
        if args.download:
            retrieve_deploy(entry)
    cli.wait_raise(_futures)
    _futures.clear()
    if args.mask:
        split_toplevel_fa()
        fts: list[cli.FuturePath] = []
        for ft in confu.as_completed(_futures):
            infile = ft.result()
            assert (mobj := re.match("[^_]+_[^_]+", infile.name))
            species = mobj.group(0)
            genus = species.split("_")[0]
            fts.append(mask.run(ft.result(), genus))
        cli.wait_raise(fts)


def split_toplevel_fa():
    pattern = "*.dna.toplevel.fa.gz"
    for fa_gz in local_db_root().rglob(pattern):
        _futures.extend(_split_toplevel_fa(fa_gz))


def _split_toplevel_fa(fa_gz: Path) -> list[cli.FuturePath]:
    fmt = "{stem}.{seqid}.fa"
    outdir = fa_gz.parent / "_work"
    return htslib.split_fa_gz(fa_gz, fmt, (r"toplevel", "chromosome"), outdir)


def retrieve_deploy(entry: DataSet):
    sp = entry["species"]
    stem = f"{sp}_{ver}" if (ver := entry.get("version", None)) else sp
    query_gff = entry["annotation"]
    sequences = entry["sequences"]
    rel_fa_gz = f"fasta/{sp}/dna/{stem}.dna.toplevel.fa.gz"
    if len(sequences) > 1:
        ft_chroms: list[cli.FuturePath] = []
        for query_fa in sequences:
            assert (mobj := re.search(r"([^.]+)", Path(query_fa).stem))
            seqid = mobj.group(1)
            tmppath = f"fasta/{sp}/dna/{stem}.dna.chromosome.{seqid}.fa.gz"
            ft_chroms.append(_retrieve_bgzip(query_fa, tmppath))
        abs_fa_gz = local_db_root() / rel_fa_gz
        fa_gz = cli.thread_submit(tools.cat, ft_chroms, abs_fa_gz)
    else:
        fa_gz = _retrieve_bgzip(sequences[0], rel_fa_gz)
    _futures.append(cli.thread_submit(htslib.faidx, fa_gz))
    gff3_gz = _retrieve_bgzip(query_gff, f"gff3/{sp}/{stem}.gff3.gz")
    _futures.append(cli.thread_submit(htslib.tabix, gff3_gz))


def _retrieve_bgzip(query: str, relpath: str):
    outfile = local_db_root() / relpath
    url = f"https://solgenomics.net/ftp/genomes/{query}"
    return tools.retrieve_compress(url, outfile)


def iter_dataset() -> Generator[DataSet, None, None]:
    with resources_data("solgenomics.toml").open("rb") as fin:
        meta = tomllib.load(fin)
    for dic in meta["dataset"]:
        if dic.get("draft", False):
            continue
        yield DataSet(dic)


def local_db_root():
    return db.path("solgenomics")


if __name__ == "__main__":
    main()
