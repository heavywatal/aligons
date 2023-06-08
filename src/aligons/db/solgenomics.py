"""https://solgenomics.net/ftp/genomes/.

{local_db_root}/fasta/{species}/dna/{stem}.chromosome.{chr}.fa.gz
{local_db_root}/gff3/{species}/{stem}.gff3.gz
"""
import concurrent.futures as confu
import logging
import re
import tomllib
from collections.abc import Generator
from pathlib import Path
from typing import TypedDict

from aligons import db
from aligons.extern import htslib
from aligons.util import cli, resources_data

from . import tools

_log = logging.getLogger(__name__)
_futures: list[confu.Future[Path]] = []


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
    parser.add_argument("-D", "--download", action="store_true")
    args = parser.parse_args(argv or None)
    for entry in iter_dataset():
        if args.download:
            retrieve_deploy(entry)
        else:
            print(entry)
    for future in confu.as_completed(_futures):
        print(future.result())


def retrieve_deploy(entry: DataSet):
    sp = entry["species"]
    stem = f"{sp}_{ver}" if (ver := entry.get("version", None)) else sp
    query_gff = entry["annotation"]
    sequences = entry["sequences"]
    rel_fa_gz = f"fasta/{sp}/dna/{stem}.dna.toplevel.fa.gz"
    if len(sequences) > 1:
        chromosomes: list[Path] = []
        for query_fa in sequences:
            assert (mobj := re.search(r"([^.]+)", Path(query_fa).stem))
            seqid = mobj.group(1)
            tmppath = f"fasta/{sp}/dna/{stem}.dna.chromosome.{seqid}.fa.gz"
            chromosomes.append(_retrieve_bgzip(query_fa, tmppath))
        fa_gz = tools.cat(chromosomes, local_db_root() / rel_fa_gz)
    else:
        fa_gz = _retrieve_bgzip(sequences[0], rel_fa_gz)
    print(fa_gz)
    _futures.append(cli.ThreadPool().submit(htslib.faidx, fa_gz))
    gff3_gz = _retrieve_bgzip(query_gff, f"gff3/{sp}/{stem}.gff3.gz")
    print(gff3_gz)
    _futures.append(cli.ThreadPool().submit(htslib.tabix, gff3_gz))


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
