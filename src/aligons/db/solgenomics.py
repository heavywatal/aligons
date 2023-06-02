"""https://solgenomics.net/ftp/genomes/.

{local_db_root}/fasta/{species}/dna/{stem}.chromosome.{chr}.fa.gz
{local_db_root}/gff3/{species}/{stem}.gff3.gz
"""
import logging
import re
import tomllib
from collections.abc import Generator
from pathlib import Path
from typing import TypedDict

from aligons import db
from aligons.util import cli, resources_data

from . import tools

_log = logging.getLogger(__name__)


class DataSet(TypedDict):
    species: str
    version: str
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
            download(entry)
        else:
            print(entry)


def download(entry: DataSet):
    sp = entry["species"]
    stem = f"{sp}_{ver}" if (ver := entry.get("version", None)) else sp
    sequences = entry["sequences"]
    for query_fa in sequences:
        if len(sequences) > 1:
            assert (mobj := re.search(r"chr([^.]+)", Path(query_fa).stem))
            chromo = mobj.group(1)
            relpath = f"fasta/{sp}/dna/{stem}.chromosome.{chromo}.fa.gz"
        else:
            relpath = f"fasta/{sp}/dna/{stem}.fa.gz"
        print(retrieve(query_fa, relpath))
    query_gff: str = entry["annotation"]
    relpath = f"gff3/{sp}/{stem}.gff3.gz"
    print(retrieve(query_gff, relpath))


def retrieve(query: str, relpath: str):
    outfile = local_db_root() / relpath
    url = f"https://solgenomics.net/ftp/genomes/{query}"
    return tools.retrieve_bgzip(url, outfile)


def iter_dataset() -> Generator[DataSet, None, None]:
    with resources_data("solgenomics.toml").open("rb") as fin:
        meta = tomllib.load(fin)
    for dic in meta["dataset"]:
        yield DataSet(dic)


def local_db_root():
    return db.path("solgenomics")


if __name__ == "__main__":
    main()
