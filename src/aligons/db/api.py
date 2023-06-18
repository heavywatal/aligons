"""Interface to working data.

prefix: {db.root}/{origin}/

- fasta/{species}/
- gff3/{species}/
"""
import functools
import logging
import re
from collections.abc import Iterable
from contextlib import suppress
from pathlib import Path

from aligons import db
from aligons.db import ensemblgenomes
from aligons.util import cli, config, fs

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-S", "--species")
    args = parser.parse_args(argv or None)
    if args.species:
        for p in _glob("*", args.species):
            print(p)
    else:
        for p in _species_dirs():
            print(p)


@functools.cache
def species_names(fmt: str = "fasta"):
    return [x.name for x in _species_dirs(fmt)]


def fasize(species: str) -> Path:
    return get_file("fasize.chrom.sizes", species)


def genome_fa(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file("*.genome.fa.gz", species, subdir)


def genome_2bit(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file("*.genome.2bit", species, subdir)


def genome_gff3(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file("*.genome.gff3.gz", species, subdir)


def list_chromosome_fa(species: str) -> Iterable[Path]:
    return fs.sorted_naturally(_glob("*.chromosome.*.fa.gz", species))


def list_chromosome_2bit(species: str) -> list[Path]:
    patt = "*.chromosome.*.2bit"
    subdir = "kmer" if config["db"]["kmer"] else ""
    return fs.sorted_naturally(_glob(patt, species, subdir))


def list_chromosome_gff3(species: str) -> Iterable[Path]:
    return fs.sorted_naturally(_glob("*.chromosome*.gff3.gz", species))


def get_file(pattern: str, species: str, subdir: str = ""):
    found = list(_glob(pattern, species, subdir))
    assert found, f"not found {pattern} in {species}/{subdir}"
    assert len(found) == 1, found
    return found[0]


def get_file_nolabel(pattern: str, species: str, subdir: str = ""):
    it = _glob(pattern, species, subdir)
    rex = re.compile(r"\.(chromosome|genome|toplevel)\.")
    found = [x for x in it if not rex.search(x.name)]
    assert found, f"not found {pattern} in {species}/{subdir}"
    assert len(found) == 1, found
    return found[0]


def sanitize_queries(target: str, queries: list[str]):
    queries = list(dict.fromkeys(queries))
    with suppress(ValueError):
        queries.remove(target)
    assert queries, target
    _log.debug(f"{queries=}")
    assert set(queries) <= set(species_names()), f"{queries} vs {species_names()}"
    return queries


def _species_dirs(fmt: str = "fasta") -> Iterable[Path]:
    for prefix in _iter_prefix():
        if not (fmt_dir := prefix / fmt).exists():
            _log.warning(f"{fmt_dir} does not exist")
            continue
        for path in fmt_dir.iterdir():
            if path.is_dir():
                yield path


def _glob(pattern: str, species: str, subdir: str = ""):
    is_enough = False  # to allow duplicated species from multiple origins
    for prefix in _iter_prefix():
        for fmt in ("fasta", "gff3"):
            d = prefix / fmt / species / subdir
            for x in fs.sorted_naturally(d.glob(pattern)):
                is_enough = True
                yield x
        if is_enough:
            break


def _iter_prefix() -> Iterable[Path]:
    return (_prefix(origin) for origin in config["db"]["origin"])


def _prefix(origin: str) -> Path:
    return db.path(origin.format(version=ensemblgenomes.version()))


if __name__ == "__main__":
    main()
