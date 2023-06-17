"""Interface to working data.

prefix: {db.root}/aligons/{db.label}/

- fasta/{species}/
- gff3/{species}/
"""
import functools
import logging
from collections.abc import Iterable
from contextlib import suppress
from pathlib import Path

from aligons import db
from aligons.util import cli, config, fs

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-S", "--species")
    args = parser.parse_args(argv or None)
    print(prefix())
    if args.species:
        for p in _glob("*", args.species):
            print(p.relative_to(prefix()))
    else:
        for p in _species_dirs():
            print(p)


def prefix() -> Path:
    return db.path("aligons") / config["db"]["label"]


@functools.cache
def species_names(fmt: str = "fasta"):
    return [x.name for x in _species_dirs(fmt)]


def fasize(species: str) -> Path:
    return _get_file("fasize.chrom.sizes", species)


def genome_fa(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return _get_file("*.genome.fa.gz", species, subdir)


def genome_2bit(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return _get_file("*.genome.2bit", species, subdir)


def genome_gff3(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return _get_file("*.genome.gff3.gz", species, subdir)


def list_chromosome_fa(species: str) -> Iterable[Path]:
    return fs.sorted_naturally(_glob("*.chromosome.*.fa.gz", species))


def list_chromosome_2bit(species: str) -> list[Path]:
    patt = "*.chromosome.*.2bit"
    subdir = "kmer" if config["db"]["kmer"] else ""
    return fs.sorted_naturally(_glob(patt, species, subdir))


def list_chromosome_gff3(species: str) -> Iterable[Path]:
    return fs.sorted_naturally(_glob("*.chromosome*.gff3.gz", species))


def sanitize_queries(target: str, queries: list[str]):
    queries = list(dict.fromkeys(queries))
    with suppress(ValueError):
        queries.remove(target)
    assert queries
    _log.debug(f"{queries=}")
    assert set(queries) <= set(species_names())
    return queries


def _species_dirs(fmt: str = "fasta") -> Iterable[Path]:
    assert (root := prefix() / fmt).exists(), root
    for path in root.iterdir():
        if path.is_dir():
            yield path


def _get_file(pattern: str, species: str, subdir: str = ""):
    found = list(_glob(pattern, species, subdir))
    assert len(found) == 1, found
    return found[0]


def _glob(pattern: str, species: str, subdir: str = ""):
    for fmt in ("fasta", "gff3"):
        d = prefix() / fmt / species / subdir
        for x in fs.sorted_naturally(d.glob(pattern)):
            yield x


if __name__ == "__main__":
    main()
