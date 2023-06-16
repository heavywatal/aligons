import logging
from collections.abc import Iterable
from contextlib import suppress
from pathlib import Path

from aligons.util import config, fs

from . import ensemblgenomes

_log = logging.getLogger(__name__)


def species_names() -> list[str]:
    return ensemblgenomes.species_names()


def species_dirs(
    fmt: str = "fasta", species: list[str] | None = None
) -> Iterable[Path]:
    return ensemblgenomes.species_dirs(fmt, species)


def fasize(species: str) -> Path:
    return ensemblgenomes.get_file("fasize.chrom.sizes", species)


def genome_fa(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return ensemblgenomes.get_file("*.genome.fa.gz", species, subdir)


def genome_2bit(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return ensemblgenomes.get_file("*.genome.2bit", species, subdir)


def genome_gff3(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return ensemblgenomes.get_file("*.genome.gff3.gz", species, subdir)


def iter_chromosome_fa(species: str) -> Iterable[Path]:
    return ensemblgenomes.glob("*.chromosome.*.fa.gz", [species])


def list_chromosome_2bit(species: str) -> list[Path]:
    patt = "*.chromosome.*.2bit"
    subdir = "kmer" if config["db"]["kmer"] else ""
    it = ensemblgenomes.glob(patt, [species], subdir)
    return fs.sorted_naturally(it)


def iter_chromosome_gff3(species: str) -> Iterable[Path]:
    return ensemblgenomes.glob("*.chromosome*.gff3.gz", [species])


def sanitize_queries(target: str, queries: list[str]):
    queries = list(dict.fromkeys(queries))
    with suppress(ValueError):
        queries.remove(target)
    assert queries
    _log.debug(f"{queries=}")
    assert set(queries) <= set(species_names())
    return queries
