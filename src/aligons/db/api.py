"""Interface to working data.

prefix: {db.root}/{origin}/

- fasta/{species}/
- gff3/{species}/
"""
import csv
import functools
import logging
import re
from collections.abc import Iterable
from contextlib import suppress
from pathlib import Path

from aligons import db
from aligons.db import ensemblgenomes, phylo, plantregmap
from aligons.util import cli, config, fs

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-S", "--species")
    parser.add_argument("-C", "--clade", default="angiospermae")
    args = parser.parse_args(argv or None)
    if args.species:
        for p in _glob("*", args.species):
            print(p)
    else:
        print_stats(args.clade)


def shorten(species: str) -> str:
    try:
        return plantregmap.shorten(species)
    except KeyError:
        return phylo.shorten(species)


@functools.cache
def species_names(fmt: str = "fasta") -> list[str]:
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


def get_file(pattern: str, species: str, subdir: str = "") -> Path:
    found = list(_glob(pattern, species, subdir))
    if not found:
        msg = f"{pattern} not found in {species}/{subdir}"
        raise FileNotFoundError(msg)
    if len(found) > 1:
        msg = f"{pattern} is not unique in {species}/{subdir}: {found}"
        raise ValueError(msg)
    return found[0]


def get_file_nolabel(pattern: str, species: str, subdir: str = "") -> Path:
    it = _glob(pattern, species, subdir)
    rex = re.compile(r"\.(chromosome|genome|toplevel)\.")
    found = [x for x in it if not rex.search(x.name)]
    if not found:
        msg = f"{pattern} not found in {species}/{subdir}"
        raise FileNotFoundError(msg)
    if len(found) > 1:
        msg = f"{pattern} is not unique in {species}/{subdir}: {found}"
        raise ValueError(msg)
    return found[0]


def sanitize_queries(target: str, queries: list[str]) -> list[str]:
    query_set = set(queries)
    with suppress(KeyError):
        query_set.remove(target)
    _log.debug(f"{query_set=}")
    not_in_db = query_set - set(species_names())
    if not_in_db:
        msg = f"{not_in_db} not in {species_names()}"
        raise ValueError(msg)
    return list(query_set)


def _species_dirs(fmt: str = "fasta") -> Iterable[Path]:
    for prefix in _iter_prefix():
        if not (fmt_dir := prefix / fmt).exists():
            _log.warning(f"{fmt_dir} does not exist")
            continue
        for path in fmt_dir.iterdir():
            if path.is_dir():
                yield path


def _glob(pattern: str, species: str, subdir: str = "") -> Iterable[Path]:
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


def print_stats(clade: str) -> None:
    newick = phylo.get_subtree([clade])
    root = phylo.parse_newick(newick)
    for pre, species in phylo.rectangulate(phylo.render_tips(root, [])):
        print(f"{pre} {species}", end="")
        try:
            chrom_sizes_ = chrom_sizes(species)
            fasize = sum(chrom_sizes_.values()) / 1e6
            nseqs = len(chrom_sizes_)
            print(f" {fasize:4.0f}Mbp {nseqs:2}", end="")
            gffsize = _gff3_size(species)
            print(f" {gffsize:4.1f}MB")
        except FileNotFoundError:
            print("")


def chrom_sizes(species: str) -> dict[str, int]:
    path = fasize(species)
    with path.open() as fin:
        reader = csv.reader(fin, delimiter="\t")
        return {x[0]: int(x[1]) for x in reader}


def _gff3_size(species: str) -> float:
    path = genome_gff3(species)
    return path.stat().st_size / 1e6


if __name__ == "__main__":
    main()
