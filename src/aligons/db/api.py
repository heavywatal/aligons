"""Interface to working data.

{db.root}/aligons/{origin}/{species}/
"""

import csv
import functools
import logging
import re
from collections.abc import Iterable
from contextlib import suppress
from pathlib import Path

from aligons.util import cli, config, fs

from . import _rsrc, ensemblgenomes, phylo, plantregmap

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-S", "--species")
    parser.add_argument("-C", "--clade", default="angiospermae")
    args = parser.parse_args(argv or None)
    if args.species:
        for p in _glob("*", args.species):
            fs.print_if_exists(p)
    else:
        print_stats(args.clade)


def shorten(species: str) -> str:
    try:
        return plantregmap.shorten(species)
    except KeyError:
        return phylo.shorten(species)


@functools.cache
def species_names() -> list[str]:
    return [x.name for x in _species_dirs()]


def fasize(species: str) -> Path:
    return get_file("fasize.chrom.sizes", species)


def genome_fa(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file("*.dna_sm.genome.fa.gz", species, subdir)


def genome_gff3(species: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file("*.genome.gff3.gz", species, subdir)


def chromosome_2bit(species: str, seqid: str) -> Path:
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file(f"*.chromosome.{seqid}.2bit", species, subdir)


def iter_chromosome_2bit(species: str) -> Iterable[Path]:
    return _glob("*.chromosome.*.2bit", species)


def get_file(pattern: str, species: str, subdir: str = "") -> Path:
    found = list(_glob(pattern, species, subdir))
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


def _species_dirs(species: str = "") -> Iterable[Path]:
    for prefix in _iter_prefix():
        for path in prefix.iterdir():
            if path.is_dir():
                if species and (species.lower() != path.name.lower()):
                    continue
                yield path


def _glob(pattern: str, species: str, subdir: str = "") -> Iterable[Path]:
    is_enough = False  # to allow duplicated species from multiple origins
    for sp_dir in _species_dirs(species):
        d = sp_dir / subdir
        for x in fs.sorted_naturally(d.glob(pattern)):
            is_enough = True
            yield x
        if is_enough:
            break


def prefix(relpath: str | Path = "") -> Path:
    return _rsrc.db_root("aligons") / relpath


def _iter_prefix() -> Iterable[Path]:
    return (prefix(origin) for origin in _iter_db_origin())


def _iter_db_origin() -> Iterable[str]:
    v = ensemblgenomes.version()
    return (origin.format(version=v) for origin in config["db"]["origin"])


def print_stats(clade: str) -> None:
    newick = phylo.get_subtree([clade])
    root = phylo.parse_newick(newick)
    for pre, species in phylo.rectangulate(phylo.render_tips(root, [])):
        words: list[str] = [pre, species]
        try:
            chrom_sizes_ = chrom_sizes(species)
            fasize = sum(chrom_sizes_.values()) / 1e6
            nseqs = len(chrom_sizes_)
            words.append(f"{fasize:4.0f}Mbp")
            words.append(f"{nseqs:2}")
            gffsize = _gff3_size(species)
            words.append(f"{gffsize:4.1f}MB")
            sources = ",".join(_sources(species))
            words.append(f"{sources}")
        except FileNotFoundError:
            pass
        _log.info(" ".join(words))


def chrom_sizes(species: str) -> dict[str, int]:
    path = fasize(species)
    with path.open() as fin:
        reader = csv.reader(fin, delimiter="\t")
        return {x[0]: int(x[1]) for x in reader}


def _gff3_size(species: str) -> float:
    path = genome_gff3(species)
    return path.stat().st_size / 1e6


def _sources(species: str) -> Iterable[str]:
    for path in _species_dirs(species):
        if path.glob("*.genome.*.gz"):
            yield _nickname(path.parent.name)


def _nickname(x: str) -> str:
    x = re.sub(r"ensembl(-compara)?-(\d+)", r"en\1\2", x)
    mapping = {
        "plantregmap": "prm",
        "solgenomics": "sol",
    }
    return mapping.get(x, x)


if __name__ == "__main__":
    main()
