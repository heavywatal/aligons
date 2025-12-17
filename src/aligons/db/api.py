"""Interface to pre-processed data.

Raw data from various sources located at `{db.root}/*` are processed and stored at
`{db.root}/aligons/{origin}/{species}/`. This module provides functions to access
those processed data files.
"""

import csv
import functools
import logging
import re
from contextlib import suppress
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator
    from pathlib import Path

from aligons.util import cli, config, fs

from . import _rsrc, phylo, plantregmap

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    """CLI for manual execution and testing."""
    parser = cli.ArgumentParser()
    parser.add_argument("-S", "--species")
    parser.add_argument("-C", "--clade", default="angiospermae")
    args = parser.parse_args(argv or None)
    if args.species:
        print_existing(args.species)
    else:
        print_stats(args.clade)


def print_existing(species: str) -> None:
    """Print existing files for the given species."""
    for p in _glob("*", species):
        fs.print_if_exists(p)


def shorten(species: str) -> str:
    """Get a short name for the given species."""
    try:
        return plantregmap.shorten(species)
    except KeyError:
        return phylo.shorten(species)


@functools.cache
def species_names() -> list[str]:
    """Collect available species names in the database."""
    return [x.name for x in _species_dirs()]


def fasize(species: str) -> Path:
    """Get `fasize.chrom.sizes` file for the given species."""
    return get_file("fasize.chrom.sizes", species)


def genome_fai(species: str) -> Path:
    """Get the genome FASTA index file for the given species."""
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file("*.dna_sm.genome.fa.gz.fai", species, subdir)


def genome_fa(species: str) -> Path:
    """Get the genome FASTA file for the given species."""
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file("*.dna_sm.genome.fa.gz", species, subdir)


def genome_gff3(species: str) -> Path:
    """Get the genome GFF3 annotation file for the given species."""
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file("*.genome.gff3.gz", species, subdir)


def chromosome_2bit(species: str, seqid: str) -> Path:
    """Get the chromosome 2bit file for the given species."""
    subdir = "kmer" if config["db"]["kmer"] else ""
    return get_file(f"*.chromosome.{seqid}.2bit", species, subdir)


def iter_chromosome_2bit(species: str) -> Iterable[Path]:
    """Iterate over all chromosome 2bit files for the given species."""
    return _glob("*.chromosome.*.2bit", species)


def get_file(pattern: str, species: str, subdir: str = "") -> Path:
    """Get a unique file matching the pattern for the given species.

    :param pattern: A glob pattern to match file names.
    :param species: Species name.
    :param subdir: A subdirectory under species directory.
    :returns: The unique matching file path:
        `{db.root}/aligons/{origin}/{species}/{subdir}/{pattern}`.
    :raises FileNotFoundError: If no matching file is found.
    :raises ValueError: If multiple matching files are found.
    """
    found = list(_glob(pattern, species, subdir))
    if not found:
        msg = f"{pattern} not found in {species}/{subdir}"
        raise FileNotFoundError(msg)
    if len(found) > 1:
        msg = f"{pattern} is not unique in {species}/{subdir}: {found}"
        raise ValueError(msg)
    return found[0]


def sanitize_queries(target: str, queries: list[str]) -> list[str]:
    """Make sure all query species exist in the database.

    :param target: Target species name.
    :param queries: Query species names. Target species may be included.
    :returns: Validated query species names without duplicates and target species.
    :raises ValueError: If any query species is not found in the database.
    """
    not_in_db: list[str] = []
    for query in queries:
        it = _species_dirs(query)
        if next(it, None) is None:
            not_in_db.append(query)
    if not_in_db:
        msg = f"{not_in_db} not in {species_names()}"
        raise ValueError(msg)
    query_set = set(queries)
    with suppress(KeyError):
        query_set.remove(target)
    _log.debug(f"{query_set=}")
    return list(query_set)


def _species_dirs(species: str = "") -> Iterator[Path]:
    """Iterate over species directories.

    :param species: Species name to filter if not empty.
    :returns: An iterator over directories: `{db.root}/aligons/{origin}/{species}/`.
    """
    if (path := prefix(species)).exists():  # {db.root}/aligons/{species}
        yield path
    for prefix_origin in _iter_prefix():
        if not prefix_origin.exists():
            continue
        for path in prefix_origin.iterdir():
            if path.is_dir():
                if species and (species.lower() != path.name.lower()):
                    continue
                yield path


def _glob(pattern: str, species: str, subdir: str = "") -> Iterable[Path]:
    """Find files matching the pattern under species directories.

    :param pattern: A glob pattern to match file names.
    :param species: Species name.
    :param subdir: A subdirectory under species directory.
    :returns: An iterator over matching file paths:
        `{db.root}/aligons/{origin}/{species}/{subdir}/{pattern}`.
    """
    is_enough = False  # to allow duplicated species from multiple origins
    for sp_dir in _species_dirs(species):
        d = sp_dir / subdir
        for x in fs.sorted_naturally(d.glob(pattern)):
            is_enough = True
            yield x
        if is_enough:
            break


def prefix(relpath: str | Path = "") -> Path:
    """Get full path under the aligons database root.

    :param relpath: Relative path under `{db.root}/aligons/`.
    :returns: Full path `{db.root}/aligons/{relpath}`.
    """
    return _rsrc.db_root("aligons") / relpath


def _iter_prefix() -> Iterable[Path]:
    """Iterate over database origin directories.

    :returns: An iterator over `{db.root}/aligons/{origin}/` paths.
    """
    return (prefix(origin) for origin in _iter_db_origin())


def _iter_db_origin() -> Iterable[str]:
    """Iterate over `config["db"]["origin"]` with version formatting."""
    for origin in config["db"]["origin"]:
        if origin.startswith("ensembl"):
            yield origin.format(version=config["ensemblgenomes"]["version"])
        elif origin.startswith("phytozome"):
            v = config["jgi"]["organism"].removeprefix("PhytozomeV")
            yield origin.format(version=v)
        else:
            yield origin


def print_stats(clade: str) -> None:
    """Print a tree with genome statistics for a given clade."""
    newick = phylo.get_subtree([clade])
    root = phylo.parse_newick(newick)
    for pre, species in phylo.rectangular(phylo.render_tips(root, [])):
        words: list[str] = [pre, species]
        try:
            chrom_sizes_ = chrom_sizes(species)
            fasize = sum(chrom_sizes_.values()) / 1e6
            n_seqs = len(chrom_sizes_)
            words.append(f"{fasize:4.0f}Mbp")
            words.append(f"{n_seqs:2}")
            gff_size = _gff3_size(species)
            words.append(f"{gff_size:4.1f}MB")
            sources = ",".join(_sources(species))
            words.append(f"{sources}")
        except FileNotFoundError:
            pass
        _log.info(" ".join(words))


def chrom_sizes(species: str) -> dict[str, int]:
    """Return a dictionary of chromosome sizes for the given species."""
    path = genome_fai(species)
    with path.open() as fin:
        reader = csv.reader(fin, delimiter="\t")
        items = {x[0]: int(x[1]) for x in reader}
    return {k: v for k, v in items.items() if not k.lower().startswith("scaffold")}


def _gff3_size(species: str) -> float:
    path = genome_gff3(species)
    return path.stat().st_size / 1e6


def _sources(species: str) -> Iterable[str]:
    for path in _species_dirs(species):
        if path.glob("*.genome.*.gz"):
            yield _nickname(path.parent.name)


def _nickname(x: str) -> str:
    x = re.sub(r"ensembl(-compara)?-(\d+)", r"en\1\2", x)
    x = re.sub(r"[pP]hytozome[-V](\d+)", r"pz\1", x)
    mapping = {
        "plantregmap": "prm",
        "solgenomics": "sol",
    }
    return mapping.get(x, x)


if __name__ == "__main__":
    main()
