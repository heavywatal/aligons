import argparse
import concurrent.futures as confu
import itertools
import logging
import os
from pathlib import Path
from typing import Iterable

from .. import cli, fs, subp
from ..extern import htslib, kent
from . import ensemblgenomes, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-c", "--checksums", action="store_true")
    parser.add_argument("-i", "--index", action="store_true")
    parser.add_argument("clade")
    args = parser.parse_args(argv or None)
    cli.dry_run = args.dry_run
    cli.logging_config(args.loglevel)
    tree = phylo.trees[args.clade]
    species = phylo.extract_labels(tree)
    if args.download:
        download(species)
    if args.checksums:
        checksums(ensemblgenomes.rglob("CHECKSUMS", species))
    if args.index:
        index(species, jobs=args.jobs)
    if not any([args.download, args.checksums, args.index]):
        for x in ensemblgenomes.rglob("CHECKSUMS", species):
            print(x.parent)


def download(species: list[str]):
    assert species
    assert not (d := set(species) - set(ensemblgenomes.species_names_all())), d
    with ensemblgenomes.FTPensemblgenomes() as ftp:
        for sp in species:
            ftp.download_fasta(sp)
            ftp.download_gff3(sp)
    # for sp in species:
    #     options = "--include *_sm.chromosome.*.fa.gz --exclude *.gz"
    #     ensemblgenomes.rsync(f"fasta/{sp}/dna", options)
    #     options = "--include *.chromosome.*.gff3.gz --exclude *.gz"
    #     ensemblgenomes.rsync(f"gff3/{sp}", options)


def checksums(files: Iterable[Path]):
    with confu.ThreadPoolExecutor(max_workers=8) as pool:
        for file in files:
            _log.info(f"{file}")
            with file.open("rt") as fin:
                lines = fin.readlines()
            pool.map(checkline, lines, itertools.repeat(file.parent))


def checkline(line: str, directory: Path):
    (e_sum, e_blocks, name) = line.split()
    if (path := directory / name).exists():
        p = subp.run(f"sum {path}", stdout=subp.PIPE, text=True, quiet=True)
        (o_sum, o_blocks, _) = p.stdout.split()
        if (e_sum.lstrip("0"), e_blocks) != (o_sum, o_blocks):
            _log.error(f"{name}")
            _log.error(f"expected: {e_sum}\t{e_blocks}")
            _log.error(f"observed: {o_sum}\t{o_blocks}")


def index(species: list[str] = [], jobs: int | None = os.cpu_count()):
    with confu.ThreadPoolExecutor(max_workers=jobs) as pool:
        pool.map(index_fasta, ensemblgenomes.species_dirs("fasta", species))
        pool.map(index_gff3, ensemblgenomes.species_dirs("gff3", species))


def index_fasta(path: Path):  # fasta/{species}
    """Create bgzipped and indexed genome.fa."""
    if path.name != "dna":
        path /= "dna"
    for chromosome in fs.sorted_naturally(path.glob(r"*.chromosome.*.fa.gz")):
        print(kent.faToTwoBit(chromosome))
    genome = htslib.create_genome_bgzip(path)
    print(genome)
    print(htslib.faidx(genome))
    print(kent.faToTwoBit(genome))
    print(kent.faSize(genome))
    return genome


def index_gff3(path: Path):  # gff3/{species}
    """Create bgzipped and indexed genome.gff3."""
    genome = htslib.create_genome_bgzip(path)
    print(genome)
    print(htslib.tabix(genome))
    return genome


if __name__ == "__main__":
    main()
