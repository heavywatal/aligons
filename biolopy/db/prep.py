import concurrent.futures as confu
import logging
import os
from pathlib import Path

from ..extern import htslib, kent
from ..util import cli, fs
from . import ensemblgenomes, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("clade")
    args = parser.parse_args(argv or None)
    cli.dry_run = args.dry_run
    cli.logging_config(args.loglevel)
    tree = phylo.newicks[args.clade]
    species = phylo.extract_names(tree)
    if args.download:
        download(species, jobs=args.jobs)
    else:
        index(species, jobs=args.jobs)


def download(species: list[str], jobs: int | None = os.cpu_count()):
    assert species
    assert not (d := set(species) - set(ensemblgenomes.species_names_all())), d
    with (
        confu.ThreadPoolExecutor(max_workers=jobs) as pool,
        ensemblgenomes.FTPensemblgenomes() as ftp,
    ):
        for sp in species:
            dir = ftp.download_fasta(sp)
            pool.submit(index_fasta, dir)
            dir = ftp.download_gff3(sp)
            pool.submit(index_gff3, dir)
    # for sp in species:
    #     options = "--include *_sm.chromosome.*.fa.gz --exclude *.gz"
    #     ensemblgenomes.rsync(f"fasta/{sp}/dna", options)
    #     options = "--include *.chromosome.*.gff3.gz --exclude *.gz"
    #     ensemblgenomes.rsync(f"gff3/{sp}", options)


def index(species: list[str] = [], jobs: int | None = os.cpu_count()):
    with confu.ThreadPoolExecutor(max_workers=jobs) as pool:
        pool.map(index_fasta, ensemblgenomes.species_dirs("fasta", species))
        pool.map(index_gff3, ensemblgenomes.species_dirs("gff3", species))


def index_fasta(path: Path):  # fasta/{species}
    """Create bgzipped and indexed genome.fa."""
    if path.name != "dna":
        path /= "dna"
    for chromosome in fs.sorted_naturally(path.glob(r"*.chromosome.*.fa.gz")):
        _log.info(str(kent.faToTwoBit(chromosome)))
    genome = htslib.create_genome_bgzip(path)
    _log.info(str(genome))
    _log.info(str(htslib.faidx(genome)))
    _log.info(str(kent.faToTwoBit(genome)))
    _log.info(str(kent.faSize(genome)))
    print(path)
    return genome


def index_gff3(path: Path):  # gff3/{species}
    """Create bgzipped and indexed genome.gff3."""
    genome = htslib.create_genome_bgzip(path)
    _log.info(str(genome))
    _log.info(str(htslib.tabix(genome)))
    print(path)
    return genome


if __name__ == "__main__":
    main()
