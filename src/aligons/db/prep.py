import logging
import re
from pathlib import Path

from aligons.extern import htslib, jellyfish, kent
from aligons.util import cli, config, fs, read_config

from . import ensemblgenomes, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("--compara", choices=ensemblgenomes.species_names())
    parser.add_argument("-C", "--clade", default="bep")
    args = parser.parse_args(argv or None)
    if args.config:
        read_config(args.config)
    if args.compara:
        with ensemblgenomes.FTPensemblgenomes() as ftp:
            dirs = ftp.download_maf(args.compara)
        pool = cli.ThreadPool()
        fts = [pool.submit(ensemblgenomes.consolidate_compara_mafs, d) for d in dirs]
        cli.wait_raise(fts)
        return
    tree = phylo.newicks[args.clade]
    species = phylo.extract_names(tree)
    if args.download:
        download(species)
    else:
        index(species)


def download(species: list[str]):
    assert species
    assert not (d := set(species) - set(ensemblgenomes.species_names_all())), d
    with ensemblgenomes.FTPensemblgenomes() as ftp:
        pool = cli.ThreadPool()
        futures: list[cli.FuturePath] = []
        for sp in species:
            fasta_dir = ftp.download_fasta(sp)
            futures.append(pool.submit(index_fasta, fasta_dir))
            gff3_dir = ftp.download_gff3(sp)
            futures.append(pool.submit(index_gff3, gff3_dir))
        cli.wait_raise(futures)
    if False:  # ensemble does not support rsync
        for sp in species:
            options = "--include *_sm.chromosome.*.fa.gz --exclude *.gz"
            ensemblgenomes.rsync(f"fasta/{sp}/dna", options)
            options = "--include *.chromosome.*.gff3.gz --exclude *.gz"
            ensemblgenomes.rsync(f"gff3/{sp}", options)


def index(species: list[str]):
    pool = cli.ThreadPool()
    futures: list[cli.FuturePath] = []
    for sp_dir in ensemblgenomes.species_dirs("fasta", species):
        futures.append(pool.submit(index_fasta, sp_dir))
    for sp_dir in ensemblgenomes.species_dirs("gff3", species):
        futures.append(pool.submit(index_gff3, sp_dir))
    if config["db"]["kmer"]:
        for sp_dir in ensemblgenomes.species_dirs("fasta", species):
            futures.append(pool.submit(softmask, sp_dir))
    cli.wait_raise(futures)


def softmask(species: Path):
    masked = jellyfish.run(species.name)
    return index_fasta(masked)


def index_fasta(path: Path):  # fasta/{species} or fasta/{species}/dna/kmer
    """Create bgzipped and indexed genome.fa."""
    if (path / "dna").exists():
        path /= "dna"
    for chromosome in fs.sorted_naturally(path.glob(r"*.chromosome.*.fa.gz")):
        kent.faToTwoBit(chromosome)
    genome = _create_genome_bgzip(path)
    htslib.faidx(genome)
    kent.faToTwoBit(genome)
    kent.faSize(genome)
    print(path)
    return genome


def index_gff3(path: Path):  # gff3/{species}
    """Create bgzipped and indexed genome.gff3."""
    genome = _create_genome_bgzip(path)
    htslib.tabix(genome)
    print(path)
    return genome


def _create_genome_bgzip(path: Path):
    """Combine chromosome files and bgzip it."""
    if (ext := path.parent.name) != "gff3":
        ext = "fa"
    files = fs.sorted_naturally(path.glob(rf"*.chromosome.*.{ext}.gz"))
    _log.debug(str(files))
    if cli.dry_run and not files:
        return Path("/dev/null")
    name = files[0].name
    (outname, count) = re.subn(rf"\.chromosome\..+\.{ext}", rf".genome.{ext}", name)
    assert count == 1
    outfile = path / outname
    return htslib.concat_bgzip(files, outfile)


if __name__ == "__main__":
    main()
