import logging
import re
from pathlib import Path

from aligons.db import tools
from aligons.extern import htslib, jellyfish, kent
from aligons.util import cli, config, fs, read_config

from . import ensemblgenomes, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("--compara", choices=phylo.list_species("angiospermae"))
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
    species = phylo.list_species(args.clade)
    prepare_ensemblgenomes(species)


def prepare_ensemblgenomes(species: list[str]):
    assert species
    species = ensemblgenomes.remove_unavailable(species)
    with ensemblgenomes.FTPensemblgenomes() as ftp:
        assert not (d := set(species) - set(ftp.available_species())), d
        pool = cli.ThreadPool()
        futures: list[cli.FuturePath] = []
        for sp in species:
            fasta_originals = ftp.download_fasta(sp)
            fasta_links = [symlink(x, sp, "fasta") for x in fasta_originals]
            futures.append(pool.submit(index_fasta, fasta_links))
            gff3_originals = ftp.download_gff3(sp)
            gff3_links = [symlink(x, sp, "gff3") for x in gff3_originals]
            futures.append(pool.submit(index_gff3, gff3_links))
        cli.wait_raise(futures)
    if config["db"]["kmer"]:
        futures.clear()
        for sp in species:
            futures.append(pool.submit(softmask, sp))
        cli.wait_raise(futures)
    if False:  # ensemble does not support rsync
        for sp in species:
            options = "--include *_sm.chromosome.*.fa.gz --exclude *.gz"
            ensemblgenomes.rsync(f"fasta/{sp}/dna", options)
            options = "--include *.chromosome.*.gff3.gz --exclude *.gz"
            ensemblgenomes.rsync(f"gff3/{sp}", options)


def symlink(path: Path, species: str, fmt: str = ""):
    if not fmt:
        fmt = "fasta" if path.name.removesuffix(".gz").endswith(".fa") else "gff3"
    assert fmt in ("fasta", "gff3"), fmt
    filename = path.name.replace("primary_assembly", "chromosome")
    link = ensemblgenomes.prefix() / fmt / species / filename
    return fs.symlink(path, link)


def index_fasta(paths: list[Path]):
    """Create bgzipped and indexed genome.fa."""
    if len(paths) == 1:
        assert "toplevel" in paths[0].name, paths[0]
        paths = split_toplevel_fa(paths[0])
    for chromosome in paths:
        kent.faToTwoBit(chromosome)
    genome = _create_genome_bgzip(paths)
    htslib.faidx(genome)
    kent.faToTwoBit(genome)
    kent.faSize(genome)
    return genome


def index_gff3(paths: list[Path]):  # gff3/{species}
    """Create bgzipped and indexed genome.gff3."""
    if len(paths) == 1:
        assert "toplevel" in paths[0].name, paths[0]
        paths = tools.split_gff(paths[0])
    genome = _create_genome_bgzip(paths)
    htslib.tabix(genome)
    return genome


def _create_genome_bgzip(files: list[Path]):
    """Combine chromosome files and bgzip it."""
    files = fs.sorted_naturally(files)
    _log.debug(str(files))
    if cli.dry_run and not files:
        return Path("/dev/null")
    name = files[0].name
    ext = files[0].with_suffix("").suffix
    (outname, count) = re.subn(rf"\.chromosome\..+{ext}", rf".genome{ext}", name)
    assert count == 1, name
    outfile = files[0].parent / outname
    return htslib.concat_bgzip(files, outfile)


def split_toplevel_fa(fa_gz: Path):
    fmt = "{stem}.{seqid}.fa.gz"
    fts = htslib.split_fa_gz(fa_gz, fmt, (r"toplevel", "chromosome"))
    return [f.result() for f in fts]


def softmask(species: str):
    masked = jellyfish.run(species)
    return index_fasta(masked)


if __name__ == "__main__":
    main()
