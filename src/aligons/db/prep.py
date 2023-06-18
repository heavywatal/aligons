import logging
from pathlib import Path

from aligons.db import tools
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
            futures.append(pool.submit(tools.index_fasta, fasta_links))
            gff3_originals = ftp.download_gff3(sp)
            gff3_links = [symlink(x, sp, "gff3") for x in gff3_originals]
            futures.append(pool.submit(tools.index_gff3, gff3_links))
        cli.wait_raise(futures)
    if config["db"]["kmer"]:
        futures.clear()
        for sp in species:
            futures.append(pool.submit(tools.softmask, sp))
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


if __name__ == "__main__":
    main()
