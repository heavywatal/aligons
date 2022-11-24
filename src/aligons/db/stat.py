import csv
import logging

from ..util import cli
from . import ensemblgenomes, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-l", "--long", action="store_true")
    parser.add_argument("-c", "--clade", default="angiospermae")
    args = parser.parse_args(argv or None)
    cli.dry_run = args.dry_run
    cli.logging_config(args.loglevel)
    newick = phylo.newicks[args.clade]
    root = phylo.parse_newick(newick)
    for pre, species in phylo.rectangulate(phylo.render_tips(root)):
        if species not in ensemblgenomes.species_names():
            print(f"{pre} {species}")
            continue
        rows = chrom_sizes(species)
        lengths = [int(row[1]) for row in rows]
        fasize = sum(lengths) / 1e6
        nseqs = len(lengths)
        gffsize = gff3_size(species)
        print(f"{pre} {species} {fasize:4.0f}Mbp {nseqs:2} {gffsize:4.1f}MB")
        if args.long:
            print("\n".join([f"{row[0]:>6} {row[1]:>11}" for row in rows]))


def chrom_sizes(species: str):
    path = ensemblgenomes.get_file("fasize.chrom.sizes", species)
    with open(path, "rt") as fin:
        reader = csv.reader(fin, delimiter="\t")
        rows = list(reader)
    return rows


def fasize(species: str):
    rows = chrom_sizes(species)
    lengths = [int(row[1]) for row in rows]
    return sum(lengths)


def gff3_size(species: str):
    path = ensemblgenomes.get_file("*.genome.gff3.gz", species)
    return path.stat().st_size / 1e6


if __name__ == "__main__":
    main()
