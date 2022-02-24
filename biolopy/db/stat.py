import csv
import logging

from ..util import cli
from . import ensemblgenomes, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-c", "--clade", default="monocot")
    args = parser.parse_args(argv or None)
    cli.dry_run = args.dry_run
    cli.logging_config(args.loglevel)
    newick = phylo.trees[args.clade]
    root = phylo.parse_newick(newick)
    for pre, species in phylo.rectangulate(phylo.render_tips(root)):
        fasize, nseqs = chrom_sizes(species)
        gffsize = gff3_size(species)
        print(f"{pre} {species} {fasize:4.0f}Mbp {nseqs:2} {gffsize:4.1f}MB")


def chrom_sizes(species: str):
    path = ensemblgenomes.get_file("fasize.chrom.sizes", species)
    with open(path, "rt") as fin:
        reader = csv.reader(fin, delimiter="\t")
        lengths = [int(row[1]) for row in reader]
    return sum(lengths) / 1e6, len(lengths)


def gff3_size(species: str):
    path = ensemblgenomes.get_file("*.genome.gff3.gz", species)
    return path.stat().st_size / 1e6


if __name__ == "__main__":
    main()
