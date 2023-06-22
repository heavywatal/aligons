import csv
import logging

from aligons.util import cli

from . import api, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-l", "--long", action="store_true")
    parser.add_argument("-C", "--clade", default="angiospermae")
    args = parser.parse_args(argv or None)
    newick = phylo.get_subtree([args.clade])
    root = phylo.parse_newick(newick)
    for pre, species in phylo.rectangulate(phylo.render_tips(root, [])):
        if species not in api.species_names():
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
    path = api.fasize(species)
    with path.open() as fin:
        reader = csv.reader(fin, delimiter="\t")
        return list(reader)


def fasize(species: str):
    rows = chrom_sizes(species)
    lengths = [int(row[1]) for row in rows]
    return sum(lengths)


def gff3_size(species: str):
    path = api.genome_gff3(species)
    return path.stat().st_size / 1e6


if __name__ == "__main__":
    main()
