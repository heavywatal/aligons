import argparse
import csv
import logging

from .. import cli
from . import ensemblgenomes

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    args = parser.parse_args(argv or None)
    cli.dry_run = args.dry_run
    cli.logging_config(args.loglevel)
    genome_size()
    genome_gff3_size()


def genome_size():
    for x in sorted(ensemblgenomes.rglob("fasize.chrom.sizes")):
        species = x.parent.parent.name
        with open(x, "rt") as fin:
            reader = csv.reader(fin, delimiter="\t")
            lengths = [int(row[1]) for row in reader]
        num_seqs = len(lengths)
        total_length = round(sum(lengths) / 1e6)
        print(f"{species:24} {total_length:5,}Mbp {num_seqs:2}")


def genome_gff3_size():
    for x in sorted(ensemblgenomes.rglob("*.genome.gff3.gz")):
        species = x.parent.name
        size = round(x.stat().st_size / 1e3)
        print(f"{species:24} {size:6,}KB")


if __name__ == "__main__":
    main()
