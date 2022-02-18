import argparse
import logging
import os

from . import cli, lastz, multiz, phast, mafs2cram, kent
from .db import ensemblgenomes, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    available_species = ensemblgenomes.species_names()
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-N", "--check-args", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("target", choices=available_species)
    parser.add_argument("clade", choices=phylo.trees.keys())
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.check_args:
        return
    phastcons(args.target, args.clade, args.jobs)


def phastcons(target: str, clade: str, jobs: int):
    pairwise = lastz.run(target, clade, jobs)
    multiple = multiz.run(pairwise, clade, jobs)
    phast.run(multiple, jobs)
    kent.run(multiple)
    mafs2cram.run(pairwise, clade, jobs)


if __name__ == "__main__":
    main()
