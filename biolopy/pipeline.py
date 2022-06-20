import logging
import os
from pathlib import Path

from .db import ensemblgenomes, phylo
from .extern import kent, lastz, mafs2cram, multiz, phast
from .util import cli

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    available_species = ensemblgenomes.species_names()
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-N", "--check-args", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("--compara", action="store_true")
    parser.add_argument("target", choices=available_species)
    parser.add_argument("clade", choices=phylo.newicks.keys())
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.config:
        cli.read_config(args.config)
    if args.check_args:
        return
    phastcons(args.target, args.clade, args.jobs, args.compara)


def phastcons(target: str, clade: str, jobs: int, compara: bool):
    if compara:
        pairwise = Path("compara") / target
    else:
        pairwise = lastz.run(target, clade, jobs)
    mafs2cram.run(pairwise, clade, jobs)
    multiple = multiz.run(pairwise, clade, jobs)
    phast.run(multiple, jobs)
    kent.run(multiple)


if __name__ == "__main__":
    main()
