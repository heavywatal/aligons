import itertools
import logging
import os
from pathlib import Path

from .db import ensemblgenomes, phylo, stat
from .extern import kent, lastz, mafs2cram, multiz, phast
from .util import cli, read_config

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    available_species = ensemblgenomes.species_names()
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-N", "--check-args", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("-t", "--tips", type=int, default=0)
    parser.add_argument("-g", "--max-bp", type=float, default=float("inf"))
    parser.add_argument("--compara", action="store_true")
    parser.add_argument("target", choices=available_species)
    parser.add_argument("clade", choices=phylo.newicks.keys())
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.config:
        read_config(args.config)
    if args.check_args:
        return
    phastcons(args.target, args.clade, args.tips, args.max_bp, args.jobs, args.compara)


def phastcons(
    target: str, clade: str, tips: int, max_bp: float, jobs: int, compara: bool
):
    if compara:
        pairwise = Path("compara") / target
    else:
        pairwise = lastz.run(target, clade, jobs)
    mafs2cram.run(pairwise, clade, jobs)
    species = phylo.extract_tip_names(phylo.newicks[clade])
    species = list(filter(lambda x: test_fasize(x, max_bp), species))
    n = tips or len(species)
    for species in itertools.combinations(species, n):
        _log.info(f"{species = }")
        multiple = multiz.run(pairwise, species, jobs)
        phast.run(multiple, jobs)
        kent.run(multiple)


def test_fasize(species: str, max_bp: float):
    bp = stat.fasize(species)
    ret = bp < max_bp
    _log.info(f"{species:30}{round(bp / 1e6):>5} Mbp {ret}")
    return ret


if __name__ == "__main__":
    main()
