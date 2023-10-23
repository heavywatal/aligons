import itertools
import logging
from pathlib import Path

from .db import api, phylo
from .extern import kent, lastz, mafs2cram, multiz, phast
from .util import cli, log_config, read_config

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    tree = phylo.get_tree()
    parser = cli.ArgumentParser()
    parser.add_argument("-N", "--check-args", action="store_true")
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("-t", "--tips", type=int, default=0)
    parser.add_argument("-g", "--max-bp", type=float, default=float("inf"))
    parser.add_argument("--compara", action="store_true")
    parser.add_argument("target", choices=phylo.extract_tip_names(tree))
    parser.add_argument("clade", choices=phylo.extract_inner_names(tree))
    args = parser.parse_args(argv or None)
    if args.config:
        read_config(args.config)
    log_config()
    if args.check_args:
        return
    phastcons(args.target, args.clade, args.tips, args.max_bp, compara=args.compara)


def phastcons(
    target: str, clade: str, tips: int, max_bp: float, *, compara: bool
) -> None:
    lst_species = phylo.list_species(clade)
    lst_species = list(filter(lambda x: test_fasize(x, max_bp), lst_species))
    if compara:  # noqa: SIM108
        pairwise = Path("compara") / target
    else:
        pairwise = lastz.run(target, lst_species)
    fts = mafs2cram.run(pairwise, lst_species)
    n = tips or len(lst_species)
    for species in itertools.combinations(lst_species, n):
        if target not in species:
            continue
        _log.info(f"{species = }")
        multiple = multiz.run(pairwise, species)
        if n == len(lst_species):
            multiple.with_name(clade).symlink_to(multiple.name)
        phast.run(multiple)
        kent.run(multiple)
    cli.wait_raise(fts)


def test_fasize(species: str, max_bp: float) -> bool:
    bp = sum(x for x in api.chrom_sizes(species).values())
    ret = bp < max_bp
    _log.info(f"{species:30}{round(bp / 1e6):>5} Mbp {ret}")
    return ret


if __name__ == "__main__":
    main()
