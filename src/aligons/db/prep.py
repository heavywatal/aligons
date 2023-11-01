import logging
from pathlib import Path

from aligons.util import cli, read_config

from . import ensemblgenomes, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("--compara", choices=phylo.list_species())
    parser.add_argument("-C", "--clade", default="bep")
    args = parser.parse_args(argv or None)
    if args.config:
        read_config(args.config)
    if args.compara:
        ensemblgenomes.download_compara(args.compara)
        return
    species = phylo.list_species(args.clade)
    ensemblgenomes.download_via_ftp(species)


if __name__ == "__main__":
    main()
