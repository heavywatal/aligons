"""https://solgenomics.net/ftp/genomes/."""
import logging
from collections.abc import Iterator
from pathlib import Path

from aligons.util import cli

from . import _rsrc, api, phylo, tools

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-M", "--mask", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    args = parser.parse_args(argv or None)
    _test_newick()
    if args.download:
        fts: list[cli.Future[Path]] = []
        for pair in iter_fetch_and_bgzip():
            fts.extend(pair)
        cli.wait_raise(fts)
    if args.mask:
        pairs = list(iter_fetch_and_bgzip())
        fts = [tools.process_genome(x) for x in pairs]
        cli.wait_raise(fts)


def iter_fetch_and_bgzip() -> Iterator[tuple[cli.Future[Path], cli.Future[Path]]]:
    for entry in _rsrc.iter_builtin_dataset("solgenomics.toml"):
        entry["url_prefix"] = "https://solgenomics.net/ftp/genomes/"
        yield tools.fetch_and_bgzip(entry, db_prefix())


def db_prefix() -> Path:
    return api.prefix("solgenomics")


def _test_newick() -> None:
    for entry in _rsrc.iter_builtin_dataset("solgenomics.toml"):
        species = entry["species"]
        if species not in phylo.list_species():
            _log.warning(f"{species} not found in phylo")


if __name__ == "__main__":
    main()
