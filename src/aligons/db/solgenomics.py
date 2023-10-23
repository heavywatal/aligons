"""https://solgenomics.net/ftp/genomes/."""
import logging
from pathlib import Path

from aligons import db
from aligons.util import cli

from . import phylo, tools

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-M", "--mask", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    args = parser.parse_args(argv or None)
    _test_newick()
    if args.download:
        fts = fetch_and_bgzip()
        cli.wait_raise(fts)
    if args.mask:
        fts = tools.index_as_completed(fetch_and_bgzip())
        cli.wait_raise(fts)


def fetch_and_bgzip() -> list[cli.FuturePath]:
    fts: list[cli.FuturePath] = []
    for entry in db.iter_builtin_dataset("solgenomics.toml"):
        entry["url_prefix"] = "https://solgenomics.net/ftp/genomes/"
        fts.extend(tools.fetch_and_bgzip(entry, db_prefix()))
    return fts


def db_prefix() -> Path:
    return db.path("solgenomics")


def _test_newick() -> None:
    for entry in db.iter_builtin_dataset("solgenomics.toml"):
        species = entry["species"]
        if species not in phylo.list_species():
            _log.warning(f"{species} not found in phylo")


if __name__ == "__main__":
    main()
