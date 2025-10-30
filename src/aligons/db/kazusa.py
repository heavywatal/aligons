"""Kazusa Genome Atlas.

https://genome.kazusa.or.jp/
"""

import logging
from collections.abc import Iterator
from pathlib import Path

from aligons.util import cli, fs

from . import _rsrc, api, tools

_log = logging.getLogger(__name__)
_HOST = "genome.kazusa.or.jp"


def db_prefix() -> Path:
    return api.prefix("kazusa")


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-M", "--mask", action="store_true")
    args = parser.parse_args(argv or None)
    if args.download or args.mask:
        fts_fa: list[cli.Future[Path]] = []
        fts_gff: list[cli.Future[Path]] = []
        for ft_fa, ft_gff in _iter_fetch_and_bgzip():
            fts_fa.append(ft_fa)
            fts_gff.append(ft_gff)
        fts: list[cli.Future[Path]] = []
        for ft in cli.as_completed(fts_fa):
            if args.mask:
                masked = tools.softmask(ft.result())
                fts.extend(tools.genome_to_twobits(masked))
            else:
                fs.print_if_exists(ft.result())
        cli.wait_raise(fts)
        cli.wait_raise(fts_gff)


def _iter_fetch_and_bgzip() -> Iterator[tuple[cli.Future[Path], cli.Future[Path]]]:
    for entry in _rsrc.iter_builtin_dataset("kazusa.toml"):
        entry["url_prefix"] = f"https://{_HOST}"
        # .genome.fasta is bgzipped and faidx-ready
        # .gff3 is not sorted and not tabix-ready
        yield tools.fetch_and_bgzip(entry, db_prefix())


if __name__ == "__main__":
    main()
