import argparse
import concurrent.futures as confu
import itertools
import logging
import re
import subprocess
from collections.abc import Iterable
from pathlib import Path

from . import cli

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("path", nargs="+", type=Path)
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    for x in sorted_naturally(args.path):
        print(x)


def is_outdated(destination: Path, source: Path | None = None):
    if not destination.exists():
        return True
    if destination.stat().st_size == 0:
        return True
    if source and destination.stat().st_ctime < source.stat().st_ctime:
        return True
    return False


def sorted_naturally(iterable: Iterable[Path]):
    return sorted(iterable, key=natural_key)


def natural_key(x: Path):
    return [try_zeropad(s) for s in re.split(r"\W", x.name)]


def name_if_path(x: str | Path):
    return x if isinstance(x, str) else x.name


def try_zeropad(s: str):
    try:
        return f"{int(s):03}"
    except (ValueError, TypeError):
        return s


def checksums(files: Iterable[Path]):
    with confu.ThreadPoolExecutor(max_workers=8) as pool:
        for file in files:
            _log.info(f"{file}")
            with file.open("rt") as fin:
                lines = fin.readlines()
            pool.map(checkline, lines, itertools.repeat(file.parent))


def checkline(line: str, directory: Path):
    (e_sum, e_blocks, name) = line.split()
    if (path := directory / name).exists():
        p = cli.run(f"sum {path}", stdout=subprocess.PIPE, text=True, quiet=True)
        (o_sum, o_blocks, _) = p.stdout.split()
        if (e_sum.lstrip("0"), e_blocks) != (o_sum, o_blocks):
            _log.error(f"{name}")
            _log.error(f"expected: {e_sum}\t{e_blocks}")
            _log.error(f"observed: {o_sum}\t{o_blocks}")


if __name__ == "__main__":
    main()
