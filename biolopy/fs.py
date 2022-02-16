import argparse
import logging
import re
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


def is_outdated(destination: Path, source: list[Path] | Path | None = None):
    if not destination.exists():
        return True
    if destination.stat().st_size == 0:
        return True
    if isinstance(source, list):
        source = newest(source)
    if source and destination.stat().st_ctime < source.stat().st_ctime:
        return True
    return False


def newest(files: list[Path]):
    champion = files[0]
    champion_ctime = 0
    for file in files:
        if (ctime := file.stat().st_ctime) > champion_ctime:
            champion = file
            champion_ctime = ctime
    return champion


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


if __name__ == "__main__":
    main()
