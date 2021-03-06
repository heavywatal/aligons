import concurrent.futures as confu
import itertools
import logging
import re
import subprocess
from collections.abc import Iterable
from pathlib import Path
from typing import TypeVar

from . import cli

StrPath = TypeVar("StrPath", str, Path)
_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.logging_argparser()
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
    return max(files, key=lambda p: p.stat().st_ctime)


def sorted_naturally(iterable: Iterable[StrPath]):
    return sorted(iterable, key=natural_key)


def natural_key(x: StrPath):
    return [try_zeropad(s) for s in re.split(r"[\W_]", name_if_path(x))]


def name_if_path(x: StrPath):
    return x if isinstance(x, str) else x.name


def try_zeropad(s: str):
    try:
        return f"{int(s):03}"
    except (ValueError, TypeError):
        return s


def checksums(file: Path):
    with confu.ThreadPoolExecutor() as pool:
        with file.open("rt") as fin:
            lines = fin.readlines()
        pool.map(checkline, lines, itertools.repeat(file.parent))


def checkline(line: str, directory: Path):
    (e_sum, e_blocks, name) = line.split()
    if (path := directory / name).exists():
        p = subprocess.run(["sum", path], stdout=subprocess.PIPE, text=True)
        (o_sum, o_blocks, _) = p.stdout.split()
        if (e_sum.lstrip("0"), e_blocks) != (o_sum, o_blocks):
            _log.error(f"{name}")
            _log.error(f"expected: {e_sum}\t{e_blocks}")
            _log.error(f"observed: {o_sum}\t{o_blocks}")


if __name__ == "__main__":
    main()
