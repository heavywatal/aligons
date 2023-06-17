import concurrent.futures as confu
import contextlib
import itertools
import logging
import os
import re
import subprocess
from collections.abc import Iterable
from pathlib import Path
from typing import TypeVar

from . import cli

StrPath = TypeVar("StrPath", str, Path)
_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("path", nargs="+", type=Path)
    args = parser.parse_args(argv or None)
    for x in sorted_naturally(args.path):
        print(x)


def symlink(path: Path, link: Path):
    if is_outdated(link, path) and not cli.dry_run:
        _log.info(f"ln -s {path} {link}")
        link.parent.mkdir(0o755, parents=True, exist_ok=True)
        link.symlink_to(path)
    _log.info(f"{link}")
    return link


def is_outdated(product: Path, source: list[Path] | Path | None = None):
    if not product.exists():
        return True
    if product.stat().st_size == 0:
        return True
    if isinstance(source, list):
        source = [x for x in source if x.exists()]
        source = newest(source) if source else None
    if source and source.exists() and product.stat().st_ctime < source.stat().st_ctime:
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
    if cli.dry_run:
        return
    with confu.ThreadPoolExecutor() as pool:
        with file.open("rt") as fin:
            lines = fin.readlines()
        pool.map(checkline, lines, itertools.repeat(file.parent))


def checkline(line: str, directory: Path):
    (e_sum, e_blocks, name) = line.split()
    if (path := directory / name).exists():
        p = subprocess.run(["/usr/bin/sum", path], stdout=subprocess.PIPE, text=True)
        (o_sum, o_blocks, _) = p.stdout.split()
        if (e_sum.lstrip("0"), e_blocks) != (o_sum, o_blocks):
            _log.error(f"{name}")
            _log.error(f"expected: {e_sum}\t{e_blocks}")
            _log.error(f"observed: {o_sum}\t{o_blocks}")


@contextlib.contextmanager
def chdir(path: Path | str):
    """Change working directory temporarily with 'with' statement."""
    previous_wd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous_wd)


if __name__ == "__main__":
    main()
