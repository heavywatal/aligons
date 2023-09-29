import concurrent.futures as confu
import contextlib
import gzip
import io
import itertools
import logging
import os
import re
import subprocess
from collections.abc import Iterable
from pathlib import Path
from typing import TypeVar
from zipfile import ZipFile

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
        if link.is_symlink():
            link.unlink()
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
    if source and source.exists() and product.stat().st_mtime < source.stat().st_mtime:
        return True
    return False


def newest(files: list[Path]):
    return max(files, key=lambda p: p.stat().st_mtime)


def sorted_naturally(iterable: Iterable[StrPath]) -> list[StrPath]:
    return sorted(iterable, key=natural_key)


def natural_key(x: StrPath) -> list[str]:
    return [try_zeropad(s) for s in re.split(r"[\W_]", name_if_path(x))]


def name_if_path(x: StrPath) -> str:
    return x if isinstance(x, str) else x.name


def try_zeropad(s: str):
    try:
        return f"{int(s):03}"
    except (ValueError, TypeError):
        return s


def zip_decompress(data: bytes) -> bytes:
    with ZipFile(io.BytesIO(data), "r") as zin:
        members = zin.namelist()
        assert len(members) == 1, members
        return zin.read(members[0])


def gzip_compress(content: bytes):
    if not is_gz(content):
        content = gzip.compress(content)
    return content


def gzip_decompress(content: bytes):
    if is_gz(content):
        content = gzip.decompress(content)
    return content


def is_gz(content: bytes):
    return content.startswith(b"\x1f\x8b")


def is_zip(content: bytes):
    return content.startswith(b"\x50\x4b")


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
        cmd = ["/usr/bin/sum", path]
        p = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True)
        try:
            (o_sum, o_blocks, _) = p.stdout.split()
        except ValueError:  # old coreutils
            (o_sum, o_blocks) = p.stdout.split()
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
