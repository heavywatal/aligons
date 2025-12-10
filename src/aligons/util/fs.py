"""File system utilities."""

import contextlib
import itertools
import logging
import os
import re
import subprocess
from collections.abc import Generator, Iterable, Sequence
from pathlib import Path
from typing import Any

from . import cli

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("path", nargs="+", type=Path)
    args = parser.parse_args(argv or None)
    for x in sorted_naturally(args.path):
        _log.info(x)


def print_if_exists(path: Path) -> Path:
    """Print the given path if it exists.

    :param path: The file or directory path to check.
    :returns: The same object passed as an argument.
    """
    if path.exists():
        _log.info(path)
    return path


def relpath(path: Path, start: Path = Path()) -> Path:
    """Get the relative path from `start` to `path`.

    :param path: Target file or directory.
    :param start: Starting directory.
    :returns: Relative path from `start` to `path`.
    """
    return Path(os.path.relpath(path, start).removeprefix("../"))


def symlink(path: Path, link: Path, *, relative: bool = False) -> Path:
    """Create a symbolic link after some preparation.

    :param path: The target file or directory to link to.
    :param link: The symbolic link to create.
    :param relative: If True, create a relative symlink.
    :returns: Path to the created symbolic link.
    """
    if is_outdated(link, path) and not cli.dry_run:
        if relative:
            path = relpath(path, link)
        _log.info(f"ln -s {path} {link}")
        if link.is_symlink():
            link.unlink()
        link.parent.mkdir(0o755, parents=True, exist_ok=True)
        link.symlink_to(path)
    print_if_exists(link)
    return link


def is_outdated(product: Path, source: Sequence[Path] | Path | None = None) -> bool:
    """Check if the product is outdated compared to the source.

    :param product: The product file or directory to check.
    :param source: The source files or directories to compare against.
    :returns: True if the product is outdated, empty, or does not exist.
    """
    if not product.exists():
        return True
    if product.stat().st_size == 0:
        return True
    if isinstance(source, Sequence):
        source = [x for x in source if x.exists()]
        source = newest(source) if source else None
    if source and source.exists():
        return product.stat().st_mtime < source.stat().st_mtime
    return False


def newest(files: Sequence[Path]) -> Path:
    """Find the most recently modified file.

    :param files: Files to check.
    :returns: The file with the maximum `st_mtime`.
    """
    return max(files, key=lambda p: p.stat().st_mtime)


def sorted_naturally[T: str | Path](iterable: Iterable[T]) -> list[T]:
    """Sort strings or Paths in numerical order.

    :param iterable: Strings or Paths.
    :returns: A list of sorted strings or Paths.
    """
    return sorted(iterable, key=_natural_key)


def _natural_key(x: str | Path) -> list[str]:
    return [_try_pad_zero(s) for s in re.split(r"[\W_]", _name_if_path(x))]


def _name_if_path(x: str | Path) -> str:
    return x if isinstance(x, str) else x.name


def _try_pad_zero(s: str) -> str:
    try:
        return f"{int(s):03}"
    except (ValueError, TypeError):
        return s


def expect_suffix(file: Path, suffix: str, *, negate: bool = False) -> None:
    """Check if the file ends with the expected suffix.

    :param file: File to check.
    :param suffix: Expected suffix.
    :param negate: If True, expect the file to NOT have the suffix.
    :raises ValueError: If `file` does not end with `suffix`.
    """
    if negate:
        if file.suffix == suffix:
            msg = f"{file}: unexpected suffix {suffix}"
            raise ValueError(msg)
    elif file.suffix != suffix:
        msg = f"{file}: expected suffix is {suffix}"
        raise ValueError(msg)


def checksums(file: Path) -> None:
    """Verify files using checksum file.

    :param file: Checksum file.
    """
    if cli.dry_run:  # pragma: no cover
        return
    with cli.ThreadPoolExecutor() as pool:
        with file.open("rt") as fin:
            lines = fin.readlines()
        pool.map(checkline, lines, itertools.repeat(file.parent))


def checkline(line: str, directory: Path) -> None:
    """Evaluate a single line of checksum file.

    :param line: A line from the checksum file.
    :param directory: Directory where the target files are located.
    :raises ValueError: If checksum does not match.
    """
    (e_sum, e_blocks, name) = line.split()
    if not (path := directory / name).exists():
        return
    cmd = ["/usr/bin/sum", path]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, check=True)  # noqa: S603
    try:
        (o_sum, o_blocks, _) = p.stdout.split()
    except ValueError:  # pragma: no cover; old coreutils
        (o_sum, o_blocks) = p.stdout.split()
    if (e_sum.lstrip("0"), e_blocks) != (o_sum.lstrip("0"), o_blocks):
        msg = f"sum {name}"
        msg += f"\nexpected: {e_sum}\t{e_blocks}"
        msg += f"\nobserved: {o_sum}\t{o_blocks}"
        _log.error(msg)


@contextlib.contextmanager
def chdir(path: Path | str) -> Generator[None, Any]:
    """Change working directory temporarily with 'with' statement."""
    previous_wd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous_wd)


if __name__ == "__main__":
    main()
