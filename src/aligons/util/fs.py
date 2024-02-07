import concurrent.futures as confu
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
        print(x)


def relpath(path: Path, start: Path = Path()) -> Path:
    assert start.is_dir() or not start.exists(), start
    return Path(os.path.relpath(path, start))


def symlink(path: Path, link: Path, *, relative: bool = False) -> Path:
    if is_outdated(link, path) and not cli.dry_run:
        if relative:
            start = link if path.is_dir() else link.parent
            path = relpath(path, start)
        _log.info(f"ln -s {path} {link}")
        if link.is_symlink():
            link.unlink()
        link.parent.mkdir(0o755, parents=True, exist_ok=True)
        link.symlink_to(path)
    _log.info(f"{link}")
    return link


def is_outdated(product: Path, source: Sequence[Path] | Path | None = None) -> bool:
    if not product.exists():
        return True
    if product.stat().st_size == 0:
        return True
    if isinstance(source, Sequence):
        source = [x for x in source if x.exists()]
        source = newest(source) if source else None
    if source and source.exists() and product.stat().st_mtime < source.stat().st_mtime:
        return True
    return False


def newest(files: Sequence[Path]) -> Path:
    return max(files, key=lambda p: p.stat().st_mtime)


def sorted_naturally[T: str | Path](iterable: Iterable[T]) -> list[T]:
    return sorted(iterable, key=natural_key)


def natural_key(x: str | Path) -> list[str]:
    return [try_zeropad(s) for s in re.split(r"[\W_]", name_if_path(x))]


def name_if_path(x: str | Path) -> str:
    return x if isinstance(x, str) else x.name


def try_zeropad(s: str) -> str:
    try:
        return f"{int(s):03}"
    except (ValueError, TypeError):
        return s


def expect_suffix(file: Path, suffix: str, *, negate: bool = False) -> None:
    if negate:
        if file.suffix == suffix:
            msg = f"{file}: unexpected suffix {suffix}"
            raise ValueError(msg)
    elif file.suffix != suffix:
        msg = f"{file}: expected suffix is {suffix}"
        raise ValueError(msg)


def checksums(file: Path) -> None:
    if cli.dry_run:  # pragma: no cover
        return
    with confu.ThreadPoolExecutor() as pool:
        with file.open("rt") as fin:
            lines = fin.readlines()
        pool.map(checkline, lines, itertools.repeat(file.parent))


def checkline(line: str, directory: Path) -> None:
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
def chdir(path: Path | str) -> Generator[None, Any, None]:
    """Change working directory temporarily with 'with' statement."""
    previous_wd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous_wd)


if __name__ == "__main__":
    main()
