"""Check output files to remove incomplete results."""

import logging
import shutil
from pathlib import Path
from typing import TYPE_CHECKING

from . import cli, fs

if TYPE_CHECKING:
    from collections.abc import Iterable

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    """CLI for checking incomplete results."""
    parser = cli.ArgumentParser()
    parser.add_argument("--delete", action="store_true")
    parser.add_argument("directory", type=Path)
    args = parser.parse_args(argv)
    if "pairwise" in args.directory.name:
        candidates = list(check_pairwise(args.directory))
    elif "multiple" in args.directory.name:
        candidates = list(check_multiple(args.directory))
    else:
        _log.error(f"Unknown directory: {args.directory}")
        return
    for directory in candidates:
        _log.info(directory)
        if args.delete:
            shutil.rmtree(directory)


def check_pairwise(directory: Path) -> Iterable[Path]:
    """Check pairwise alignment output files for missing or small files.

    :param directory: Directory to check.
    """
    return _check_file_size(directory, "sing.maf", 256 * 1024)


def check_multiple(directory: Path) -> Iterable[Path]:
    """Check multiple alignment output files for missing or small files.

    :param directory: Directory to check.
    """
    return _check_file_size(directory, "multiz.maf", 1024 * 1024)


def _check_file_size(directory: Path, filename: str, threshold: int) -> Iterable[Path]:
    """Check files for missing or small files.

    :param directory: Directory to check.
    :param filename: Filename to check for.
    :param threshold: Size threshold in bytes.
    """
    for dirpath, _subdirs, _files in directory.walk():
        if not dirpath.name.startswith("chromosome"):
            continue
        file_path = dirpath / filename
        if not file_path.exists():
            _log.warning(f"{file_path} missing")
            yield print_ls_long(dirpath)
        elif (size := file_path.stat().st_size) < threshold:
            _log.warning(f"{file_path} ({size} bytes)")
            yield print_ls_long(dirpath)


def print_ls_long(directory: Path) -> Path:
    """Print ls -l output for a directory.

    :param directory: Directory to list.
    """
    for x in fs.sorted_naturally(directory.iterdir()):
        size = x.stat().st_size if x.exists() else -1
        _log.info(f"{size:>10} {x.name}")
    return directory


if __name__ == "__main__":
    main()
