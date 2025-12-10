"""Wrappers for [`subprocess`](https://docs.python.org/3/library/subprocess.html)."""

import logging
import os
import re
import shlex
import subprocess
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from subprocess import PIPE, CompletedProcess, Popen
from typing import IO, Any, Final

from . import cli

type StrPath = str | Path
type Args = list[StrPath]
type _CMD = Sequence[StrPath] | str
type FILE = IO[Any] | int | None

CalledProcessError = subprocess.CalledProcessError
STDOUT: Final = subprocess.STDOUT

_log = logging.getLogger(__name__)


def popen(  # noqa: PLR0913
    args: _CMD,
    *,
    if_: bool = True,
    executable: StrPath | None = None,
    stdin: FILE = None,
    stdout: FILE = None,
    shell: bool = False,
    quiet: bool = False,
) -> Popen[bytes]:
    """Wrap [`subprocess.Popen`](https://docs.python.org/3/library/subprocess.html#subprocess.Popen).

    :param args: Command and arguments to execute.
    :param if_: Passed to `prepare_args()`.
    :param executable: Passed to `subprocess.Popen()`.
    :param stdin: Passed to `subprocess.Popen()`.
    :param stdout: Passed to `subprocess.Popen()`.
    :param shell: Passed to `subprocess.Popen()`.
    :param quiet: Decrease log level to `DEBUG` if `True`.
    :return: A `Popen` object.
    """
    # kwargs hinders type inference of output type [str | bytes]
    (args, cmd) = prepare_args(args, if_=if_)
    if quiet:
        _log.debug(cmd)
    else:
        _log.info(cmd)
    return Popen(  # noqa: S603
        args,
        executable=executable,
        stdin=stdin,
        stdout=stdout,
        shell=shell,
    )


def run(  # noqa: PLR0913
    args: _CMD,
    *,
    if_: bool = True,
    executable: StrPath | None = None,
    stdin: FILE = None,
    input: bytes | None = None,  # noqa: A002
    stdout: FILE = None,
    stderr: FILE = None,
    shell: bool = False,
    cwd: Path | None = None,
    text: bool | None = None,
    quiet: bool = False,
    check: bool = True,
) -> CompletedProcess[Any]:
    """Wrap [`subprocess.run`](https://docs.python.org/3/library/subprocess.html#subprocess.run).

    :param args: Command and arguments to execute.
    :param if_: Passed to `prepare_args()`.
    :param executable: Passed to `subprocess.run()`.
    :param stdin: Passed to `subprocess.run()`.
    :param input: Passed to `subprocess.run()`.
    :param stdout: Passed to `subprocess.run()`.
    :param stderr: Passed to `subprocess.run()`.
    :param shell: Passed to `subprocess.run()`.
    :param cwd: Passed to `subprocess.run()`.
    :param text: Passed to `subprocess.run()`.
    :param quiet: Decrease log level to `DEBUG` if `True`.
    :param check: Passed to `subprocess.run()`.
    :return: A `CompletedProcess` object.
    """
    # kwargs hinders type inference of output type [str | bytes]
    (args, cmd) = prepare_args(args, if_=if_)
    if quiet:
        _log.debug(cmd)
    else:
        _log.info(cmd)
    if shell:
        args = cmd
    return subprocess.run(  # noqa: S603
        args,
        executable=executable,
        stdin=stdin,
        input=input,
        stdout=stdout,
        stderr=stderr,
        shell=shell,
        cwd=cwd,
        text=text,
        check=check,
    )


def prepare_args(args: _CMD, *, if_: bool) -> tuple[Sequence[str | Path], str]:
    """Prepare argument list and concatenated command string.

    :param args: Command and arguments to execute.
    :param if_: Returns a no-op command if `False`.
    :return: The first element is the split arguments for execution;
             the second element is the concatenated command string for logging.
    """
    if isinstance(args, str):
        cmd = args.rstrip()
        args = shlex.split(cmd)
    else:
        cmd = " ".join(str(x) for x in args)
    if cli.dry_run or not if_:
        cmd = re.sub("^", "# ", cmd, flags=re.MULTILINE)
        args = ["true"]
    return (args, cmd)


def optargs(conf: Mapping[str, Any], prefix: str = "--") -> list[str]:
    """Convert a dictionary to a list of command-line options.

    :param conf: A dictionary of options.
    :param prefix: A prefix for each option.
    :return: A list of command-line options.
    """
    return list(_iter_optargs(conf, prefix))


def _iter_optargs(conf: Mapping[str, Any], prefix: str = "--") -> Iterator[str]:
    for key, value in conf.items():
        if value is None or value is False:
            continue
        if value is True:
            yield f"{prefix}{key}"
        else:
            yield f"{prefix}{key}={value}"


def open_(file: Path, mode: str, *, if_: bool = True) -> IO[Any]:
    """Conditionally open a file.

    :param file: A file path.
    :param mode: Passed to `Path.open()`.
    :param if_: Open `os.devnull` if `False`.
    :return: An opened file object.
    """
    if cli.dry_run or not if_:
        file = Path(os.devnull)
    return file.open(mode)


def gzip(data: FILE | bytes, outfile: Path, *, if_: bool = True) -> Path:
    """Compress data with `zstd --format=gzip`.

    :param data: A file-like object or bytes to compress.
    :param outfile: Output file path.
    :param if_: Do nothing if `False`.
    :return: The output file path.
    """
    assert outfile.suffix == ".gz", outfile
    args = ["zstd", "--format=gzip"]
    if_ = if_ and bool(data)
    with open_(outfile, "wb", if_=if_) as fout:
        if isinstance(data, bytes):
            run(args, input=data, stdout=fout, if_=if_, quiet=True)
        else:
            run(args, stdin=data, stdout=fout, if_=if_, quiet=True)
    return outfile


def popen_zcat(
    infile: Path | Sequence[Path], *, stdout: FILE = PIPE, if_: bool = True
) -> Popen[bytes]:
    """View compressed file(s) with `zstdcat` or `unzip`.

    :param infile: Input file path(s).
    :param stdout: Passed to `popen()`.
    :param if_: Do nothing if `False`.
    :return: A `Popen` object.
    """
    return popen(_zcat_args(infile), stdout=stdout, if_=if_, quiet=True)


def run_zcat(
    infile: Path | Sequence[Path],
    outfile: Path | None = None,
    *,
    stdout: FILE = PIPE,
    if_: bool = True,
) -> CompletedProcess[bytes]:
    """View compressed file(s) with `zstdcat` or `unzip`.

    :param infile: Input file path(s).
    :param outfile: Output file path.
    :param stdout: Passed to `run()`. Ignored if `outfile` is given.
    :param if_: Do nothing if `False`.
    :return: A `CompletedProcess` object.
    """
    args = _zcat_args(infile)
    if outfile is not None:
        if args[0] == "unzip":
            with open_(outfile, "wb", if_=if_) as fout:
                return run(args, stdout=fout, if_=if_, quiet=True)
        args = [*args, "-o", outfile]
        stdout = None
    return run(args, stdout=stdout, if_=if_, quiet=True)


def _zcat_args(infile: Path | Sequence[Path]) -> Args:
    if isinstance(infile, Path):
        infile = [infile]
    if infile[0].suffix == ".zip":
        return ["unzip", "-p", *infile]
    return ["zstdcat", "-q", *infile]


def popen_sd(
    pattern: str,
    repl: str = "",
    *,
    stdin: FILE = PIPE,
    stdout: FILE = PIPE,
    if_: bool = True,
) -> Popen[bytes]:
    """Open [`sd`](https://github.com/chmln/sd) process.

    :param pattern: Pattern to search.
    :param repl: Replacement string.
    :param stdin: Passed to `popen()`.
    :param stdout: Passed to `popen()`.
    :param if_: Passed to `popen()`
    :return: A `Popen` object.
    """
    if pattern:
        return popen(["sd", pattern, repl], stdin=stdin, stdout=stdout, if_=if_)
    assert not repl, repl
    return popen(["cat"], stdin=stdin, stdout=stdout, quiet=True, if_=if_)


def run_sd(
    pattern: str,
    repl: str = "",
    *,
    stdin: FILE = PIPE,
    stdout: FILE = PIPE,
    if_: bool = True,
) -> CompletedProcess[bytes]:
    """Run [`sd`](https://github.com/chmln/sd) command.

    :param pattern: Pattern to search.
    :param repl: Replacement string.
    :param stdin: Passed to `run()`.
    :param stdout: Passed to `run()`.
    :param if_: Passed to `run()`.
    :return: A `CompletedProcess` object.
    """
    return run(["sd", pattern, repl], stdin=stdin, stdout=stdout, if_=if_)
