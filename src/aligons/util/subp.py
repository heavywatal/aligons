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
):  # kwargs hinders type inference of output type [str | bytes]
    (args, cmd) = prepare_args(args, if_=if_)
    if quiet:
        _log.debug(cmd)
    else:
        _log.info(cmd)
    return Popen(
        args,
        executable=executable,
        stdin=stdin,
        stdout=stdout,
        shell=shell,  # noqa: S603
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
):  # kwargs hinders type inference of output type [str | bytes]
    (args, cmd) = prepare_args(args, if_=if_)
    if quiet:
        _log.debug(cmd)
    else:
        _log.info(cmd)
    if shell:
        args = cmd
    return subprocess.run(
        args,
        executable=executable,
        stdin=stdin,
        input=input,
        stdout=stdout,
        stderr=stderr,
        shell=shell,  # noqa: S603
        cwd=cwd,
        text=text,
        check=check,
    )


def prepare_args(args: _CMD, *, if_: bool):
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
    if cli.dry_run or not if_:
        file = Path(os.devnull)
    return file.open(mode)


def gzip(data: FILE | bytes, outfile: Path, *, if_: bool = True) -> Path:
    assert outfile.suffix == ".gz", outfile
    args = ["zstd", "--format=gzip", "-T2"]
    if_ = if_ and bool(data)
    with open_(outfile, "wb", if_=if_) as fout:
        if isinstance(data, bytes):
            run(args, input=data, stdout=fout, if_=if_, quiet=True)
        else:
            run(args, stdin=data, stdout=fout, if_=if_, quiet=True)
    return outfile


def popen_zcat(infile: Path, stdout: FILE = PIPE, *, if_: bool = True) -> Popen[bytes]:
    return popen(_zcat_args(infile), stdout=stdout, if_=if_, quiet=True)


def run_zcat(
    infile: Path, stdout: FILE = PIPE, *, if_: bool = True
) -> CompletedProcess[bytes]:
    return run(_zcat_args(infile), stdout=stdout, if_=if_, quiet=True)


def _zcat_args(infile: Path) -> Args:
    if infile.suffix == ".zip":
        return ["unzip", "-p", infile]
    return ["zstdcat", "-T2", infile]


def popen_sd(
    pattern: str, repl: str = "", *, stdin: FILE = PIPE, stdout: FILE = PIPE
) -> Popen[bytes]:
    if pattern:
        return popen(["sd", pattern, repl], stdin=stdin, stdout=stdout)
    return popen(["cat"], stdin=stdin, stdout=PIPE, quiet=True)
