import logging
import re
import shlex
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import IO, Any, Final, TypeAlias

from aligons.util import ConfDict

from . import cli

StrPath: TypeAlias = str | Path
Args: TypeAlias = list[StrPath]
_CMD: TypeAlias = Sequence[StrPath] | str
_FILE: TypeAlias = IO[Any] | int | None

CalledProcessError = subprocess.CalledProcessError
PIPE: Final = subprocess.PIPE
STDOUT: Final = subprocess.STDOUT

_log = logging.getLogger(__name__)


def popen(
    args: _CMD,
    *,
    if_: bool = True,
    stdin: _FILE = None,
    stdout: _FILE = None,
    quiet: bool = False,
):  # kwargs hinders type inference of output type [str | bytes]
    (args, cmd) = prepare_args(args, if_=if_)
    if quiet:
        _log.debug(cmd)
    else:
        _log.info(cmd)
    return subprocess.Popen(args, stdin=stdin, stdout=stdout)  # noqa: S603


def run(  # noqa: PLR0913
    args: _CMD,
    *,
    if_: bool = True,
    executable: StrPath | None = None,
    stdin: _FILE = None,
    input: bytes | None = None,  # noqa: A002
    stdout: _FILE = None,
    stderr: _FILE = None,
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
        args = ["sleep", "0"]
    return (args, cmd)


def optjoin(values: ConfDict, prefix: str = "--"):
    return "".join([optstr(k, v, prefix) for k, v in values.items()])


def optstr(key: str, value: Any, prefix: str = "--"):
    if value is None or value is False:
        return ""
    if value is True:
        return f" {prefix}{key}"
    return f" {prefix}{key}={value}"
