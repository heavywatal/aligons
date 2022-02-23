import argparse
import logging
import re
import shlex
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import IO, Any, Final, TypeAlias

StrPath: TypeAlias = str | Path[str]
Args: TypeAlias = list[StrPath]
_CMD: TypeAlias = Sequence[StrPath] | str
_FILE: TypeAlias = IO[Any] | int | None

CalledProcessError = subprocess.CalledProcessError
PIPE: Final = subprocess.PIPE
STDOUT: Final = subprocess.STDOUT
dry_run = False

_log = logging.getLogger(__name__)


def popen(
    args: _CMD,
    stdin: _FILE = None,
    stdout: _FILE = None,
    quiet: bool = False,
):
    return popen_if(True, args, stdin=stdin, stdout=stdout, quiet=quiet)


def run(
    args: _CMD,
    executable: StrPath | None = None,
    stdin: _FILE = None,
    stdout: _FILE = None,
    stderr: _FILE = None,
    shell: bool = False,
    cwd: Path | None = None,
    text: bool | None = None,
    quiet: bool = False,
):
    return run_if(
        True,
        args,
        executable=executable,
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
        shell=shell,
        cwd=cwd,
        text=text,
        quiet=quiet,
    )


def popen_if(
    cond: bool,
    args: _CMD,
    stdin: _FILE = None,
    stdout: _FILE = None,
    quiet: bool = False,
):  # kwargs hinders type inference of output type [str | bytes]
    (args, cmd) = prepare_args(args, cond)
    if quiet:
        _log.debug(cmd)
    else:
        _log.info(cmd)
    return subprocess.Popen(args, stdin=stdin, stdout=stdout)


def run_if(
    cond: bool,
    args: _CMD,
    executable: StrPath | None = None,
    stdin: _FILE = None,
    stdout: _FILE = None,
    stderr: _FILE = None,
    shell: bool = False,
    cwd: Path | None = None,
    text: bool | None = None,
    quiet: bool = False,
):  # kwargs hinders type inference of output type [str | bytes]
    (args, cmd) = prepare_args(args, cond)
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
        stdout=stdout,
        stderr=stderr,
        shell=shell,
        cwd=cwd,
        text=text,
        check=True,
    )


def prepare_args(args: _CMD, cond: bool = True):
    if isinstance(args, str):
        cmd = args.rstrip()
        args = shlex.split(args)
    else:
        cmd = " ".join(str(x) for x in args)
    if dry_run or not cond:
        cmd = re.sub("^", "# ", cmd, flags=re.MULTILINE)
        args = ["sleep", "0"]
    return (args, cmd)


def logging_argparser(options: str = "v") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_mutually_exclusive_group()
    if "v" in options:
        group.add_argument(
            "-v", "--verbose", action="count", default=0, dest="loglevel"
        )
    if "d" in options:
        group.add_argument(
            "-d", "--debug", action="store_const", const=logging.DEBUG, dest="loglevel"
        )
    if "q" in options:
        group.add_argument(
            "-q", "--quiet", action="store_const", const=logging.ERROR, dest="loglevel"
        )
    return parser


def logging_config(level: int | None = None):
    if level is not None and level < 10:
        level = _from_verbosity(level)
    logging.basicConfig(level=level, handlers=[ConsoleHandler()])
    logging.logThreads = False
    logging.logProcesses = False
    logging.logMultiprocessing = False


class ConsoleHandler(logging.StreamHandler):  # type: ignore
    def format(self, record: logging.LogRecord):
        if record.levelno < logging.WARNING:
            return record.msg
        return super().format(record)


def _from_verbosity(level: int):
    if level == 0:
        return logging.WARNING
    elif level == 1:
        return logging.INFO
    elif level == 2:
        return logging.DEBUG
    else:
        return logging.NOTSET


def main():
    parser = argparse.ArgumentParser(parents=[logging_argparser("vdq")])
    args = parser.parse_args()
    logging_config(args.loglevel)
    _log = logging.getLogger(__name__)
    _log.debug("debug message")
    _log.info("info message")
    _log.warning("warning message")
    _log.error("error message")
    _log.critical("critical message")


if __name__ == "__main__":
    main()
