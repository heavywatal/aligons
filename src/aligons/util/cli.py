import argparse
import concurrent.futures as confu
import logging
import os
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any, TypeAlias

FuturePath: TypeAlias = confu.Future[Path]

dry_run = False

_log = logging.getLogger(__name__)

_verbosity_to_level = {
    -2: logging.CRITICAL,
    -1: logging.ERROR,
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}


class ArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args: Any, **kwargs: Any):
        super().__init__(*args, **kwargs)
        group = self.add_mutually_exclusive_group()
        group.add_argument("-v", "--verbose", action=ConfigLogging, const=1)
        group.add_argument("-q", "--quiet", action=ConfigLogging, const=-1)
        self.add_argument("-n", "--dry-run", action="store_true")
        self.add_argument("-j", "--jobs", type=int, default=os.cpu_count())

    def parse_args(
        self,
        args: Sequence[str] | None = None,
        namespace: argparse.Namespace | None = None,
    ):
        res = super().parse_args(args, namespace)
        assert res is not None
        global dry_run  # noqa: PLW0603
        dry_run = res.dry_run
        verbosity = res.verbosity
        level = _verbosity_to_level.get(verbosity, logging.NOTSET)
        logging.basicConfig(level=level, handlers=[ConsoleHandler()])
        logging.logThreads = False
        logging.logProcesses = False
        logging.logMultiprocessing = False
        ThreadPool(res.jobs)
        return res


class ConfigLogging(argparse.Action):
    def __init__(
        self,
        option_strings: list[str],
        dest: str,
        nargs: int = 0,
        default: int = 0,
        **kwargs: Any,
    ):
        dest = "verbosity"
        super().__init__(option_strings, dest, nargs=nargs, default=default, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,  # noqa: ARG002
        namespace: argparse.Namespace,
        values: str | Sequence[Any] | None,
        option_string: str | None = None,
    ):
        assert not values, values
        assert option_string
        value = getattr(namespace, self.dest, 0)
        setattr(namespace, self.dest, max(value + self.const, -2))


class ConsoleHandler(logging.StreamHandler):  # type: ignore[reportMissingTypeArgument]
    def format(self, record: logging.LogRecord):  # noqa: A003
        if record.levelno < logging.WARNING:
            return record.msg
        return super().format(record)


class ThreadPool:
    _instance = None

    def __new__(cls, max_workers: int | None = None):
        if cls._instance is None:
            cls._instance = confu.ThreadPoolExecutor(max_workers)
        elif max_workers is not None:
            maxw = cls._instance._max_workers  # noqa: SLF001
            _log.warning(f"max_workers = {maxw}; ignored {max_workers}")
        return cls._instance


def thread_submit(fn: Callable[..., Any], /, *args: Any, **kwargs: Any):
    return ThreadPool().submit(fn, *args, **kwargs)


def wait_raise(futures: Iterable[confu.Future[Any]]):
    for f in confu.as_completed(futures):
        f.result()


def main(argv: list[str] | None = None):
    parser = ArgumentParser()
    args = parser.parse_args(argv)
    level = _log.getEffectiveLevel()
    print(f"{args = }")
    print(f"{dry_run = }")
    print(f"{level = }")
    print(f"{logging.getLevelName(level) = }")
    _log.debug("debug message")
    _log.info("info message")
    _log.warning("warning message")
    _log.error("error message")
    _log.critical("critical message")


if __name__ == "__main__":
    main()
