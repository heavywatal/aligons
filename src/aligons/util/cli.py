"""Command-line interface utilities."""

import argparse
import datetime
import logging
import sys
import threading
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from typing import TYPE_CHECKING, Any, override

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Sequence
    from pathlib import Path

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
    """Wrapper of `argparse.ArgumentParser` with common arguments."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        """Add common arguments for verbosity, dry-run, and concurrency."""
        super().__init__(*args, **kwargs)
        group = self.add_mutually_exclusive_group()
        group.add_argument("-v", "--verbose", action="count", default=1)
        group.add_argument("-q", "--quiet", action=DecreaseVerbosity, nargs=0)
        self.add_argument("-n", "--dry-run", action="store_true")
        self.add_argument("-j", "--jobs", type=int, default=1)

    @override
    def parse_args(  # type: ignore[reportIncompatibleMethodOverride]
        self,
        args: Sequence[str] | None = None,
        namespace: argparse.Namespace | None = None,
    ) -> argparse.Namespace:
        res = super().parse_args(args, namespace or argparse.Namespace())
        global dry_run  # noqa: PLW0603
        dry_run = res.dry_run
        level = _verbosity_to_level.get(res.verbose, logging.NOTSET)
        logging.basicConfig(
            level=level,
            handlers=[_stdout_handler(), _stderr_handler(), _file_handler()],
        )
        logging.logThreads = False
        logging.logProcesses = False
        logging.logMultiprocessing = False
        ThreadPool(res.jobs)
        return res


class DecreaseVerbosity(argparse.Action):
    """Custom action to decrease verbosity level."""

    @override
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str | Sequence[Any] | None,
        option_string: str | None = None,
    ) -> None:
        dest = "verbose"
        count = getattr(namespace, dest, 0)
        setattr(namespace, dest, max(count - 1, -3))


def _stdout_handler() -> logging.Handler:
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.DEBUG)
    handler.addFilter(lambda rec: rec.levelno < logging.WARNING)
    handler.setFormatter(logging.Formatter("%(msg)s"))
    return handler


def _stderr_handler() -> logging.Handler:
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.WARNING)
    return handler


def _file_handler() -> logging.Handler:
    today = _now().date()
    logfile = f".aligons.{today}.log"
    handler = logging.FileHandler(logfile, "a", delay=True)
    handler.setLevel(logging.WARNING)
    return handler


def _now() -> datetime.datetime:
    tz = datetime.datetime.now(datetime.UTC).astimezone().tzinfo
    return datetime.datetime.now(tz=tz)


class ThreadPool:
    """Singleton factory for ThreadPoolExecutor."""

    _instance = None

    def __new__(cls, max_workers: int | None = None) -> ThreadPoolExecutor:
        """Create a new ThreadPoolExecutor or return the existing one.

        :param max_workers: Passed to `ThreadPoolExecutor` constructor.
            The instance is replaced if it already exists and this value is not None.
        :returns: The singleton ThreadPoolExecutor instance.
        """
        if cls._instance is None:
            cls._instance = ThreadPoolExecutor(max_workers)
        elif max_workers is not None:
            max_w = cls._instance._max_workers  # noqa: SLF001
            _log.warning(f"Replaced ThreadPool with {max_workers=} (was {max_w})")
            cls._instance = ThreadPoolExecutor(max_workers)
        return cls._instance


def thread_submit(fn: Callable[..., Any], /, *args: Any, **kwargs: Any) -> Future[Any]:
    """Submit a callable to the singleton ThreadPoolExecutor.

    :param fn: Function to execute in a separate thread.
    :param args: Positional arguments to pass to the function.
    :param kwargs: Keyword arguments to pass to the function.
    :returns: Future representing the given call.
    """
    if threading.current_thread() != threading.main_thread():
        _log.warning("submit() from non-main thread may cause deadlock.")
    return ThreadPool().submit(fn, *args, **kwargs)


def wait_raise(futures: Iterable[Future[Any]]) -> None:
    """Wait for futures to complete and raise any exceptions."""
    for f in as_completed(futures):
        f.result()


def result(x: Path | Future[Path]) -> Path:
    """Syntactic sugar for getting result from a Future."""
    if isinstance(x, Future):
        return x.result()
    return x


def main(argv: list[str] | None = None) -> None:
    """Test this module."""
    parser = ArgumentParser()
    args = parser.parse_args(argv)
    level = _log.getEffectiveLevel()
    _log.info(f"{args = }")
    _log.info(f"{dry_run = }")
    _log.info(f"{level = }")
    _log.info(f"{logging.getLevelName(level) = }")
    _log.debug("debug message")
    _log.info("info message")
    _log.warning("warning message")
    _log.error("error message")
    _log.critical("critical message")


if __name__ == "__main__":
    main()
