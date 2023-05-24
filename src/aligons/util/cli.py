import argparse
import logging
from typing import Any, Sequence

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

    def parse_args(
        self, args: list[str] | None = None, namespace: argparse.Namespace | None = None
    ):
        res = super().parse_args(args, namespace)
        assert res is not None
        global dry_run
        dry_run = res.dry_run
        verbosity = res.verbosity
        level = _verbosity_to_level.get(verbosity, logging.NOTSET)
        logging.basicConfig(level=level, handlers=[ConsoleHandler()])
        logging.logThreads = False
        logging.logProcesses = False
        logging.logMultiprocessing = False
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
        super().__init__(option_strings, "verbosity", nargs=0, default=0, **kwargs)

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str | Sequence[Any] | None,
        option_string: str | None = None,
    ):
        value = getattr(namespace, self.dest, 0)
        setattr(namespace, self.dest, max(value + self.const, -2))


class ConsoleHandler(logging.StreamHandler):  # type: ignore
    def format(self, record: logging.LogRecord):
        if record.levelno < logging.WARNING:
            return record.msg
        return super().format(record)


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
