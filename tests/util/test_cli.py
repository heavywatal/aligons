import logging
import time
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pytest

from aligons.util import cli


def test_arg_parser():
    parser = cli.ArgumentParser("vdq")
    assert parser.parse_args([]).verbose == 1
    assert parser.parse_args(["-v"]).verbose == 2  # noqa: PLR2004
    assert parser.parse_args(["-vv"]).verbose == 3  # noqa: PLR2004
    assert parser.parse_args(["-q"]).verbose == 0


def test_logging_config(caplog: pytest.LogCaptureFixture):
    cli.main([])
    assert "debug" not in caplog.text
    assert "info" not in caplog.text
    assert "warning" in caplog.text
    assert "error" in caplog.text
    assert "critical" in caplog.text
    caplog.clear()
    caplog.set_level(logging.INFO)
    cli.main([])
    assert "debug" not in caplog.text
    assert "info" in caplog.text
    caplog.clear()
    caplog.set_level(logging.DEBUG)
    cli.main([])
    assert "debug" in caplog.text


def test_threadpool(caplog: pytest.LogCaptureFixture):
    jobs = 4
    seconds = 0.02
    acceptance = 3
    cli.ThreadPool(jobs)
    start = time.time()
    fts = [cli.thread_submit(time.sleep, seconds) for _ in range(jobs)]
    cli.wait_raise(fts)
    elapsed = time.time() - start
    assert elapsed < seconds * acceptance
    cli.ThreadPool(42)
    assert "ignored" in caplog.text
