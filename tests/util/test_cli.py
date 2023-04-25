import logging

import pytest
from aligons.util import cli


def test_argparser():
    parser = cli.ArgumentParser("vdq")
    assert parser.parse_args([]).verbosity == 0
    assert parser.parse_args(["-v"]).verbosity == 1
    assert parser.parse_args(["-vv"]).verbosity == 2
    assert parser.parse_args(["-q"]).verbosity == -1


def test_from_verbosity():
    # pyright: reportPrivateUsage=false
    assert cli._verbosity_to_level.get(-1) == logging.ERROR
    assert cli._verbosity_to_level.get(0) == logging.WARNING
    assert cli._verbosity_to_level.get(1) == logging.INFO
    assert cli._verbosity_to_level.get(2) == logging.DEBUG
    assert cli._verbosity_to_level.get(42, logging.NOTSET) == logging.NOTSET


def test_logging_config(caplog: pytest.LogCaptureFixture):
    cli.main([])
    assert "debug" not in caplog.text
    assert "info" not in caplog.text
    assert "warning" in caplog.text
    assert "error" in caplog.text
    assert "critical" in caplog.text
    caplog.set_level(logging.INFO)
    cli.main([])
    assert "debug" not in caplog.text
    assert "info" in caplog.text
    caplog.set_level(logging.DEBUG)
    cli.main([])
    assert "debug" in caplog.text
