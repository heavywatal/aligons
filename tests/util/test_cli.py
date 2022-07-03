import logging

import pytest
from aligons.util import cli


def test_argparser():
    parser = cli.logging_argparser("vdq")
    assert parser.parse_args([]).loglevel == 0
    assert parser.parse_args(["-v"]).loglevel == 1
    assert parser.parse_args(["-vv"]).loglevel == 2
    assert parser.parse_args(["-d"]).loglevel == logging.DEBUG
    assert parser.parse_args(["-q"]).loglevel == logging.ERROR


def test_from_verbosity():
    # pyright: reportPrivateUsage=false
    assert cli._from_verbosity(0) == logging.WARNING
    assert cli._from_verbosity(1) == logging.INFO
    assert cli._from_verbosity(2) == logging.DEBUG
    assert cli._from_verbosity(42) == logging.NOTSET


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
