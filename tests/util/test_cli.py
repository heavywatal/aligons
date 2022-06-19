import logging
from biolopy.util import cli

_log = logging.getLogger(__name__)


def test_argparser():
    parser = cli.logging_argparser("vdq")
    assert parser.parse_args([]).loglevel == 0
    assert parser.parse_args(["-v"]).loglevel == 1
    assert parser.parse_args(["-vv"]).loglevel == 2
    assert parser.parse_args(["-d"]).loglevel == logging.DEBUG
    assert parser.parse_args(["-q"]).loglevel == logging.ERROR
    # pyright: reportPrivateUsage=false
    assert cli._from_verbosity(0) == logging.WARNING
    assert cli._from_verbosity(1) == logging.INFO
    assert cli._from_verbosity(2) == logging.DEBUG


def test_update_nested():
    x = {"both": 1, "x": 1}
    y = {"both": 2, "y": 2}
    cli.update_nested(x, y)
    assert x == {"both": 2, "x": 1, "y": 2}
    assert y == {"both": 2, "y": 2}
