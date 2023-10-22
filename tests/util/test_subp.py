import logging

import pytest
from aligons.util import ConfDict, cli, empty_options, subp

hello = "printf hello"


@pytest.mark.parametrize("cmd", [hello, hello.split()])
def test_run(cmd: str | list[str]):
    cli.dry_run = False
    assert subp.run(cmd, if_=False, stdout=subp.PIPE).stdout == b""
    assert subp.run(cmd, if_=True, stdout=subp.PIPE).stdout == b"hello"
    assert subp.run(cmd, stdout=subp.PIPE).stdout == b"hello"
    assert subp.run(cmd, shell=True, stdout=subp.PIPE).stdout == b"hello"  # noqa: S604


@pytest.mark.parametrize("cmd", [hello, hello.split()])
def test_popen(cmd: str | list[str]):
    cli.dry_run = False
    pf = subp.popen(cmd, if_=False, stdout=subp.PIPE)
    assert pf.communicate()[0] == b""
    pt = subp.popen(cmd, if_=True, stdout=subp.PIPE)
    assert pt.communicate()[0] == b"hello"
    assert subp.popen(cmd, stdout=subp.PIPE).communicate()[0] == b"hello"


def test_subp_log(caplog: pytest.LogCaptureFixture):
    cmd = hello
    subp.run(cmd, quiet=False)
    subp.run(cmd, quiet=True)
    assert not caplog.text

    caplog.set_level(logging.INFO)
    subp.run(cmd, quiet=False)
    assert cmd in caplog.text
    caplog.clear()
    subp.run(cmd, quiet=False, if_=False)
    assert f"# {cmd}" in caplog.text
    caplog.clear()
    subp.popen(cmd, quiet=False, if_=False).communicate()
    assert f"# {cmd}" in caplog.text
    caplog.clear()
    subp.run(cmd, quiet=True)
    assert not caplog.text

    caplog.clear()
    subp.popen(cmd, quiet=False).communicate()
    assert cmd in caplog.text
    caplog.clear()
    subp.popen(cmd, quiet=True).communicate()
    assert not caplog.text

    caplog.set_level(logging.DEBUG)
    subp.run(cmd, quiet=True)
    assert cmd in caplog.text
    caplog.clear()
    subp.popen(cmd, quiet=True).communicate()
    assert cmd in caplog.text


def test_optjoin():
    values = {
        "key": "value",
        "zero": 0,
        "one": 1,
        "true": True,
        "false": False,
        "none": None,
    }
    assert subp.optjoin(ConfDict(values)) == " --key=value --zero=0 --one=1 --true"
    assert subp.optjoin(ConfDict({"key": "value"}), "-") == " -key=value"
    assert not subp.optjoin(empty_options)
