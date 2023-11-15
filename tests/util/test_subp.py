import gzip
import logging
from pathlib import Path

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
    with subp.popen(cmd, if_=False, stdout=subp.PIPE) as p:
        assert p.communicate()[0] == b""
    with subp.popen(cmd, if_=True, stdout=subp.PIPE) as p:
        assert p.communicate()[0] == b"hello"
    with subp.popen(cmd, stdout=subp.PIPE) as p:
        assert p.communicate()[0] == b"hello"


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
    with subp.popen(cmd, quiet=False, if_=False) as p:
        p.communicate()
    assert f"# {cmd}" in caplog.text
    caplog.clear()
    subp.run(cmd, quiet=True)
    assert not caplog.text

    caplog.clear()
    with subp.popen(cmd, quiet=False) as p:
        p.communicate()
    assert cmd in caplog.text
    caplog.clear()
    with subp.popen(cmd, quiet=True) as p:
        p.communicate()
    assert not caplog.text

    caplog.set_level(logging.DEBUG)
    subp.run(cmd, quiet=True)
    assert cmd in caplog.text
    caplog.clear()
    with subp.popen(cmd, quiet=True) as p:
        p.communicate()
    assert cmd in caplog.text


def test_optargs():
    values = {
        "key": "value",
        "zero": 0,
        "one": 1,
        "true": True,
        "false": False,
        "none": None,
    }
    expected = ["--key=value", "--zero=0", "--one=1", "--true"]
    assert subp.optargs(ConfDict(values)) == expected
    assert subp.optargs(ConfDict({"key": "value"}), "-") == ["-key=value"]
    assert not subp.optargs(empty_options)


def test_open(tmp_path: Path):
    hello_txt = tmp_path / "hello"
    with subp.open_(hello_txt, "wt", if_=False) as fout:
        fout.write("hello")
    assert not hello_txt.exists()
    with subp.open_(hello_txt, "wt", if_=True) as fout:
        fout.write("hello")
    with hello_txt.open("rt") as fin:
        assert fin.read() == "hello"


def test_gzip(tmp_path: Path):
    hello_gz = tmp_path / "hello.gz"
    subp.gzip(b"hello", hello_gz, if_=False)
    assert not hello_gz.exists()
    subp.gzip(b"hello", hello_gz, if_=True)
    with subp.popen_zcat(hello_gz) as zcat:
        content, _ = zcat.communicate()
    assert content == b"hello"
    assert subp.run_zcat(hello_gz).stdout == b"hello"
    with gzip.open(hello_gz, "rb") as fin:
        assert fin.read() == content
    hello_gz.unlink()
    with subp.popen(hello, stdout=subp.PIPE) as phello:
        subp.gzip(phello.stdout, hello_gz)
    with gzip.open(hello_gz, "rt") as fin:
        assert fin.read().startswith("hello")


def test_sd():
    try:
        subp.run(["sd", "--version"], stdout=subp.subprocess.DEVNULL)
    except FileNotFoundError as e:
        pytest.skip(e.strerror)
    with (
        subp.popen(hello, stdout=subp.PIPE) as phello,
        subp.popen_sd("$", "ween", stdin=phello.stdout) as sd,
    ):
        result, _ = sd.communicate()
    assert result == b"helloween"
    with (
        subp.popen(hello, stdout=subp.PIPE) as phello,
        subp.popen_sd("", stdin=phello.stdout) as sd,
    ):
        result, _ = sd.communicate()
    assert result == b"hello"
