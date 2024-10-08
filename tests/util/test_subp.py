import logging
import zipfile
from pathlib import Path

import pytest

from aligons.util import cli, subp

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
    assert subp.optargs(values) == expected
    assert subp.optargs({"key": "value"}, "-") == ["-key=value"]
    assert not subp.optargs({})


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
    content = b"hello"
    hello_gz = tmp_path / "hello.gz"
    subp.gzip(content, hello_gz, if_=False)
    assert not hello_gz.exists()
    subp.gzip(content, hello_gz, if_=True)
    assert subp.run_zcat(hello_gz).stdout == content
    hello_txt = hello_gz.with_suffix(".txt")
    assert subp.run_zcat(hello_gz, hello_txt).stdout is None
    with hello_txt.open("rb") as fin:
        assert fin.read() == content
    hello_gz.unlink()
    with subp.popen(hello, stdout=subp.PIPE) as p_hello:
        subp.gzip(p_hello.stdout, hello_gz)
    with subp.popen_zcat(hello_gz) as zcat:
        res, _ = zcat.communicate()
    assert res == content
    hello_zip = tmp_path / "hello.zip"
    with zipfile.ZipFile(hello_zip, "w") as fout:
        fout.write(hello_txt)
    hello_zip_txt = hello_zip.with_suffix(".zip.txt")
    assert subp.run_zcat(hello_zip, hello_zip_txt).stdout is None
    with hello_zip_txt.open("rb") as fin:
        assert fin.read() == content


def test_sd():
    with (
        subp.popen(hello, stdout=subp.PIPE) as p_hello,
        subp.popen_sd("", stdin=p_hello.stdout) as sd,
    ):
        result, _ = sd.communicate()
    assert result == b"hello"
    try:
        subp.run(["sd", "--version"], stdout=subp.subprocess.DEVNULL)
    except FileNotFoundError as e:
        pytest.skip(e.strerror)
    with (
        subp.popen(hello, stdout=subp.PIPE) as p_hello,
        subp.popen_sd("$", "ween", stdin=p_hello.stdout) as sd,
    ):
        result, _ = sd.communicate()
    assert result == b"helloween"
    with subp.popen(hello, stdout=subp.PIPE) as p_hello:
        result = subp.run_sd("", stdin=p_hello.stdout).stdout
    assert result == b"hello"
