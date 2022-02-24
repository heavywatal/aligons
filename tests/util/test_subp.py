import pytest
from biolopy.util import subp

hello = "printf hello"


@pytest.mark.parametrize("cmd", [hello, hello.split()])
def test_run(cmd: str | list[str]):
    assert subp.run_if(False, cmd, stdout=subp.PIPE).stdout == b""
    assert subp.run_if(True, cmd, stdout=subp.PIPE).stdout == b"hello"
    assert subp.run(cmd, stdout=subp.PIPE).stdout == b"hello"


@pytest.mark.parametrize("cmd", [hello, hello.split()])
def test_popen(cmd: str | list[str]):
    assert subp.popen_if(False, cmd, stdout=subp.PIPE).communicate()[0] == b""
    assert subp.popen_if(True, cmd, stdout=subp.PIPE).communicate()[0] == b"hello"
    assert subp.popen(cmd, stdout=subp.PIPE).communicate()[0] == b"hello"
