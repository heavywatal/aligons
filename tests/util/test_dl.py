import email
import io
import logging
import urllib.request
import urllib.response
from pathlib import Path

import pytest
from aligons.util import dl

url_netloc = "localhost"
url_path = "/path/to/file.txt"
url_full = f"http://{url_netloc}{url_path}"


@pytest.fixture(scope="module")
def tmp_path_module(tmp_path_factory: pytest.TempPathFactory):
    return tmp_path_factory.mktemp("dl")


@pytest.fixture()
def monkey_url(tmp_path_module: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.chdir(tmp_path_module)
    monkeypatch.setattr(urllib.request, "urlopen", mock_urlopen)
    return url_full


def mock_urlopen(url: str):
    fp = io.BytesIO(b"content")
    msg = email.message_from_string("Content-Type: text/plain\n")
    return urllib.response.addinfourl(fp, msg, url, 200)


def _test_response(res: dl.Response, path: Path):
    assert res.url == url_full
    assert res.path == path
    assert res.path.is_file()
    assert res.content == b"content"
    assert res.text == "content"


def test_get(monkey_url: str, caplog: pytest.LogCaptureFixture):
    caplog.set_level(logging.INFO)
    outfile = Path(Path(url_path).name)
    _test_response(dl.get(monkey_url), outfile)
    assert url_full in caplog.text
    caplog.clear()
    _test_response(dl.get(monkey_url), outfile)  # cached
    assert url_full not in caplog.text
    assert caplog.text.count(str(outfile)) == 1


def test_get_outfile(monkey_url: str):
    outfile = Path("outfile.txt")
    _test_response(dl.get(monkey_url, outfile), outfile)


def test_mirror(monkey_url: str):
    outfile = Path(f"{url_netloc}{url_path}")
    _test_response(dl.mirror(monkey_url), outfile)


def test_mirror_outdir(monkey_url: str):
    outdir = Path("outdir")
    outfile = outdir / f"{url_netloc}{url_path}"
    _test_response(dl.mirror(monkey_url, outdir), outfile)
