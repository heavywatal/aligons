import logging
import os
from collections.abc import Mapping
from ftplib import FTP
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import requests
import tomli_w

from aligons.util import cli, fs, tomllib

_log = logging.getLogger(__name__)


class LazySession:
    def __init__(self, url: str = "", data: Mapping[str, str] = {}):
        self._url = url
        self._data = data
        self._req_session = None

    def _lazy_init(self):
        if self._req_session is None:
            self._req_session = requests.Session()
            if self._url:
                _log.info(f"POST {self._url}")
                response = self._req_session.post(self._url, data=self._data)
                response.raise_for_status()

    def get(self, url: str, **kwargs: Any):
        self._lazy_init()
        assert self._req_session
        _log.info(f"GET {url}")
        return self._req_session.get(url, **kwargs)

    def fetch(self, url: str, outfile: Path | None = None):
        if outfile is None:
            urlp = urlparse(url)
            outfile = Path(Path(urlp.path).name)
        return Response(self, url, outfile)

    def mirror(self, url: str, outdir: Path = Path()):
        urlp = urlparse(url)
        return self.fetch(url, outdir / (urlp.netloc + urlp.path))


class Response:
    def __init__(self, session: LazySession, url: str, path: Path):
        self._session = session
        self._url = url
        self._path = path
        self._content = b""

    def _fetch(self):
        response = self._session.get(self._url)
        response.raise_for_status()
        self._content = response.content
        if fs.is_gz(self._content) ^ (self._path.suffix == ".gz"):
            msg = f"gzip mismatch: '{self._url}' content vs filename '{self._path}'"
            raise ValueError(msg)
        _log.info(f"{self._path}")
        self._path.parent.mkdir(0o755, parents=True, exist_ok=True)
        with self._path.open("wb") as fout:
            fout.write(self._content)

    @property
    def url(self) -> str:
        return self._url

    @property
    def path(self) -> Path:
        if not self._path.exists() and not cli.dry_run:
            self._fetch()
        return self._path

    @property
    def content(self) -> bytes:
        if not self._content and not cli.dry_run:
            _log.info(f"{self.path}")
            with self.path.open("rb") as fin:
                self._content = fin.read()
        return self._content

    @property
    def text(self) -> str:
        return self.content.decode()


_global_session = LazySession()


def fetch(url: str, outfile: Path | None = None) -> Response:
    return _global_session.fetch(url, outfile)


def mirror(url: str, outdir: Path = Path()) -> Response:
    return _global_session.mirror(url, outdir)


class LazyFTP(FTP):
    def __init__(self, host: str, slug: str, prefix: Path, timeout: float = 0):
        _log.info("LazyFTP()")
        self.host = host
        self.slug = slug
        self.prefix = prefix
        self._size_cache: dict[str, int] = {}
        self._size_cache_toml = self.prefix / ".ftp_size_cache.toml"
        super().__init__(
            timeout=timeout or None  # pyright: ignore[reportGeneralTypeIssues]
        )

    def quit(self):  # noqa: A003
        _log.info(f"os.chdir({self.orig_wd})")
        os.chdir(self.orig_wd)
        _log.info("ftp.quit()")
        resp = super().quit()
        _log.info(resp)
        return resp

    def _lazy_init(self):
        if self.sock is not None:
            return
        self.orig_wd = Path.cwd()
        _log.debug(f"ftp.connect({self.host})")
        _log.info(self.connect(self.host))
        _log.debug("ftp.login()")
        _log.info(self.login())
        _log.info(f"ftp.cwd({self.slug})")
        _log.info(self.cwd(self.slug))
        _log.info(f"os.chdir({self.prefix})")
        self.prefix.mkdir(0o755, parents=True, exist_ok=True)
        os.chdir(self.prefix)  # for RETR only

    def nlst_cache(self, relpath: str):
        cache = self.prefix / relpath / ".ftp_nlst_cache"
        if cache.exists():
            _log.info(f"{cache=}")
            with cache.open("r") as fin:
                names = fin.read().rstrip().splitlines()
            lst = [str(Path(relpath) / x) for x in names]
        else:
            self._lazy_init()
            _log.info(f"ftp.nlst({relpath})")
            lst = self.nlst(relpath)  # mlsd is better but rarely supported
            cache.parent.mkdir(0o755, parents=True, exist_ok=True)
            with cache.open("w") as fout:
                fout.write("\n".join([Path(x).name for x in lst]) + "\n")
        return lst

    def size(self, filename: str):
        if not self._size_cache and self._size_cache_toml.exists():
            with self._size_cache_toml.open("rb") as fin:
                self._size_cache = tomllib.load(fin)
        size = self._size_cache.get(filename, 0)
        if not size:
            self._lazy_init()
            size = super().size(filename)
            assert size
            self._size_cache[filename] = size
            with self._size_cache_toml.open("wb") as fout:
                tomli_w.dump(self._size_cache, fout)
        return size

    def retrieve(self, path: str, *, checksize: bool = False):
        outfile = self.prefix / path
        size_obs = outfile.stat().st_size if outfile.exists() else 0
        if checksize:
            size_exp = self.size(path)
            assert size_exp
            is_to_run = size_obs != size_exp
        else:
            size_exp = 0
            is_to_run = size_obs == 0
        if is_to_run and not cli.dry_run:
            outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
            with outfile.open("wb") as fout:
                cmd = f"RETR {path}"
                self._lazy_init()
                _log.info(f"ftp.retrbinary({cmd}) {size_exp}")
                _log.info(self.retrbinary(cmd, fout.write))
            if size_exp and (size_obs := outfile.stat().st_size) != size_exp:
                _log.warning(f"{outfile} {size_obs} != {size_exp}")
        _log.info(f"{outfile}")
        return outfile
