"""HTTP and FTP download utilities with lazy access support."""

import logging
import tomllib
from ftplib import FTP
from pathlib import Path
from typing import TYPE_CHECKING, Any
from urllib.parse import urlparse

if TYPE_CHECKING:
    from collections.abc import Mapping

import requests
import tomli_w

from . import cli

_log = logging.getLogger(__name__)


class LazySession:
    """Wrapper of `requests.Session` that supports lazy initialization."""

    def __init__(self, url: str = "", data: Mapping[str, str] = {}) -> None:
        """Set properties without actually opening the session.

        :param url: URL to `POST` to initialize the session.
        :param data: Data to `POST` to initialize the session.
        """
        self._url = url
        self._data = data
        self._req_session = None

    def get(self, url: str, **kwargs: Any) -> requests.Response:
        """Send a GET request after initializing the session if needed.

        :param url: URL for the GET request.
        :param kwargs: Additional arguments passed to `requests.Session.get()`.
        :returns: Response object.
        """
        if self._req_session is None:
            self._req_session = requests.Session()
            if self._url:
                _log.info(f"POST {self._url}")
                response = self._req_session.post(self._url, data=self._data)
                response.raise_for_status()
        _log.info(f"GET {url}")
        return self._req_session.get(url, **kwargs)

    def fetch(self, url: str, outfile: Path | None = None) -> Response:
        """Create a lazy Response without actually downloading the resource.

        :param url: URL of the resource.
        :param outfile: Local path to save the file. `url` is used if `None`.
        :returns: A `Response` object.
        """
        if outfile is None:
            urlp = urlparse(url)
            outfile = Path(Path(urlp.path).name)
        return Response(self, url, outfile)

    def mirror(self, url: str, outdir: Path = Path()) -> Response:
        """Create a lazy Response that preserves the directory structure.

        :param url: URL of the resource.
        :param outdir: Output directory to save the file.
            Relative path from here corresponds to the URL path.
        :returns: A `Response` object.
        """
        urlp = urlparse(url)
        return self.fetch(url, outdir / (urlp.netloc + urlp.path))


class Response:
    """Wrapper of requests.Response that supports lazy access to the resource.

    HTTP session is opened when the file path is requested and not existing.
    The content is read from the file upon request.
    """

    def __init__(self, session: LazySession, url: str, path: Path) -> None:
        """Set properties without actually opening the session.

        :param session: LazySession to use for downloading.
        :param url: URL of the resource.
        :param path: Local path to save the file.
        """
        self._session = session
        self._url = url
        self._path = path
        self._content = b""

    @property
    def url(self) -> str:
        """URL of the resource."""
        return self._url

    @property
    def path(self) -> Path:
        """Path to the fetched file.

        Trigger fetching if the file does not exist and not in dry-run mode.
        """
        if not self._path.exists() and not cli.dry_run:
            self._fetch()
        return self._path

    @property
    def content(self) -> bytes:
        """Byte content of the fetched file. Empty in dry-run mode."""
        if cli.dry_run:
            return self._content
        return self.content_force

    @property
    def content_force(self) -> bytes:
        """Read the content from the file if available."""
        if not self._content and self.path.exists():
            _log.info(self.path)
            with self.path.open("rb") as fin:
                self._content = fin.read()
        return self._content

    @property
    def text(self) -> str:
        """String content of the fetched file. Empty in dry-run mode."""
        return self.content.decode()

    @property
    def text_force(self) -> str:
        """Read the string content from the file if available."""
        return self.content_force.decode()

    def _fetch(self) -> None:
        with self._session.get(self._url, stream=True) as response:
            response.raise_for_status()
            self._path.parent.mkdir(0o755, parents=True, exist_ok=True)
            iter_content = response.iter_content(chunk_size=2**28)
            chunk0 = next(iter_content)
            with self._path.open("wb") as fout:
                fout.write(chunk0)
                for chunk in iter_content:
                    fout.write(chunk)
            _log.info(self._path)


_global_session = LazySession()


def fetch(url: str, outfile: Path | None = None) -> Response:
    """Get a lazy Response from the global LazySession.

    :param url: URL of the resource.
    :param outfile: Local path to save the file. `url` is used if `None`.
    :returns: A `Response` object.
    """
    return _global_session.fetch(url, outfile)


def mirror(url: str, outdir: Path = Path()) -> Response:
    """Get a lazy Response that preserves the directory structure.

    :param url: URL of the resource.
    :param outdir: Output directory to save the file.
        Relative path from here corresponds to the URL path.
    :returns: A `Response` object.
    """
    return _global_session.mirror(url, outdir)


class LazyFTP(FTP):
    """Wrapper of `ftplib.FTP` that supports lazy initialization and caching."""

    def __init__(self, host: str, slug: str, prefix: Path, timeout: float = 0) -> None:
        """Set properties without connecting.

        :param host: URL of the FTP server.
        :param slug: Directory path on the FTP server.
        :param prefix: Local directory path to store downloaded files.
        :param timeout: Timeout for FTP connection in seconds.
        """
        _log.info("LazyFTP()")
        self.host = host
        self.slug = slug
        self.prefix = prefix
        self._size_cache: dict[str, int] = {}
        self._size_cache_toml = self.prefix / ".ftp_size_cache.toml"
        super().__init__(timeout=timeout or None)

    def quit(self) -> str:
        """Quit FTP session with logging."""
        _log.info("ftp.quit()")
        resp = super().quit()
        _log.info(resp)
        return resp

    def _lazy_init(self) -> None:
        if self.sock is not None:
            return
        _log.debug(f"ftp.connect({self.host})")
        _log.info(self.connect(self.host))
        _log.debug("ftp.login()")
        _log.info(self.login())
        _log.info(f"ftp.cwd({self.slug})")
        _log.info(self.cwd(self.slug))
        self.prefix.mkdir(0o755, parents=True, exist_ok=True)

    def nlst_cache(self, relpath: str) -> list[str]:
        """Call `FTP.nlst()` with caching.

        :param relpath: Relative directory path on the FTP server.
        :returns: List of file and directory names in the given directory.
        """
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

    def size(self, filename: str) -> int:
        """Get a file size with caching.

        :param filename: Relative file path on the FTP server.
        :returns: File size in bytes.
        """
        if not self._size_cache and self._size_cache_toml.exists():
            with self._size_cache_toml.open("rb") as fin:
                self._size_cache = tomllib.load(fin)
        size = self._size_cache.get(filename, 0)
        if not size:
            self._lazy_init()
            size = super().size(filename)
            if size is None:
                _log.warning(f"size not available for {filename}")
                size = -1
            self._size_cache[filename] = size
            with self._size_cache_toml.open("wb") as fout:
                tomli_w.dump(self._size_cache, fout)
        return size

    def retrieve(self, path: str, *, checksize: bool = False) -> Path:
        """Retrieve a file from the FTP server with caching.

        :param path: Relative file path on the FTP server.
        :param checksize: If True, re-download if the file size differs.
        :returns: Local path to the downloaded file.
        """
        outfile = self.prefix / path
        size_obs = outfile.stat().st_size if outfile.exists() else 0
        if checksize:
            size_exp = self.size(path)
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
            _log.warning(f"{outfile} ({size_obs=} != {size_exp=})")
        _log.info(outfile)
        return outfile
