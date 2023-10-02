import logging
import os
import urllib.request
from ftplib import FTP
from pathlib import Path
from urllib.parse import urlparse

import tomli_w

from aligons.util import cli, tomllib

_log = logging.getLogger(__name__)


def retrieve_content(
    url: str, outfile: Path | None = None, *, force: bool = False
) -> bytes:
    if outfile is None:
        urlp = urlparse(url)
        outfile = Path(urlp.netloc + urlp.path)
    _log.info(f"{outfile}")
    if cli.dry_run and not force:
        content = b""
    elif not outfile.exists():
        _log.info(url)
        with urllib.request.urlopen(url) as response:  # noqa: S310
            content = response.read()
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        with outfile.open("wb") as fout:
            fout.write(content)
    else:
        with outfile.open("rb") as fin:
            content = fin.read()
    return content


class LazyFTP(FTP):
    def __init__(self, host: str, slug: str, prefix: Path, timeout: float = 0):
        _log.info("LazyFTP()")
        self.host = host
        self.slug = slug
        self.prefix = prefix
        self.sizes: dict[str, int] = {}
        self.size_cache_toml = self.prefix / ".ftp_size_cache.toml"
        if self.size_cache_toml.exists():
            with self.size_cache_toml.open("rb") as fin:
                self.sizes = tomllib.load(fin)
        super().__init__(
            timeout=timeout or None  # pyright: ignore[reportGeneralTypeIssues]
        )

    def quit(self):  # noqa: A003
        if self.sizes:
            with self.size_cache_toml.open("wb") as fout:
                tomli_w.dump(self.sizes, fout)
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

    def size_cache(self, relpath: str):
        size = self.sizes.get(relpath, 0)
        if not size:
            self._lazy_init()
            size = self.size(relpath)
            assert size
            self.sizes[relpath] = size
        return size

    def retrieve(self, path: str, *, checksize: bool = False):
        outfile = self.prefix / path
        size_obs = outfile.stat().st_size if outfile.exists() else 0
        if checksize:
            size_exp = self.size_cache(path)
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
