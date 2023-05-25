"""http://plantdhs.org
"""
import logging
import re
import urllib.request
import zipfile
from pathlib import Path

from aligons import db
from aligons.extern import htslib
from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("pattern", nargs="?", default="*")
    args = parser.parse_args(argv or None)
    if args.download:
        for query in iter_download_queries():
            file = download(query)
            _log.info(f"{file}")
            unzipped = unzip_if_is_zipfile(file)
            if unzipped and unzipped.suffix == ".gff":
                # TODO: sort
                bgz = htslib.bgzip(unzipped)
                _log.info(f"{bgz}")
                tabix = htslib.tabix(bgz)
                _log.info(f"{tabix}")
    for x in fs.sorted_naturally(glob(args.pattern)):
        print(x)


def glob(pattern: str):
    return local_db_root().glob(pattern)


def download(query: str):
    host = "bioinfor.yzu.edu.cn"
    url = f"https://{host}/download/plantdhs/{query}"
    path = local_db_root() / query
    cmd = f"wget {url} -O {path}"
    subp.run_if(fs.is_outdated(path), cmd, shell=True)
    return path


def unzip_if_is_zipfile(path: Path):
    outfile = path.with_suffix("")
    if path.name.endswith(".gff.zip"):
        with zipfile.ZipFile(path, "r") as fin:
            members = fin.namelist()
            assert len(members) == 1
            assert not Path(members[0]).is_absolute()
            gzipped = outfile.with_suffix(outfile.suffix + ".gz")
            if fs.is_outdated(outfile, path) and fs.is_outdated(gzipped, path):
                fin.extract(members[0], path.parent)
        return outfile
    if path.name.endswith(".gff.gz"):
        tabix = path.with_suffix(path.suffix + ".csi")
        subp.run_if(fs.is_outdated(tabix, path), ["gunzip", path])
        return outfile
    return None


def iter_download_queries():
    for query in iter_download_queries_all():
        if re.search(r"Rice|TIGR7", query):
            yield query


def iter_download_queries_all():
    content = download_page()
    for mobj in re.finditer(r"/download/plantdhs/([^\"']+)", content):
        yield mobj[1]


def download_page():
    cache = local_db_root() / "download.html"
    if cache.exists():
        with cache.open("r") as fin:
            content = fin.read()
    else:
        local_db_root().mkdir(0o755, exist_ok=True)
        host = "plantdhs.org"
        url = f"http://{host}/Download"
        _log.info(url)
        response = urllib.request.urlopen(url)  # noqa: S310
        content = response.read().decode()
        with cache.open("w") as fout:
            fout.write(content)
    return content


def local_db_root():
    return db.path("plantdhs")


if __name__ == "__main__":
    main()
