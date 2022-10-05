"""http://plantdhs.org
"""
import logging
import re
import urllib.request
import zipfile
from pathlib import Path

from ..extern import htslib
from ..util import cli, config, fs, subp

_log = logging.getLogger(__name__)
LOCAL_DB_ROOT: Path = config["db"]["root"] / "plantdhs"


def main(argv: list[str] | None = None):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("pattern", nargs="?", default="*")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.download:
        for query in iter_download_queries():
            file = download(query)
            _log.info(f"{file}")
            if unzipped := unzip_if_is_zipfile(file):
                if unzipped.suffix == ".gff":
                    # TODO: sort
                    bgz = htslib.bgzip(unzipped)
                    _log.info(f"{bgz}")
                    tabix = htslib.tabix(bgz)
                    _log.info(f"{tabix}")
    for x in fs.sorted_naturally(glob(args.pattern)):
        print(x)


def glob(pattern: str):
    return LOCAL_DB_ROOT.glob(pattern)


def download(query: str):
    host = "bioinfor.yzu.edu.cn"
    url = f"https://{host}/download/plantdhs/{query}"
    path = LOCAL_DB_ROOT / query
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
    elif path.name.endswith(".gff.gz"):
        tabix = path.with_suffix(path.suffix + ".csi")
        subp.run_if(fs.is_outdated(tabix, path), ["gunzip", path])
        return outfile


def iter_download_queries():
    for query in iter_download_queries_all():
        if re.search(r"Rice|TIGR7", query):
            yield query


def iter_download_queries_all():
    content = download_page()
    for mobj in re.finditer(r"/download/plantdhs/([^\"']+)", content):
        yield mobj[1]


def download_page():
    cache = LOCAL_DB_ROOT / "download.html"
    if cache.exists():
        with cache.open("r") as fin:
            content = fin.read()
    else:
        LOCAL_DB_ROOT.mkdir(0o755, exist_ok=True)
        host = "plantdhs.org"
        url = f"http://{host}/Download"
        _log.info(url)
        response = urllib.request.urlopen(url)
        content = response.read().decode()
        with cache.open("w") as fout:
            fout.write(content)
    return content


if __name__ == "__main__":
    main()
