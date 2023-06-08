"""http://plantdhs.org."""
import logging
import re

from aligons import db
from aligons.util import cli, fs

from . import tools

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
    for x in fs.sorted_naturally(glob(args.pattern)):
        print(x)


def glob(pattern: str):
    return local_db_root().glob(pattern)


def download(query: str):
    url = f"https://bioinfor.yzu.edu.cn/download/plantdhs/{query}"
    outfile = local_db_root() / query
    if outfile.name.endswith(".gff.zip"):
        outfile = outfile.with_suffix(".gz")
    return tools.retrieve_compress(url, outfile)


def iter_download_queries():
    for query in iter_download_queries_all():
        if re.search(r"Rice|TIGR7", query):
            yield query


def iter_download_queries_all():
    content = download_page()
    for mobj in re.finditer(r"/download/plantdhs/([^\"']+)", content):
        yield mobj[1]


def download_page():
    url = "http://plantdhs.org/Download"
    cache = local_db_root() / "download.html"
    return tools.retrieve_cache(url, cache).decode()


def local_db_root():
    return db.path("plantdhs")


if __name__ == "__main__":
    main()
