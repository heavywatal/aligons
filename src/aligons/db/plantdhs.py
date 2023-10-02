"""http://plantdhs.org."""
import logging
import re

from aligons import db
from aligons.extern import htslib
from aligons.util import cli, dl, fs

from . import tools

_log = logging.getLogger(__name__)
_HOST = "plantdhs.org"


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("pattern", nargs="?", default="*")
    args = parser.parse_args(argv or None)
    if args.download:
        cli.wait_raise(retrieve_deploy(q) for q in iter_download_queries())
    else:
        for x in fs.sorted_naturally(db_prefix().rglob(args.pattern)):
            print(x)


def retrieve_deploy(query: str):
    url = f"https://bioinfor.yzu.edu.cn/download/plantdhs/{query}"
    rawfile = db.path_mirror(_HOST) / query
    outfile = db_prefix() / query
    if outfile.name.endswith(".zip"):
        outfile = outfile.with_suffix(".gz")
    content = dl.get(url, rawfile).content
    if outfile.suffix == ".gz":
        future = cli.thread_submit(tools.compress, content, outfile)
    else:
        future = cli.thread_submit(fs.symlink, rawfile, outfile)
    return cli.thread_submit(htslib.try_index, future)


def iter_download_queries():
    for query in iter_download_queries_all():
        if re.search(r"Rice|TIGR7", query):
            yield query


def iter_download_queries_all():
    content = download_page()
    for mobj in re.finditer(r"/download/plantdhs/([^\"']+)", content):
        yield mobj[1]


def download_page() -> str:
    url = f"http://{_HOST}/Download"
    cache = db.path_mirror(_HOST) / "Download.html"
    return dl.get(url, cache).text


def db_prefix():
    return db.path("plantdhs")


if __name__ == "__main__":
    main()
