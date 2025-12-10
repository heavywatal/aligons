"""PlantDHS: a database for DNase I hypersensitive sites in plants.

<http://plantdhs.org>.
"""

import logging
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path

from aligons.util import cli, dl, fs

from . import _rsrc, api, tools

_log = logging.getLogger(__name__)
_HOST = "plantdhs.org"


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("pattern", nargs="?", default="*")
    args = parser.parse_args(argv or None)
    if args.download:
        cli.wait_raise(retrieve_deploy(q) for q in iter_download_queries())
    else:
        for x in fs.sorted_naturally(db_prefix().rglob(args.pattern)):
            fs.print_if_exists(x)


def retrieve_deploy(query: str) -> cli.Future[Path]:
    url = f"https://bioinfor.yzu.edu.cn/download/plantdhs/{query}"
    rawfile = _rsrc.db_root(_HOST) / query
    outfile = db_prefix() / query
    if outfile.name.endswith(".zip"):
        outfile = outfile.with_suffix(".gz")
    response = dl.fetch(url, rawfile)
    return cli.thread_submit(tools.bgzip_or_symlink, response, outfile)


def iter_download_queries() -> Iterable[str]:
    for query in iter_download_queries_all():
        if re.search(r"Rice|TIGR7", query):
            yield query


def iter_download_queries_all() -> Iterable[str]:
    content = download_page()
    for mobj in re.finditer(r"/download/plantdhs/([^\"']+)", content):
        yield mobj[1]


def download_page() -> str:
    url = f"http://{_HOST}/Download"
    cache = _rsrc.db_root(_HOST) / "Download.html"
    return dl.fetch(url, cache).text


def db_prefix() -> Path:
    return api.prefix("plantdhs")


if __name__ == "__main__":
    main()
