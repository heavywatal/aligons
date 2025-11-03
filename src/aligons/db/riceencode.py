"""Rice Encyclopedia of DNA Elements.

<http://glab.hzau.edu.cn/RiceENCODE/>
"""

import logging
import re
from collections.abc import Iterable
from pathlib import Path

from aligons.util import cli, dl, fs

from . import _rsrc, api, tools

_log = logging.getLogger(__name__)
_HOST = "glab.hzau.edu.cn"


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
    slug = f"RiceENCODE/download/{query}"
    url = f"http://{_HOST}/{slug}"
    rawfile = _rsrc.db_root(_HOST) / slug
    outfile = db_prefix() / query
    if outfile.name.endswith(".bigWig"):
        outfile = outfile.with_suffix(".bw")
    elif outfile.name.endswith(".bed"):
        outfile = outfile.with_suffix(".bed.gz")
    response = dl.fetch(url, rawfile)
    return cli.thread_submit(tools.bgzip_or_symlink, response, outfile)


def iter_download_queries() -> Iterable[str]:
    content = download_page()
    content = re.sub(r"<!--.+?-->", "", content)
    content = re.sub(r"(/MH63/OpenChromatin/\w+).bw", r"\1.bigWig", content)
    for mobj in re.finditer(r"href=\"\.\./download/([^\"]+)", content):
        yield mobj[1]


def download_page() -> str:
    slug = "RiceENCODE/pages/download.html"
    url = f"http://{_HOST}/{slug}"
    cache = _rsrc.db_root(_HOST) / slug
    return dl.fetch(url, cache).text


def db_prefix() -> Path:
    return api.prefix("riceencode")


if __name__ == "__main__":
    main()
