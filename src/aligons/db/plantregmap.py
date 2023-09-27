"""http://plantregmap.gao-lab.org/."""
import logging
import re
from collections.abc import Iterator
from pathlib import Path

from aligons import db
from aligons.extern import htslib
from aligons.util import cli, fs

from . import ftplib, tools

_log = logging.getLogger(__name__)
_HOST = "plantregmap.gao-lab.org"

_longer = {
    "Osj": "Oryza_sativa_Japonica_Group",
    "Sly": "Solanum_lycopersicum",
}


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("pattern", nargs="?", default="*")
    args = parser.parse_args(argv or None)
    if args.download:
        with FTPplantregmap() as ftp:
            sp = "Osj"
            ftp.ls_cache(_longer[sp])
            ftp.download_pairwise_alignments(sp)
            ftp.download_multiple_alignments(_longer[sp])
            ftp.download_conservation(_longer[sp])
        cli.wait_raise(retrieve_deploy(q) for q in iter_download_queries())
    else:
        for x in fs.sorted_naturally(db_prefix().rglob(args.pattern)):
            print(x)


def retrieve_deploy(query: str):
    url = f"http://{_HOST}/download_ftp.php?{query}"
    relpath = query.split("/", 1)[1]
    rawfile = db.path_mirror(_HOST) / relpath
    outfile = db_prefix() / relpath
    if outfile.suffix in (".bed", ".gff"):
        outfile = outfile.with_suffix(outfile.suffix + ".gz")
    elif outfile.name.endswith(".gtf.gz"):
        outfile = outfile.with_suffix("").with_suffix(".gff.gz")
    content = tools.retrieve_content(url, rawfile)
    if outfile.suffix == ".gz":
        future = cli.thread_submit(tools.compress, content, outfile)
    else:
        future = cli.thread_submit(fs.symlink, rawfile, outfile)
    return cli.thread_submit(htslib.try_index, future)


def iter_download_queries():
    for query in iter_download_queries_all():
        if re.search(r"Oryza_sativa_Japonica|Solanum_lycopersicum", query):
            yield query


def iter_download_queries_all():
    content = download_php()
    for mobj in re.finditer(r"download_ftp\.php\?([^\"']+)", content):
        yield mobj[1]


def download_php():
    url = f"http://{_HOST}/download.php"
    cache = db.path_mirror(_HOST) / "download.php.html"
    return tools.retrieve_content(url, cache, force=True).decode()


def db_prefix():
    return db.path("plantregmap")


def rglob(pattern: str, species: str = ".") -> Iterator[Path]:
    for species_dir in db_prefix().iterdir():
        if re.search(species, species_dir.name, re.IGNORECASE):
            yield from species_dir.rglob(pattern)


class FTPplantregmap(ftplib.LazyFTP):
    def __init__(self):
        host = "ftp.cbi.pku.edu.cn"
        super().__init__(
            host,
            "/pub/database/PlantRegMap",
            db.path_mirror(host) / "plantregmap",
            timeout=65535,
        )

    def ls_cache(self, species: str = ""):
        self.nlst_cache("")
        self.nlst_cache("08-download")
        self.nlst_cache("08-download/FTP")
        self.nlst_cache("08-download/FTP/pairwise_alignments")
        if species:
            self.nlst_cache("08-download/{species}")

    def download_pairwise_alignments(self, species: str):
        relpath = f"08-download/FTP/pairwise_alignments/{species}"
        nlst = self.nlst_cache(relpath)
        return [self.retrieve(x, checksize=True) for x in nlst]

    def download_multiple_alignments(self, species: str):
        relpath = f"08-download/{species}/multiple_alignments"
        nlst = self.nlst_cache(relpath)
        return [self.retrieve(x, checksize=True) for x in nlst]

    def download_conservation(self, species: str):
        relpath = f"08-download/{species}/sequence_conservation"
        nlst = self.nlst_cache(relpath)
        return [self.retrieve(x, checksize=True) for x in nlst]


if __name__ == "__main__":
    main()
