"""http://plantregmap.gao-lab.org/
"""
import logging
import re
import urllib.request

from ..util import cli, config, fs, subp

_log = logging.getLogger(__name__)
LOCAL_DB_ROOT = config["db"]["root"] / "plantregmap"
HOST = "plantregmap.gao-lab.org"


def main(argv: list[str] | None = None):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("species", nargs="?", default=".")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.download:
        for query in iter_download_queries():
            print(download(query))
    for x in fs.sorted_naturally(rglob("*.*", args.species)):
        print(x)


def rglob(pattern: str, species: str = "."):
    for dir in LOCAL_DB_ROOT.iterdir():
        if re.search(species, dir.name, re.IGNORECASE):
            for x in dir.rglob(pattern):
                yield x


def download(query: str):
    url = f"http://{HOST}/download_ftp.php?{query}"
    path = LOCAL_DB_ROOT / query.split("/", 1)[1]
    if path.suffix in (".bed", ".gff", ".txt"):
        path = path.parent / (path.name + ".gz")
        cmd = f"wget -O- {url}"
    elif path.name.endswith(".gtf.gz"):
        path = path.with_suffix("").with_suffix(".gff.gz")
        cmd = f"wget -O- {url}"
        cmd += " | gunzip -c"
    else:
        cmd = f"wget -O {path} {url}"
    if path.with_suffix("").suffix in (".bed", ".gff"):
        cmd += f" | bgzip -@2 -c >{path}"
        cmd += f"; tabix --csi {path}"
    elif path.with_suffix("").suffix in (".txt"):
        cmd += f" | gzip -c >{path}"
    path.parent.mkdir(0o755, parents=True, exist_ok=True)
    subp.run_if(not path.exists(), cmd, shell=True)
    return path


def iter_download_queries():
    for query in iter_download_queries_all():
        if re.search(r"Oryza_sativa_Japonica|Solanum_lycopersicum", query):
            yield query


def iter_download_queries_all():
    content = download_php()
    for mobj in re.finditer(r"download_ftp\.php\?([^\"']+)", content):
        yield mobj[1]


def download_php():
    cache = LOCAL_DB_ROOT / "download.php.html"
    if cache.exists():
        with cache.open("r") as fin:
            content = fin.read()
    else:
        url = f"http://{HOST}/download.php"
        _log.info(url)
        response = urllib.request.urlopen(url)
        content = response.read().decode()
        with cache.open("w") as fout:
            fout.write(content)
    return content


if __name__ == "__main__":
    main()
