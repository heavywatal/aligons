"""https://plants.ensembl.org/
"""
import argparse
import functools
import logging
import os
import re
from ftplib import FTP
from pathlib import Path

from .. import cli, fs
from . import name

_log = logging.getLogger(__name__)
LOCAL_DB_ROOT = Path("~/db/ensemblgenomes/plants").expanduser()
VERSION = os.environ["ENSEMBLGENOMES_VERSION"]
PREFIX = LOCAL_DB_ROOT / f"release-{VERSION}"


def main(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-V", "--versions", action="store_true")
    parser.add_argument("-a", "--all", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-c", "--checksums", action="store_true")
    parser.add_argument(
        "-f", "--fasta", action="store_const", const="fasta", dest="format"
    )
    parser.add_argument(
        "-g", "--gff3", action="store_const", const="gff3", dest="format"
    )
    parser.add_argument("--name", action="store_true")
    parser.add_argument("species", nargs="*")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.download:
        download(args.species)
        return
    if args.checksums:
        fs.checksums(PREFIX.rglob("CHECKSUMS"))
        return
    if args.versions:
        for x in sorted(list_versions()):
            print(x)
        return
    if args.all and not args.format:
        species = list_all_species()
    else:
        species = list_species()
    if args.species:
        species = name.filter_by_shortname(species, args.species)
    if not args.format:
        for sp in species:
            print(sp)
        return
    for sp in species:
        for x in list_files(sp, args.format, all=args.all):
            if args.name:
                print(x.name)
            else:
                print(x)


def list_versions():
    _log.debug(f"{LOCAL_DB_ROOT=}")
    return LOCAL_DB_ROOT.glob("release-*")


@functools.cache
def list_all_species():
    cache = PREFIX / "species.tsv"
    if not cache.exists():
        assert cache.parent.exists(), f"{cache.parent} exists"
        with FTPensemblgenomes() as ftp:
            species = ftp.list_species()
        with open(cache, "wt") as fout:
            fout.write("\n".join(species) + "\n")
    _log.debug(f"{cache=}")
    with open(cache, "r") as fin:
        return [x.rstrip() for x in fin.readlines()]


@functools.cache
def list_species(format: str = "fasta"):
    path = PREFIX / format
    _log.debug(f"{path=}")
    return [x.name for x in path.iterdir()]


def list_files(species: str = "", format: str = "fasta", all: bool = False):
    patt = "*" if all else {"fasta": "*.fa.gz", "gff3": "*.gff3.gz"}[format]
    for x in rglob(patt, species, format):
        if x.is_file():
            yield x


def get_file(pattern: str, species: str = "", format: str = "fasta"):
    found = list(rglob(pattern, species, format))
    _log.debug(f"{found=}")
    assert len(found) == 1
    return found[0]


def rglob(pattern: str, species: str = "", format: str = "fasta"):
    path = PREFIX / format / species
    _log.debug(f"({path}).rglob({pattern})")
    return path.rglob(pattern)


def expand_shortnames(shortnames: list[str]):
    return name.filter_by_shortname(list_species(), shortnames)


def download(species: list[str]):
    assert species
    assert not (diff := set(species) - set(list_all_species())), diff
    os.chdir(PREFIX)
    with FTPensemblgenomes() as ftp:
        ftp.cwd_prefix()
        for sp in species:
            ftp.download_fasta(sp)
            ftp.download_gff3(sp)
    # for sp in species:
    #     options = "--include *_sm.chromosome.*.fa.gz --exclude *.gz"
    #     rsync(f"fasta/{sp}/dna", options)
    #     options = "--include *.chromosome.*.gff3.gz --exclude *.gz"
    #     rsync(f"gff3/{sp}", options)


class FTPensemblgenomes(FTP):
    def __init__(self):
        host = "ftp.ensemblgenomes.org"
        _log.info(f"FTP({host})")
        super().__init__(host)
        _log.debug("ftp.login()")
        _log.info(self.login())
        _log.debug("ftp.getwelcome()")
        _log.info(self.getwelcome())

    def cwd_prefix(self, subdir: str = ""):
        prefix = f"/pub/plants/release-{VERSION}"
        path = f"{prefix}/{subdir}".rstrip("/")
        _log.info(f"ftp.cwd({path})")
        _log.info(self.cwd(path))

    def list_species(self, format: str = "fasta"):
        self.cwd_prefix(format)
        _log.info("ftp.nlst()")
        return self.nlst()

    def download_fasta(self, species: str):
        pattern = r"/CHECKSUMS|/README|_sm\.chromosome\..+fa\.gz$"
        for x in self.nlst_search(f"fasta/{species}/dna", pattern):
            print(self.retrieve(x))

    def download_gff3(self, species: str):
        pattern = r"/CHECKSUMS|/README|\.chromosome\..+gff3\.gz$"
        for x in self.nlst_search(f"gff3/{species}", pattern):
            print(self.retrieve(x))

    def nlst_search(self, dir: str, pattern: str):
        _log.info(f"ftp.nlst({dir})")  # ensembl does not support mlsd
        rex = re.compile(pattern)
        for x in self.nlst(dir):
            if rex.search(x):
                yield x

    def retrieve(self, path: str):
        outfile = Path(path)
        if not outfile.exists() and not cli.dry_run:
            outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
            with open(outfile, "wb") as fout:
                cmd = f"RETR {path}"
                _log.info(cmd)
                _log.info(self.retrbinary(cmd, fout.write))
        return outfile

    def quit(self):
        _log.debug("ftp.quit()")
        _log.info(ret := super().quit())
        return ret


def rsync(relpath: str, options: str = ""):
    server = "ftp.ensemblgenomes.org"
    remote_prefix = f"rsync://{server}/all/pub/plants/release-{VERSION}"
    src = f"{remote_prefix}/{relpath}/"
    dst = PREFIX / relpath
    return cli.run(f"rsync -auv {options} {src} {dst}")


if __name__ == "__main__":
    main()
