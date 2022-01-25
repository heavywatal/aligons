import logging
import os
from ftplib import FTP
from pathlib import Path

from . import name
from .. import cli

_log = logging.getLogger(__name__)
LOCAL_DB_ROOT = Path("~/db/ensemblgenomes/plants").expanduser()
VERSION = os.environ["ENSEMBLGENOMES_VERSION"]


def prefix(version: str = VERSION):
    return LOCAL_DB_ROOT / f"release-{version}"


def list_versions():
    _log.debug(f"{LOCAL_DB_ROOT=}")
    return LOCAL_DB_ROOT.glob("release-*")


def list_ftp_species(version: str, format: str = "fasta"):
    server = "ftp.ensemblgenomes.org"
    path = f"pub/plants/release-{version}/{format}"
    _log.info(f"FTP({server})")
    with FTP(server) as ftp:
        _log.debug("ftp.login()")
        _log.info(ftp.login())
        _log.debug("ftp.getwelcome()")
        _log.info(ftp.getwelcome())
        _log.info(f"ftp.cwd({path})")
        _log.info(ftp.cwd(path))
        _log.info("ftp.nlst()")
        return ftp.nlst()


def list_all_species(version: str = VERSION):
    cache = prefix(version) / "species.tsv"
    if not cache.exists():
        assert cache.parent.exists(), f"{cache.parent} exists"
        species = list_ftp_species(version)
        with open(cache, "wt") as fout:
            fout.write("\n".join(species) + "\n")
    _log.debug(f"{cache=}")
    with open(cache, "r") as fin:
        return (x.rstrip() for x in fin.readlines())


def list_species(version: str = VERSION, format: str = "fasta"):
    path = prefix(version) / format
    _log.debug(f"{path=}")
    return (x.name for x in path.iterdir())


def list_files(
    species: str = "", format: str = "fasta", version: str = VERSION, all: bool = False
):
    patt = "*" if all else {"fasta": "*.fa.gz", "gff3": "*.gff3.gz"}[format]
    for x in rglob(patt, species, format, version):
        if x.is_file():
            yield x


def get_file(
    pattern: str, species: str = "", format: str = "fasta", version: str = VERSION
):
    found = list(rglob(pattern, species, format, version))
    assert len(found) == 1
    return found[0]


def rglob(
    pattern: str, species: str = "", format: str = "fasta", version: str = VERSION
):
    path = prefix(version) / format / species
    _log.debug(f"({path}).rglob({pattern})")
    return path.rglob(pattern)


def expand_shortnames(shortnames: list[str]):
    return name.filter_by_shortname(list_species(), shortnames)


def main(argv: list[str] | None = None):
    import argparse

    parser = argparse.ArgumentParser(parents=[cli.logging_argparser("v")])
    parser.add_argument("-r", "--version", default=VERSION)
    parser.add_argument("-V", "--versions", action="store_true")
    parser.add_argument("-a", "--all", action="store_true")
    parser.add_argument(
        "-f", "--fasta", action="store_const", const="fasta", dest="format"
    )
    parser.add_argument(
        "-g", "--gff3", action="store_const", const="gff3", dest="format"
    )
    parser.add_argument("-n", "--name", action="store_true")
    parser.add_argument("species", nargs="*")
    args = parser.parse_args()
    cli.logging_config(args.loglevel)
    if args.versions:
        for x in sorted(list_versions()):
            print(x)
        return
    if args.all and not args.format:
        species = list_all_species(args.version)
    else:
        species = list_species(args.version)
    if args.species:
        species = name.filter_by_shortname(species, args.species)
    if not args.format:
        for sp in species:
            print(sp)
        return
    for sp in species:
        for x in list_files(sp, args.format, args.version, all=args.all):
            if args.name:
                print(x.name)
            else:
                print(x)


if __name__ == "__main__":
    main()
