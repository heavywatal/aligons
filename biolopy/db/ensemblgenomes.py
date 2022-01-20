#!/usr/bin/env python3
import os
import sys
from ftplib import FTP
from pathlib import Path

from . import name

LOCAL_DB_ROOT = Path("~/db/ensemblgenomes/plants").expanduser()
VERSION = os.environ["ENSEMBLGENOMES_VERSION"]


def prefix(version: str = VERSION):
    return LOCAL_DB_ROOT / f"release-{version}"


def list_versions():
    return LOCAL_DB_ROOT.glob("release-*")


def list_ftp_species(version: str, format: str = "fasta"):
    with FTP("ftp.ensemblgenomes.org") as ftp:
        path = f"pub/plants/release-{version}/{format}"
        print(ftp.login(), file=sys.stderr)
        print(ftp.getwelcome(), file=sys.stderr)
        print(ftp.cwd(path), file=sys.stderr)
        return ftp.nlst()


def list_all_species(version: str = VERSION):
    cache = prefix(version) / "species.tsv"
    if not cache.exists():
        assert cache.parent.exists(), f"{cache.parent} exists"
        species = list_ftp_species(version)
        with open(cache, "wt") as fout:
            fout.write("\n".join(species) + "\n")
        print(f"Wrote {cache}", file=sys.stderr)
    with open(cache, "r") as fin:
        return [x.rstrip() for x in fin.readlines()]


def list_species(version: str = VERSION, format: str = "fasta"):
    return [x.name for x in (prefix(version) / format).iterdir()]


def list_files(species: str, version: str = VERSION, format: str = "fasta"):
    patt = {"fasta": "**/*.fa.gz", "gff3": "**/*.gff3.gz"}
    return (prefix(version) / format / species).glob(patt[format])


def rglob(pattern: str, species: str = "", format: str = "fasta"):
    if species:
        path = prefix() / format / species
    else:
        path = prefix() / format
    return path.rglob(pattern)


def single_path(pattern: str, species: str = "", format: str = "fasta"):
    found = list(rglob(pattern, species, format))
    assert len(found) == 1
    return found[0]


def main(argv: list[str] | None = None):
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", default=VERSION)
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
    if args.versions:
        for x in list_versions():
            print(x)
        return
    if args.all:
        species = list_all_species(args.version)
    else:
        species = list_species(args.version)
    if args.species:
        species = [x for x in species if name.shorten(x) in args.species]
    if args.format:
        for sp in species:
            for x in list_files(sp, args.version, args.format):
                if args.name:
                    print(x.name)
                else:
                    print(x)
    else:
        for sp in species:
            print(sp)


if __name__ == "__main__":
    main()
