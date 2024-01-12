import logging
import re
from pathlib import Path

import polars as pl

from aligons.db import api
from aligons.util import cli, dl, fs

_log = logging.getLogger(__name__)
_galaxy_domain = "depot.galaxyproject.org"
_galaxy_prefix = f"https://{_galaxy_domain}/singularity/"
_galaxy_apps = {
    "bedtools": [],
    "jbrowse2": ["jbrowse"],
    "jellyfish": [],
    "last": ["maf-convert"],
    "lastz": [],
    "multiz": ["multiz", "roast", "maf_project"],
    "phast": ["phastCons", "phyloP", "phyloFit", "phyloBoot", "msa_view"],
    "repeatmasker": ["RepeatMasker", "/usr/local/share/RepeatMasker/famdb.py"],
    "samtools": ["samtools", "bgzip", "tabix"],
    "seqkit": [],
    "trf": [],
    "ucsc-axtchain": ["axtChain"],
    "ucsc-axtsort": ["axtSort"],
    "ucsc-axttomaf": ["axtToMaf"],
    "ucsc-chainmergesort": ["chainMergeSort"],
    "ucsc-chainnet": ["chainNet"],
    "ucsc-chainprenet": ["chainPreNet"],
    "ucsc-fasize": ["faSize"],
    "ucsc-fatotwobit": ["faToTwoBit"],
    "ucsc-netsyntenic": ["netSyntenic"],
    "ucsc-nettoaxt": ["netToAxt"],
    "ucsc-wigtobigwib": ["wigToBigWig"],
}


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-a", "--all", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("--prefix", type=Path, default=Path())
    args = parser.parse_args(argv or None)
    if args.all:
        tsv = galaxy_index()
        print(tsv.open("rt").read())
    elif args.download:
        pull_galaxy(args.prefix)


def pull_galaxy(prefix: Path) -> None:
    tsv = galaxy_index()
    table = pl.read_csv(tsv, separator="\t")
    imgdir = prefix / "biocontainers"
    bindir = prefix / "bin"
    if not cli.dry_run:
        imgdir.mkdir(parents=True, exist_ok=True)
        bindir.mkdir(parents=True, exist_ok=True)
    for x in latest_apps(table):
        url = f"{_galaxy_prefix}{x}"
        sif = dl.fetch(url, imgdir / x).path
        key = sif.name.split(":", 1)[0]
        cmds = _galaxy_apps[key] or [key]
        for command in cmds:
            sh = make_sh(sif, command, bindir)
            _log.info(f"{sh}")


def make_sh(sif: Path, command: str = "", outdir: Path = Path()) -> Path:
    if not command:
        command = sif.name.split(":", 1)[0]
    name = Path(command).name
    outfile = outdir / name
    apptainer_exec = f"apptainer exec {sif.absolute()} {command} $@"
    _log.debug(apptainer_exec)
    if not cli.dry_run:
        with outfile.open("w") as fout:
            fout.write(f"#!/bin/sh\n{apptainer_exec}\n")
        outfile.chmod(0o775)
    return outfile


def latest_apps(table: pl.DataFrame) -> list[str]:
    return (
        table.filter(pl.col("app").is_in(_galaxy_apps.keys()))
        .select(name=pl.concat_str(["app", "tag"], separator=":"))
        .to_series()
        .to_list()
    )


def galaxy_index() -> Path:
    cache_html = _cache_dir() / "singularity.html"
    cache_tsv = cache_html.with_suffix(".tsv")
    _log.info(f"{cache_tsv}")
    if fs.is_outdated(cache_tsv, cache_html):
        content = dl.fetch(_galaxy_prefix, cache_html).content
        table = _parse_galaxy_index_html(content)
        table.write_csv(cache_tsv, separator="\t")
    return cache_tsv


def _parse_galaxy_index_html(content: bytes) -> pl.DataFrame:
    tsv = re.sub(rb"(</a>| ) +", rb"\t", content)
    cols = ["anchor", "time", "size"]
    raw = pl.read_csv(
        tsv, separator="\t", skip_rows=6, has_header=False, new_columns=cols
    )
    return (
        raw.filter(pl.col("anchor").str.starts_with("<a "))
        .with_columns(
            pl.col("time").str.to_datetime("%d-%b-%Y %H:%M"),
            href=pl.col("anchor")
            .str.extract(r"href=\"([^\"]+)")
            .str.replace("%3A", ":"),
        )
        .filter(~pl.col("href").str.starts_with("mulled"))
        .with_columns(
            app=pl.col("href").str.extract("(.+):"),
            tag=pl.col("href").str.extract(":(.+)"),
        )
        .sort(["app", "time"])
        .group_by("app")
        .last()
        .with_columns(date=pl.col("time").dt.date())
        .select(["app", "tag", "date", "size"])
        .sort("app")
    )


def _cache_dir() -> Path:
    return api.prefix(_galaxy_domain)


if __name__ == "__main__":
    main()
