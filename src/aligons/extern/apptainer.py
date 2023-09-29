import logging
import re
from pathlib import Path

import polars as pl

from aligons import db
from aligons.db import tools
from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)
_galaxy_domain = "depot.galaxyproject.org"
_galaxy_prefix = f"https://{_galaxy_domain}/singularity/"
_galaxy_apps = [
    "bedtools",
    "jbrowse2",
    "jellyfish",
    "lastz",
    "multiz",
    "phast",
    "repeatmasker",
    "samtools",
    "seqkit",
    "trf",
]


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-a", "--all", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    args = parser.parse_args(argv or None)
    tsv = galaxy_index()
    if args.all:
        print(tsv.open("rt").read())
    else:
        table = pl.read_csv(tsv, separator="\t")
        for x in latest_apps(table):
            url = f"{_galaxy_prefix}{x}"
            if args.download:
                wget_nc(url)
            print(x)


def wget_nc(url: str):
    subp.run(["wget", "-nc", url])


def latest_apps(table: pl.DataFrame):
    return (
        table.filter(pl.col("app").is_in(_galaxy_apps))
        .select(
            pl.concat_str([pl.col("app"), pl.col("tag")], separator=":").alias("name")
        )
        .to_series()
        .to_list()
    )


def galaxy_index() -> Path:
    cache_html = _cache_dir() / "singularity.html"
    cache_tsv = cache_html.with_suffix(".tsv")
    _log.info(f"{cache_tsv}")
    if fs.is_outdated(cache_tsv, cache_html):
        content = tools.retrieve_content(_galaxy_prefix, cache_html)
        table = _parse_galaxy_index_html(content)
        table.write_csv(cache_tsv, separator="\t")
    return cache_tsv


def _parse_galaxy_index_html(content: bytes):
    tsv = re.sub(rb"(</a>| ) +", rb"\t", content)
    cols = ["anchor", "time", "size"]
    raw = pl.read_csv(
        tsv, separator="\t", skip_rows=6, has_header=False, new_columns=cols
    )
    return (
        raw.filter(pl.col("anchor").str.starts_with("<a "))
        .with_columns(
            pl.col("time").str.to_datetime("%d-%b-%Y %H:%M"),
            pl.col("anchor")
            .str.extract(r"href=\"([^\"]+)")
            .str.replace("%3A", ":")
            .alias("href"),
        )
        .filter(~pl.col("href").str.starts_with("mulled"))
        .with_columns(
            pl.col("href").str.extract("(.+):").alias("app"),
            pl.col("href").str.extract(":(.+)").alias("tag"),
        )
        .sort(["app", "time"])
        .group_by("app")
        .last()
        .with_columns(pl.col("time").dt.date().alias("date"))
        .select(["app", "tag", "date", "size"])
        .sort("app")
    )


def _cache_dir():
    return db.path_mirror(_galaxy_domain)


if __name__ == "__main__":
    main()
