"""Experimental support for Apptainer containers to run external tools."""

import logging
import re
from pathlib import Path

import polars as pl

from aligons.db import api
from aligons.util import cli, config, dl, fs

_log = logging.getLogger(__name__)
_galaxy_domain = "depot.galaxyproject.org"
_galaxy_prefix = f"https://{_galaxy_domain}/singularity/"
_galaxy_apps = config["apptainer"]["galaxy_apps"]


def main(argv: list[str] | None = None) -> None:
    """CLI for downloading Apptainer images from the Galaxy."""
    parser = cli.ArgumentParser()
    parser.add_argument("-a", "--all", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("--prefix", type=Path, default=Path())
    args = parser.parse_args(argv or None)
    if args.all:
        tsv = galaxy_index()
        _log.info(tsv.open("rt").read())
    elif args.download:
        pull_galaxy(args.prefix)


def pull_galaxy(prefix: Path) -> None:
    """Download selected Apptainer images from the Galaxy."""
    tsv = galaxy_index()
    table = pl.read_csv(tsv, separator="\t").lazy()
    img_dir = prefix / "biocontainers"
    bindir = prefix / "bin"
    if not cli.dry_run:
        img_dir.mkdir(parents=True, exist_ok=True)
        bindir.mkdir(parents=True, exist_ok=True)
    for x in latest_apps(table):
        key = x.split(":", 1)[0]
        if key.startswith("ucsc-") or key in ("phast", "repeatmasker"):
            continue
        url = f"{_galaxy_prefix}{x}"
        sif = dl.fetch(url, img_dir / x).path
        cmds = _galaxy_apps[key]
        name = cmds[0] if cmds else key
        sh = make_sh(sif, name, bindir)
        _log.info(sh)
        for command in cmds:
            if command != sh.name:
                fs.symlink(Path(sh.name), bindir / command)


def make_sh(sif: Path, command: str = "", outdir: Path = Path()) -> Path:
    """Generate a shell script to run an Apptainer container.

    :param sif: Apptainer image file.
    :param command: Command name. The image name is used if empty.
    :param outdir: Output directory for the shell script.
    :returns: The generated shell script.
    """
    if not command:
        command = sif.name.split(":", 1)[0]
    name = Path(command).name
    outfile = outdir / name
    content = (
        "#!/bin/sh\n"
        "APP=${0##*/}\n"
        "SRC_DIR=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)\n"
        "IMG_DIR=${SRC_DIR}/../biocontainers\n"
        f"apptainer exec ${{IMG_DIR}}/{sif.name} $APP $@\n"
    )
    _log.debug(content)
    if not cli.dry_run:
        with outfile.open("w") as fout:
            fout.write(content)
        outfile.chmod(0o775)
    return outfile


def latest_apps(table: pl.LazyFrame) -> list[str]:
    """Filter lines with the latest versions of the selected apps.

    :param table: LazyFrame from `_parse_galaxy_index_html()`.
    :returns: A list of `{app}:{tag}` strings.
    """
    return (
        table.filter(pl.col("app").is_in(_galaxy_apps.keys()))
        .select(name=pl.concat_str(["app", "tag"], separator=":"))
        .collect()
        .to_series()
        .to_list()
    )


def galaxy_index() -> Path:
    """Fetch and cache the Galaxy Apptainer image index.

    :returns: Cached TSV file: `{_cache_dir()}/singularity.tsv`.
    """
    cache_html = _cache_dir() / "singularity.html"
    cache_tsv = cache_html.with_suffix(".tsv")
    if fs.is_outdated(cache_tsv, cache_html):
        content = dl.fetch(_galaxy_prefix, cache_html).content
        table = _parse_galaxy_index_html(content)
        table.collect().write_csv(cache_tsv, separator="\t")
    return fs.print_if_exists(cache_tsv)


def _parse_galaxy_index_html(content: bytes) -> pl.LazyFrame:
    tsv = re.sub(rb"(</a>| ) +", rb"\t", content)
    cols = ["anchor", "time", "size"]
    raw = pl.read_csv(
        tsv, separator="\t", skip_rows=6, has_header=False, new_columns=cols
    ).lazy()
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
