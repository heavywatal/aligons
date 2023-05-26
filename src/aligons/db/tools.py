import gzip
import io
import logging
import re
import urllib.request
from pathlib import Path
from zipfile import ZipFile

import polars as pl

from aligons.extern import htslib
from aligons.util import cli, fs

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    split_gff(args.infile)


def retrieve_bgzip(url: str, outfile: Path):
    _log.info(url)
    if not cli.dry_run and fs.is_outdated(outfile):
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        with urllib.request.urlopen(url) as response:  # noqa: S310
            content = response.read()
        if url.endswith(".zip") and outfile.suffix != ".zip":
            content = zip_decompress(content)
        if htslib.to_be_bgzipped(outfile.name):
            assert outfile.suffix == ".gz"
            if url.endswith(".gz"):
                content = gzip.decompress(content)
            content = htslib.bgzip_compress(content)
        with outfile.open("wb") as fout:
            fout.write(content)
        if htslib.to_be_tabixed(outfile.name):
            htslib.tabix(outfile)
    return outfile


def retrieve_cache(url: str, cache: Path) -> bytes:
    if cache.exists():
        with cache.open("rb") as fin:
            return fin.read()
    cache.parent.mkdir(0o755, parents=True, exist_ok=True)
    _log.info(url)
    response = urllib.request.urlopen(url)  # noqa: S310
    content = response.read()
    with cache.open("wb") as fout:
        fout.write(content)
    return content


def zip_decompress(data: bytes) -> bytes:
    with ZipFile(io.BytesIO(data), "r") as zin:
        members = zin.namelist()
        assert len(members) == 1
        return zin.read(members[0])


def split_gff(path: Path):
    body = read_gff_body(path)
    stem = path.stem.removesuffix(".gff").removesuffix(".gff3")
    _log.debug(f"{stem=}")
    for sequence_region in read_gff_sequence_region(path):
        assert (mobj := re.search(r"##sequence-region\s+(\S+)", sequence_region))
        seqid = mobj.group(1)
        outfile = path.parent / f"{stem}.chromosome.{seqid}.gff3.gz"
        print(outfile)
        subdf = body.filter(pl.col("seqid") == seqid)
        assert subdf.height
        if cli.dry_run or not fs.is_outdated(outfile, path):
            continue
        with gzip.open(outfile, "wt") as fout:
            fout.write("##gff-version 3\n")
            fout.write(sequence_region)
            fout.write(subdf.write_csv(has_header=False, separator="\t"))


def read_gff_sequence_region(path: Path):
    lines: list[str] = []
    with gzip.open(path, "rt") as fin:
        gff_version = next(fin)
        assert gff_version.startswith("##gff-version")
        for line in fin:
            if not line.startswith("#"):
                break
            if line.startswith("##sequence-region"):
                lines.append(line)
    assert lines
    return lines


def read_gff_body(path: Path):
    return pl.read_csv(
        path,
        separator="\t",
        comment_char="#",
        has_header=False,
        dtypes=[pl.Utf8],
        new_columns=[
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ],
    )


if __name__ == "__main__":
    main()
