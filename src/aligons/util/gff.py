import gzip
import io
import logging
import re
from pathlib import Path

import polars as pl

from aligons.util import cli, fs

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    split_by_seqid(args.infile)


def split_by_seqid(path: Path) -> list[Path]:
    regions = _read_sequence_region(path)
    stem = path.stem.removesuffix(".gff").removesuffix(".gff3")
    if cli.dry_run:
        return [path.with_name(f"{stem}.chromosome.{x}.gff3.gz") for x in regions]
    files: list[Path] = []
    _log.debug(f"{stem=}")
    body = _read_body(path)
    for name, data in body.group_by("seqid", maintain_order=True):
        seqid = str(name)
        if re.search(r"scaffold|contig", seqid):
            _log.debug(f"ignoring {seqid} in {path}")
            continue
        outfile = path.with_name(f"{stem}.chromosome.{seqid}.gff3.gz")
        files.append(outfile)
        _log.info(f"{outfile}")
        if not fs.is_outdated(outfile, path):
            continue
        with gzip.open(outfile, "wt") as fout:
            fout.write("##gff-version 3\n")
            fout.write(regions.get(seqid, ""))
            fout.write(data.sort(["start"]).write_csv(has_header=False, separator="\t"))
    return files


def sort(content: bytes) -> bytes:
    return _extract_header(content) + _sort_body(content)


def _read_sequence_region(path: Path) -> dict[str, str]:
    lines: list[str] = []
    if not path.exists():
        return {}
    with gzip.open(path, "rt") as fin:
        for line in fin:
            if not line.startswith("#"):
                break
            lines.append(line)
    if not lines or not lines[0].startswith("##gff-version"):
        _log.warning(f"{path}:invalid GFF without ##gff-version")
    else:
        lines.pop(0)
    regions: dict[str, str] = {}
    comments: list[str] = []
    # solgenomics has dirty headers without space: ##sequence-regionSL4.0ch01
    pattern = re.compile(r"(##sequence-region)\s*(.+)", re.S)
    for line in lines:
        if mobj := pattern.match(line):
            value = mobj.group(2)
            regions[value.split()[0]] = " ".join(mobj.groups())
        else:
            comments.append(line)
    if not regions:
        _log.info(f"{path}:unfriendly GFF without ##sequence-region")
    if comments:
        ignored = "".join(comments)
        _log.warning(f"{path}:ignoring comments in GFF:\n{ignored}")
    return regions


def _extract_header(content: bytes) -> bytes:
    if m := re.match(rb"##gff-version.+?(?=^[^#])", content, re.M | re.S):
        return m.group(0)
    _log.warning("invalid GFF without ##gff-version")
    return b"##gff-version 3\n"


def _sort_body(content: bytes) -> bytes:
    bio = io.BytesIO()
    (
        _read_body(content)
        .sort(["seqid", "start"])
        .write_csv(bio, has_header=False, separator="\t")
    )
    return bio.getvalue()


def _read_body(source: Path | str | bytes) -> pl.DataFrame:
    if isinstance(source, bytes):
        source = re.sub(rb"\n\n+", rb"\n", source)
    return pl.read_csv(
        source,
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
