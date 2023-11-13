import gzip
import io
import logging
import re
from pathlib import Path
from typing import IO

import polars as pl

from . import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    print(_read_sequence_region(args.infile))


def split_with_fasize(path: Path, fasize: Path) -> list[Path]:
    regions: dict[str, str] = {}
    with fasize.open("rt") as fin:
        for line in fin:
            seqid, length = line.split()
            regions[seqid] = f"##sequence-region {seqid} 1 {length}\n"
    return split_with_hint(path, regions)


def split_with_hint(path: Path, regions: dict[str, str]) -> list[Path]:
    stem = path.stem.removesuffix(".gff").removesuffix(".gff3")
    _log.debug(f"{stem=}")
    regions_gff = _read_sequence_region(path)
    files: list[Path] = []
    body = None
    for seqid, seq_region in regions.items():
        if re.search(r"scaffold|contig", seqid):
            _log.debug(f"ignoring {seqid} in {path}")
            continue
        if (in_gff := regions_gff.get(seqid, "")) and seq_region != in_gff:
            _log.warning(f"'{seq_region}' != '{in_gff}'")
        outfile = path.with_name(f"{stem}.chromosome.{seqid}.gff3.gz")
        files.append(outfile)
        _log.info(f"{outfile}")
        if not fs.is_outdated(outfile, path) or cli.dry_run:
            continue
        if body is None:
            body = _read_body(path).lazy()
        data = body.filter(pl.col("seqid") == seqid).sort(["start"]).collect()
        with gzip.open(outfile, "wt") as fout:
            fout.write("##gff-version 3\n")
            fout.write(seq_region)
            fout.write(data.write_csv(include_header=False, separator="\t"))
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
        .filter(pl.col("start") < pl.col("end"))  # SL2.40: end 338074 < begin 338105
        .sort(["seqid", "start"])
        .write_csv(bio, include_header=False, separator="\t")
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


def sort_clean_chromosome(infile: Path, stdout: IO[bytes]) -> None:
    # TODO: jbrowse2 still needs billzt/gff3sort precision?
    cmd1 = f"zstdgrep -v '^#' {infile!s}"
    cmd2 = "grep -v '\tchromosome\t'"
    cmd3 = "sort -k4,4n"
    with (
        subp.popen(cmd1, stdout=subp.PIPE, quiet=True) as v_hash,
        subp.popen(cmd2, stdin=v_hash.stdout, stdout=subp.PIPE, quiet=True) as v_chr,
        subp.popen(cmd3, stdin=v_chr.stdout, stdout=stdout, quiet=True) as sort,
    ):
        sort.communicate()


if __name__ == "__main__":
    main()
