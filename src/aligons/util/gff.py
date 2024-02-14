import io
import logging
import re
import typing
from collections.abc import Iterable, Iterator
from pathlib import Path

import polars as pl

from . import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    print(_read_sequence_region(args.infile))


class GFF:
    def __init__(self, source: Path) -> None:
        self.header = _read_header(source)
        self.body = _read_body(source)

    def write(self, iobytes: typing.IO[bytes]) -> None:
        iobytes.writelines(self.header)
        iobytes = typing.cast(io.BytesIO, iobytes)
        self.body.write_csv(iobytes, include_header=False, separator="\t")


def sort(content: bytes) -> bytes:
    return _extract_header(content) + _sort_body(content)


def collect_sequence_region(infiles: Iterable[Path]) -> list[bytes]:
    lines: list[bytes] = []
    for file in infiles:
        for line in _read_sequence_region(file).values():
            lines.append(line)  # noqa: PERF402
    return lines


def iter_split(
    path: Path, regions: dict[str, bytes]
) -> Iterator[tuple[Path, io.BytesIO]]:
    stem = path.stem.removesuffix(".gff").removesuffix(".gff3")
    _log.debug(f"{stem=}")
    regions_gff = _read_sequence_region(path)
    files: list[Path] = []
    body = None
    for seqid, seq_region_exp in regions.items():
        if re.search(r"scaffold|contig", seqid):
            _log.debug(f"ignoring {seqid} in {path}")
            continue
        seq_region = regions_gff.get(seqid, seq_region_exp)
        if seq_region.split() != seq_region_exp.split():
            _log.warning(f"'{seq_region}' != '{seq_region_exp}'")
        outfile = path.with_name(f"{stem}.chromosome.{seqid}.gff3.gz")
        files.append(outfile)
        _log.info(f"{outfile}")
        if not fs.is_outdated(outfile, path) or cli.dry_run:
            continue
        if body is None:
            body = _read_body(path).lazy()
        data = body.filter(pl.col("seqid") == seqid).sort(["start"]).collect()
        buffer = io.BytesIO()
        buffer.write(b"##gff-version 3\n")
        buffer.write(seq_region)
        data.write_csv(buffer, include_header=False, separator="\t")
        yield outfile, buffer


def _read_sequence_region(path: Path) -> dict[str, bytes]:
    if not path.exists():
        return {}
    lines = _read_header(path)
    lines.pop(0)  # gff-version
    regions: dict[str, bytes] = {}
    comments: list[bytes] = []
    # solgenomics has dirty headers without space: ##sequence-regionSL4.0ch01
    pattern = re.compile(rb"##sequence-region\s*(.+)", re.S)
    for line in lines:
        if mobj := pattern.match(line):
            value = mobj.group(1)
            seqid = value.split(maxsplit=1)[0].decode()
            regions[seqid] = mobj.group(0)
        else:
            comments.append(line)
    if not regions:
        _log.info(f"{path}:unfriendly GFF without ##sequence-region")
    if comments:
        ignored = b"".join(comments)
        _log.warning(f"{path}:ignoring comments in GFF:\n{ignored}")
    return regions


def _read_header(source: Path) -> list[bytes]:
    lines: list[bytes] = []
    with subp.popen_zcat(source) as zcat:
        assert zcat.stdout, source
        for line in zcat.stdout:
            if not line.startswith(b"#"):
                break
            lines.append(line)
    if not lines or not lines[0].startswith(b"##gff-version"):
        _log.warning(f"{source}:invalid GFF without ##gff-version")
        lines = [b"##gff-version 3\n", *lines]
    return lines


def _extract_header(content: bytes) -> bytes:
    if m := re.match(rb"#.+?(?=^[^#])", content, re.M | re.S):
        header = m.group(0)
    else:
        header = b""
    if not header.startswith(b"##gff-version"):
        _log.warning("invalid GFF without ##gff-version")
        header = b"##gff-version 3\n" + header
    return header


def _sort_body(content: bytes) -> bytes:
    bio = io.BytesIO()
    (
        _read_body(content)
        .lazy()
        .filter(pl.col("start") < pl.col("end"))  # SL2.40: end 338074 < begin 338105
        .sort([pl.col("seqid").str.pad_start(4, "0"), "start"])
        .collect()
        .write_csv(bio, include_header=False, separator="\t")
    )
    return bio.getvalue()


def _read_body(source: Path | str | bytes) -> pl.DataFrame:
    if isinstance(source, bytes):
        source = re.sub(rb"\n\n+", rb"\n", source)
    # scan_csv does not read compressed files
    return pl.read_csv(
        source,
        separator="\t",
        comment_prefix="#",
        has_header=False,
        infer_schema_length=0,
        dtypes={"start": pl.UInt64, "end": pl.UInt64},
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


def sort_clean_chromosome(infile: Path, stdout: typing.IO[bytes]) -> None:
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
