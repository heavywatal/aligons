import io
import logging
import re
import typing
from pathlib import Path

import polars as pl

from . import cli, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    x = GFF(args.infile)
    print(x.sanitize().to_string())


class GFF:
    def __init__(self, source: Path) -> None:
        self.header: list[bytes] = _read_header(source)
        self.body: pl.LazyFrame = _read_body(source)
        self._workaround(source)

    def sanitize(self) -> "GFF":
        self.body = (
            self.body.sort([pl.col("seqid").str.pad_start(4, "0"), "start"])
            # remove chromosome for visualization
            .filter(pl.col("type") != "chromosome")
        )
        return self

    def seqid_replace(self, pattern: str, repl: str) -> "GFF":
        self.body = self.body.with_columns(pl.col("seqid").str.replace(pattern, repl))
        return self

    def write(self, iobytes: typing.IO[bytes]) -> None:
        iobytes.writelines(self.header)
        iobytes = typing.cast(io.BytesIO, iobytes)
        self.body.collect().write_csv(iobytes, include_header=False, separator="\t")

    def to_string(self) -> str:
        buffer = io.BytesIO()
        self.write(buffer)
        return buffer.getvalue().decode()

    def _workaround(self, source: Path) -> None:
        if (
            source.name.startswith("PGSC_DM_V403")
            and source.resolve().parent.parent.name == "spuddb.uga.edu"
        ):
            self.seqid_replace("^ST4.03ch", "chr")
        if source.name.startswith("ITAG2.3"):
            # SL2.40ch01: begin 338105 > end 338074
            self.body = self.body.filter(pl.col("start") < pl.col("end"))


def _read_header(source: Path) -> list[bytes]:
    lines: list[bytes] = []
    with subp.popen_zcat(source) as zcat:
        assert zcat.stdout, source
        for line in zcat.stdout:
            if not line.startswith(b"#"):
                break
            lines.append(line)
    if source.name.startswith("ITAG4.1"):
        # solgenomics SL4.0 has dirty headers without space
        _log.info(f"{source}:\n{b''.join(lines)}")
        rex = re.compile(rb"^(##sequence-region)(\S)")
        lines = [rex.sub(rb"\1 \2", x) for x in lines]
    if not lines or not lines[0].startswith(b"##gff-version"):
        _log.warning(f"{source}:invalid GFF without ##gff-version")
        lines = [b"##gff-version 3\n", *lines]
    return lines


def _read_body(source: Path | bytes | io.BytesIO) -> pl.LazyFrame:
    if isinstance(source, bytes):
        source = re.sub(rb"\n\n+", rb"\n", source)
    # scan_csv does not read compressed files
    if isinstance(source, Path) and source.suffix == ".zip":
        source = subp.run_zcat(source).stdout
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
    ).lazy()


if __name__ == "__main__":
    main()
