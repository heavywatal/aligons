import gzip
import logging
import re
from pathlib import Path

import polars as pl

from ..util import cli, fs

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path)
    args = parser.parse_args(argv or None)
    split_gff(args.infile)


def split_gff(path: Path):
    df = read_gff_body(path)
    stem = path.stem.removesuffix(".gff").removesuffix(".gff3")
    _log.debug(f"{stem=}")
    for sequence_region in read_gff_sequence_region(path):
        assert (mobj := re.search(r"##sequence-region\s+(\S+)", sequence_region))
        seqid = mobj.group(1)
        outfile = path.parent / f"{stem}.chromosome.{seqid}.gff3.gz"
        print(outfile)
        subdf = df.filter(pl.col("seqid") == seqid)
        assert subdf.height
        if cli.dry_run or not fs.is_outdated(outfile, path):
            continue
        with gzip.open(outfile, "wt") as fout:
            fout.write("##gff-version 3\n")
            fout.write(sequence_region)
            fout.write(subdf.write_csv(has_header=False, separator="\t"))
    return


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
