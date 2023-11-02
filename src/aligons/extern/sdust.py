"""sdust: symmetric DUST algorithm for finding low-complexity sequences.

https://github.com/lh3/sdust

src: {basename}.fa(.gz)
dst: {basename}.fa.sdust.bed.gz
"""
import logging
from pathlib import Path

from aligons.util import cli, fs, subp

from . import htslib

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path, nargs="+")
    args = parser.parse_args(argv or None)
    for infile in args.infile:
        cli.thread_submit(run, infile)


def run(in_fa: Path, outdir: Path | None = None) -> Path:
    if outdir is None:
        outdir = in_fa.parent
    out_bed = outdir / (in_fa.name.removeprefix(".gz") + ".sdust.bed.gz")
    is_to_run = fs.is_outdated(out_bed, in_fa) and not cli.dry_run
    p = subp.popen(["sdust", in_fa], if_=is_to_run, stdout=subp.PIPE)
    if is_to_run:
        htslib.bgzip(p.stdout, out_bed)
    _log.info(f"{out_bed}")
    return out_bed


if __name__ == "__main__":
    main()
