"""sdust: symmetric DUST algorithm for finding low-complexity sequences.

<https://github.com/lh3/sdust>
"""

import logging
from pathlib import Path

from aligons.util import cli, fs, subp

from . import htslib

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    """CLI for manual execution and testing."""
    parser = cli.ArgumentParser()
    parser.add_argument("infile", type=Path, nargs="+")
    args = parser.parse_args(argv or None)
    for infile in args.infile:
        cli.thread_submit(run, infile)


def run(in_fa: Path, outdir: Path | None = None) -> Path:
    """Run `sdust`.

    :param in_fa: Input FASTA: `{basename}.fa(.gz)`.
    :param outdir: Output directory. The same directory as `in_fa` if `None`.
    :returns: Output BED: `{basename}.fa.sdust.bed.gz`.
    """
    if outdir is None:
        outdir = in_fa.parent
    out_bed = outdir / (in_fa.name.removeprefix(".gz") + ".sdust.bed.gz")
    is_to_run = fs.is_outdated(out_bed, in_fa) and not cli.dry_run
    with subp.popen(["sdust", in_fa], stdout=subp.PIPE, if_=is_to_run) as p:
        htslib.bgzip(p.stdout, out_bed, if_=is_to_run)
    return fs.print_if_exists(out_bed)


if __name__ == "__main__":
    main()
