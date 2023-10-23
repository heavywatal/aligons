"""https://bioinf.shenwei.me/seqkit."""
import gzip
import logging
from collections.abc import Iterable
from pathlib import Path

from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("file", type=Path)
    args = parser.parse_args(argv or None)
    split(args.file)


def split(path: Path, *, compress: bool = True) -> Path:
    """https://bioinf.shenwei.me/seqkit/usage/#split.

    dir/seq.fa.gz -> dir/_work/seq.part_{id}.fasta.gz
    """
    outdir = path.with_name("_work")
    args: subp.Args = ["seqkit", "split", "--by-id", "-2"]
    if compress:
        args.extend(["-e", ".gz"])
    args.extend(["-O", outdir, path])
    subp.run(args, if_=fs.is_outdated(outdir, path))
    return outdir


def seq_line_width(infile: Path, width: int) -> bytes:
    """https://bioinf.shenwei.me/seqkit/usage/#seq."""
    args: subp.Args = ["seqkit", "seq", "--line-width", str(width), infile]
    return subp.run(args, stdout=subp.PIPE).stdout


def read_fasta_line_width(infile: Path) -> int:
    if infile.suffix == ".gz":
        with gzip.open(infile, "rt") as fin:
            return _fasta_line_width(fin)
    with infile.open("rt") as fin:
        return _fasta_line_width(fin)


def _fasta_line_width(lines: Iterable[str]) -> int:
    for line in lines:
        if not line.startswith(">"):
            return len(line.rstrip())
    _log.warning("invalid fasta")
    return 60


if __name__ == "__main__":
    main()
