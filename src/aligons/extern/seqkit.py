"""https://bioinf.shenwei.me/seqkit
"""
import logging
from pathlib import Path

from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("file", type=Path)
    args = parser.parse_args(argv or None)
    split(args.file)


def split(path: Path):
    """https://bioinf.shenwei.me/seqkit/usage/#split
    dir/seq.fa.gz -> dir/.seqkit/seq.part_{id}.fasta.gz
    """
    outdir = path.parent / ".seqkit"
    args: subp.Args = ["seqkit", "split"]
    args.extend(["-2", "-i", "-e", ".gz", "-O", outdir, path])
    subp.run_if(fs.is_outdated(outdir, path), args)
    return outdir


if __name__ == "__main__":
    main()
