"""A fast multi-threaded k-mer counter.

<https://github.com/gmarcais/Jellyfish>

src: *.fa.gz
dst: kmer/*.fa.gz
"""

import logging
from pathlib import Path
from typing import Any

import polars as pl
import tomli_w

from aligons.db import api
from aligons.util import cli, config, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("species", type=str)
    args = parser.parse_args(argv or None)
    fs.print_if_exists(run(args.species))


def run(species: str) -> Path:
    genome = api.genome_fa(species)
    jf = count(genome)
    dumpfile = dump(jf)
    histo_file = histo(jf)
    threshold = calc_threshold(histo_file)
    log_config(histo_file, threshold)
    return mask_genome(genome, dumpfile, threshold)


def count(infile: Path) -> Path:
    opts = config["jellyfish"]["count"]
    outfile = infile.with_name("kmer") / "mer_counts.jf"
    args: subp.Args = ["jellyfish", "count", "--canonical"]
    args.extend(["-m", str(opts["mer_len"])])
    args.extend(["-s", str(opts["size"])])
    args.extend(["-t", str(opts["threads"])])
    args.extend(["-L", str(opts["lower_count"])])
    args.extend(["-o", outfile])
    args.append("/dev/stdin")
    is_to_run = fs.is_outdated(outfile, infile) and not cli.dry_run
    if is_to_run:
        outfile.parent.mkdir(0o755, exist_ok=True)
    with subp.popen_zcat(infile, if_=is_to_run) as zcat:
        subp.run(args, stdin=zcat.stdout, if_=is_to_run)
    return outfile


def dump(jf_file: Path) -> Path:
    opts = config["jellyfish"]["dump"]
    outfile = jf_file.with_suffix(".dump.fa")
    args: subp.Args = ["jellyfish", "dump"]
    args.extend(["-L", str(opts["lower_count"])])
    args.extend(["-o", outfile])
    args.append(jf_file)
    is_to_run = fs.is_outdated(outfile, jf_file)
    subp.run(args, if_=is_to_run)
    hard_mask = "N" * config["jellyfish"]["count"]["mer_len"]
    if is_to_run and not cli.dry_run:
        # append NNN to prevent dCNS from reverting n to N
        with outfile.open("a") as fout:
            fout.write(f">65535\n{hard_mask}\n")
    return outfile


def histo(jf_file: Path) -> Path:
    opts = config["jellyfish"]["histo"]
    outfile = jf_file.with_suffix(".histo")
    args: subp.Args = ["jellyfish", "histo"]
    args.extend(["--high", str(opts["high"])])
    args.extend(["-o", outfile])
    args.append(jf_file)
    subp.run(args, if_=fs.is_outdated(outfile, jf_file))
    return outfile


def calc_threshold(histo_file: Path) -> int:
    if cli.dry_run:
        return 65535
    cols = ["x", "y"]
    daf = pl.read_csv(histo_file, has_header=False, new_columns=cols, separator=" ")
    y = daf["y"]
    ddy = y.diff(null_behavior="drop").diff(null_behavior="drop")
    x = daf["x"][: len(ddy)]
    i = (ddy**2 + x**2).arg_min()
    assert i, histo_file
    return int(x[i])


def log_config(histo_file: Path, freq: int) -> None:
    config_log = histo_file.with_name("config.toml")
    opts: dict[str, Any] = {}
    opts["jellyfish"] = config["jellyfish"]
    opts["dCNS"] = {"freq": freq}
    _log.debug(tomli_w.dumps(opts))
    if fs.is_outdated(config_log, histo_file) and not cli.dry_run:
        with config_log.open("wb") as fout:
            tomli_w.dump(opts, fout)


def mask_genome(infile: Path, kmer_fa: Path, freq: int = 50) -> Path:
    """https://github.com/baoxingsong/dCNS."""
    dump_lower_count = config["jellyfish"]["dump"]["lower_count"]
    if freq < dump_lower_count:
        _log.warning(f"threshold frequency: {freq} < {dump_lower_count=}")
    outfile = kmer_fa.with_name(infile.name)
    args: subp.Args = ["dCNS", "maskGenome", "-d"]
    args.extend(["-i", "/dev/stdin"])
    args.extend(["-o", "/dev/stdout"])
    args.extend(["-k", kmer_fa])
    args.extend(["-f", str(freq)])
    is_to_run = fs.is_outdated(outfile, [infile, kmer_fa])
    with (
        subp.popen_zcat(infile, if_=is_to_run) as zcat,
        subp.popen(args, stdin=zcat.stdout, stdout=subp.PIPE, if_=is_to_run) as pd,
    ):
        subp.gzip(pd.stdout, outfile, if_=is_to_run)
    return outfile


if __name__ == "__main__":
    main()
