"""Soft-mask fasta sequences based on kmer frequencies.

src: {db.api.prefix}/fasta/{species}/*.fa.gz
dst: {db.api.prefix}/fasta/{species}/kmer/*.fa.gz
"""
import concurrent.futures as confu
import gzip
import logging
from pathlib import Path

import polars as pl
import tomli_w

from aligons.db import api
from aligons.util import cli, config, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("species", type=str)
    args = parser.parse_args(argv or None)
    print(run(args.species))


def run(species: str) -> list[Path]:
    genome = api.genome_fa(species)
    jf = count(genome)
    dumpfile = dump(jf)
    histofile = histo(jf)
    threshold = calc_threshold(histofile)
    log_config(histofile, threshold)
    threads = config["jellyfish"]["count"]["threads"]
    with confu.ThreadPoolExecutor(max_workers=threads) as pool:
        fts = [
            pool.submit(mask_genome, chromosome, dumpfile, threshold)
            for chromosome in api.list_chromosome_fa(species)
        ]
    return [f.result() for f in fts]


def count(infile: Path):
    opts = config["jellyfish"]["count"]
    outfile = infile.parent / "kmer" / "mer_counts.jf"
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
    p = subp.popen(args, if_=is_to_run, stdin=subp.PIPE)
    if is_to_run:
        assert p.stdin
        with gzip.open(infile, "rb") as fin:
            p.stdin.write(fin.read())
    p.communicate()
    return outfile


def dump(jffile: Path):
    opts = config["jellyfish"]["dump"]
    outfile = jffile.with_suffix(".dump.fa")
    args: subp.Args = ["jellyfish", "dump"]
    args.extend(["-L", str(opts["lower_count"])])
    args.extend(["-o", outfile])
    args.append(jffile)
    is_to_run = fs.is_outdated(outfile, jffile)
    subp.run(args, if_=is_to_run)
    hard_mask = "N" * config["jellyfish"]["count"]["mer_len"]
    if is_to_run and not cli.dry_run:
        # append NNN to prevent dCNS from reverting n to N
        with outfile.open("a") as fout:
            fout.write(f">65535\n{hard_mask}\n")
    return outfile


def histo(jffile: Path):
    opts = config["jellyfish"]["histo"]
    outfile = jffile.with_suffix(".histo")
    args: subp.Args = ["jellyfish", "histo"]
    args.extend(["--high", str(opts["high"])])
    args.extend(["-o", outfile])
    args.append(jffile)
    subp.run(args, if_=fs.is_outdated(outfile, jffile))
    return outfile


def calc_threshold(histofile: Path) -> int:
    if cli.dry_run:
        return 65535
    cols = ["x", "y"]
    daf = pl.read_csv(histofile, has_header=False, new_columns=cols, separator=" ")
    y = daf["y"]
    ddy = y.diff(null_behavior="drop").diff(null_behavior="drop")
    x = daf["x"][: len(ddy)]
    i = (ddy**2 + x**2).arg_min()
    assert i
    return int(x[i])


def log_config(histofile: Path, freq: int):
    config_log = histofile.parent / "config.toml"
    opts = {}
    opts["jellyfish"] = config["jellyfish"]
    opts["dCNS"] = {"freq": freq}
    _log.debug(tomli_w.dumps(opts))
    if fs.is_outdated(config_log, histofile) and not cli.dry_run:
        with config_log.open("wb") as fout:
            tomli_w.dump(opts, fout)


def mask_genome(infile: Path, kmer_fa: Path, freq: int = 50):
    """https://github.com/baoxingsong/dCNS."""
    dump_lower_count = config["jellyfish"]["dump"]["lower_count"]
    if freq < dump_lower_count:
        _log.warning(f"threshold frequency: {freq} < {dump_lower_count=}")
    output = kmer_fa.parent / infile.name
    tmp_fa = output.with_suffix("")
    args: subp.Args = ["dCNS", "maskGenome", "-d"]
    args.extend(["-i", "/dev/stdin"])
    args.extend(["-o", tmp_fa])
    args.extend(["-k", kmer_fa])
    args.extend(["-f", str(freq)])
    is_to_run = fs.is_outdated(output, [infile, kmer_fa])
    p = subp.popen(args, if_=is_to_run, stdin=subp.PIPE)
    if is_to_run and not cli.dry_run:
        with gzip.open(infile, "rb") as fin:
            p.communicate(fin.read())
    subp.run(["gzip", "-f", tmp_fa], if_=is_to_run)
    return output


if __name__ == "__main__":
    main()
