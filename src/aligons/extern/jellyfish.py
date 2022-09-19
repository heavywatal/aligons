"""Soft-mask fasta sequences based on kmer frequencies

src: {ensemblgenomes.prefix}/fasta/{species}/dna/*.fa.gz
dst: {ensemblgenomes.prefix}/fasta/{species}/dna/kmer/*.fa.gz
"""
import concurrent.futures as confu
import gzip
import logging
from pathlib import Path

import numpy as np
import tomli_w

from ..db import ensemblgenomes
from ..util import cli, config, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("species", type=str)
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    outdir = run(args.species)
    print(outdir)


def run(species: str):
    genome = ensemblgenomes.get_file("*.genome.fa.gz", species)
    jf = count(genome)
    dumpfile = dump(jf)
    histofile = histo(jf)
    threshold = calc_threshold(histofile)
    log_config(histofile, threshold)
    threads = config["jellyfish"]["count"]["threads"]
    with confu.ThreadPoolExecutor(max_workers=threads) as pool:
        for chromosome in ensemblgenomes.glob("*.chromosome.*.fa.gz", [species]):
            pool.submit(mask_genome, chromosome, dumpfile, threshold)
    return jf.parent


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
    p = subp.popen_if(is_to_run, args, stdin=subp.PIPE)
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
    subp.run_if(is_to_run, args)
    hard_mask = "N" * config["jellyfish"]["count"]["mer_len"]
    if is_to_run and not cli.dry_run:
        # append NNN to prevent dCNS from reverting n to N
        with open(outfile, "a") as fout:
            fout.write(f">65535\n{hard_mask}\n")
    return outfile


def histo(jffile: Path):
    opts = config["jellyfish"]["histo"]
    outfile = jffile.with_suffix(".histo")
    args: subp.Args = ["jellyfish", "histo"]
    args.extend(["--high", str(opts["high"])])
    args.extend(["-o", outfile])
    args.append(jffile)
    subp.run_if(fs.is_outdated(outfile, jffile), args)
    return outfile


def calc_threshold(histofile: Path):
    if cli.dry_run:
        return 65535
    arr = np.loadtxt(histofile, delimiter=" ", dtype=int)
    y = arr[:, 1]
    ddy = np.diff(np.diff(y))  # type: ignore
    x = arr[: len(ddy), 0]
    i = np.argmin(ddy**2 + x**2)  # type: ignore
    return int(arr[i, 0])


def log_config(histofile: Path, freq: int):
    config_log = histofile.parent / "config.toml"
    opts = {}
    opts["jellyfish"] = config["jellyfish"]
    opts["dCNS"] = {"freq": freq}
    _log.debug(tomli_w.dumps(opts))
    if fs.is_outdated(config_log, histofile) and not cli.dry_run:
        with open(config_log, "wb") as fout:
            tomli_w.dump(opts, fout)


def mask_genome(input: Path, kmer_fa: Path, freq: int = 50):
    """https://github.com/baoxingsong/dCNS"""
    dump_lower_count = config["jellyfish"]["dump"]["lower_count"]
    if freq < dump_lower_count:
        _log.warning(f"threshold frequency: {freq} < {dump_lower_count=}")
    output = kmer_fa.parent / input.name
    tmp_fa = output.with_suffix("")
    args: subp.Args = ["dCNS", "maskGenome", "-d"]
    args.extend(["-i", "/dev/stdin"])
    args.extend(["-o", tmp_fa])
    args.extend(["-k", kmer_fa])
    args.extend(["-f", str(freq)])
    is_to_run = fs.is_outdated(output, [input, kmer_fa])
    p = subp.popen_if(is_to_run, args, stdin=subp.PIPE)
    if is_to_run and not cli.dry_run:
        with gzip.open(input, "rb") as fin:
            p.communicate(fin.read())
    subp.run_if(is_to_run, ["gzip", "-f", tmp_fa])
    return output


if __name__ == "__main__":
    main()
