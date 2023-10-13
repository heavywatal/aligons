"""Pairwise genome alignment.

src: {db.api.prefix}/fasta/{species}/*.fa.gz
dst: ./pairwise/{target}/{query}/{chromosome}/sing.maf

https://lastz.github.io/lastz/
"""
import concurrent.futures as confu
import gzip
import logging
from pathlib import Path

from aligons.db import api
from aligons.util import cli, config, fs, read_config, subp

from . import kent

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("target", choices=api.species_names())
    parser.add_argument("query", nargs="*")
    args = parser.parse_args(argv or None)
    if args.config:
        read_config(args.config)
    run(args.target, args.query)


def run(target: str, queries: list[str]):
    queries = api.sanitize_queries(target, queries)
    futures: list[confu.Future[Path]] = []
    for query in queries:
        pa = PairwiseAlignment(target, query)
        futures.extend(pa.run())
    cli.wait_raise(futures)
    return Path("pairwise") / target


class PairwiseAlignment:
    def __init__(self, target: str, query: str):
        self._target = target
        self._query = query
        self._target_sizes = api.fasize(target)
        self._query_sizes = api.fasize(query)
        self._outdir = Path("pairwise") / target / query

    def run(self):
        pool = cli.ThreadPool()
        if not cli.dry_run:
            self._outdir.mkdir(0o755, parents=True, exist_ok=True)
        target_chromosomes = api.list_chromosome_2bit(self._target)
        query_chromosomes = api.list_chromosome_2bit(self._query)
        flists: list[list[confu.Future[Path]]] = [
            [pool.submit(self.align_chr, t, q) for q in query_chromosomes]
            for t in target_chromosomes
        ]
        return [pool.submit(self.wait_integrate, futures) for futures in flists]

    def align_chr(self, t2bit: Path, q2bit: Path):
        axtgz = lastz(t2bit, q2bit, self._outdir)
        return kent.axt_chain(t2bit, q2bit, axtgz)

    def wait_integrate(self, futures: list[confu.Future[Path]]):
        return self.integrate([f.result() for f in futures])

    def integrate(self, chains: list[Path]):
        pre_chain = kent.merge_sort_pre(chains, self._target_sizes, self._query_sizes)
        syntenic_net = kent.chain_net_syntenic(
            pre_chain, self._target_sizes, self._query_sizes
        )
        sing_maf = syntenic_net.parent / "sing.maf"
        kent.net_to_maf(syntenic_net, pre_chain, sing_maf, self._target, self._query)
        if sing_maf.exists():
            print(sing_maf)
        return sing_maf


def lastz(t2bit: Path, q2bit: Path, outdir: Path):
    target_label = t2bit.stem.rsplit("dna_sm.", 1)[1]
    query_label = q2bit.stem.rsplit("dna_sm.", 1)[1]
    subdir = outdir / target_label
    if not cli.dry_run:
        subdir.mkdir(0o755, exist_ok=True)
    axtgz = subdir / f"{query_label}.axt.gz"
    cmd = f"lastz {t2bit} {q2bit} --format=axt"
    cmd += subp.optjoin(config["lastz"])
    is_to_run = fs.is_outdated(axtgz, [t2bit, q2bit])
    lastz = subp.run(cmd, stdout=subp.PIPE, if_=is_to_run)
    if is_to_run and not cli.dry_run:
        with gzip.open(axtgz, "wb") as fout:
            fout.write(lastz.stdout)
    return axtgz


if __name__ == "__main__":
    main()
