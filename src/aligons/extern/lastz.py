"""Pairwise genome alignment.

src: {db.api.prefix}/fasta/{species}/*.fa.gz
dst: ./pairwise/{target}/{query}/{chromosome}/sing.maf

https://lastz.github.io/lastz/
"""
import concurrent.futures as confu
import logging
from collections.abc import Iterator
from pathlib import Path

from aligons.db import api
from aligons.util import cli, config, fs, subp

from . import kent

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("target", choices=api.species_names())
    parser.add_argument("query", nargs="*")
    args = parser.parse_args(argv or None)
    run(args.target, args.query)


def run(target: str, queries: list[str]) -> Path:
    queries = api.sanitize_queries(target, queries)
    futures: list[confu.Future[Path]] = []
    for query in queries:
        pa = PairwiseAlignment(target, query)
        futures.extend(pa.run())
    cli.wait_raise(futures)
    return Path("pairwise") / target


class PairwiseAlignment:
    def __init__(self, target: str, query: str) -> None:
        self._target = target
        self._query = query
        self._target_sizes = api.fasize(target)
        self._query_sizes = api.fasize(query)
        self._outdir = Path("pairwise") / target / query

    def run(self) -> list[cli.FuturePath]:
        pool = cli.ThreadPool()
        if not cli.dry_run:
            self._outdir.mkdir(0o755, parents=True, exist_ok=True)
        target_chromosomes = list(iter_chromosome_2bit(self._target))
        query_chromosomes = list(iter_chromosome_2bit(self._query))
        flists: list[list[confu.Future[Path]]] = [
            [pool.submit(self.align_chr, t, q) for q in query_chromosomes]
            for t in target_chromosomes
        ]
        return [pool.submit(self.wait_integrate, futures) for futures in flists]

    def align_chr(self, target_2bit: Path, query_2bit: Path) -> Path:
        axtgz = lastz(target_2bit, query_2bit, self._outdir)
        return kent.axtChain(axtgz, target_2bit, query_2bit)

    def wait_integrate(self, futures: list[confu.Future[Path]]) -> Path:
        return self.integrate([f.result() for f in futures])

    def integrate(self, chains: list[Path]) -> Path:
        chain = kent.chainMergeSort(chains)
        net, _qnet = kent.chain_net(chain, self._target_sizes, self._query_sizes)
        sing_maf = net.with_name("sing.maf")
        kent.net_to_maf(net, chain, sing_maf, self._target, self._query)
        if sing_maf.exists():
            print(sing_maf)
        return sing_maf


def lastz(t2bit: Path, q2bit: Path, outdir: Path) -> Path:
    target_label = t2bit.stem.rsplit("dna_sm.", 1)[1]
    query_label = q2bit.stem.rsplit("dna_sm.", 1)[1]
    subdir = outdir / target_label
    if not cli.dry_run:
        subdir.mkdir(0o755, exist_ok=True)
    axtgz = subdir / f"{query_label}.axt.gz"
    args = ["lastz", t2bit, q2bit, "--format=axt", *subp.optargs(config["lastz"])]
    is_to_run = fs.is_outdated(axtgz, [t2bit, q2bit])
    with subp.popen(args, stdout=subp.PIPE, if_=is_to_run) as p:
        return subp.gzip(p.stdout, axtgz, if_=is_to_run)


def iter_chromosome_2bit(species: str) -> Iterator[Path]:
    fts = [
        cli.thread_submit(kent.faToTwoBit, f) for f in api.list_chromosome_fa(species)
    ]
    for future in fts:
        yield future.result()


if __name__ == "__main__":
    main()
