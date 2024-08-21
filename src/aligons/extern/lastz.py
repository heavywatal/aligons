"""Pairwise genome alignment.

src: {db.api.prefix}/{species}/*.fa.gz
dst: ./pairwise/{target}/{query}/{chromosome}/sing.maf

https://lastz.github.io/lastz/
"""

import logging
from pathlib import Path
from typing import Any

from aligons.db import api
from aligons.util import cli, config, fs, maf, subp

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
    futures: list[cli.Future[Path]] = []
    for query in queries:
        pa = PairwiseGenomeAlignment(target, query)
        futures.extend(pa.run())
    cli.wait_raise(futures)
    return Path("pairwise") / target


class PairwiseChromosomeAlignment:
    def __init__(self, bed: Path, queries: list[str]) -> None:
        bed_df = maf.read_bed(bed)
        rows = bed_df.collect().iter_rows(named=True)
        row0 = next(rows)
        self._target: str = row0["name"]
        self._target_sizes = api.fasize(self._target)
        self._t2bit = api.chromosome_2bit(self._target, row0["chrom"])
        self._target_dir = Path("pairwise") / self._target
        self._queries: dict[str, Path] = {}
        for row in rows:
            query = row["name"]
            if query not in queries:
                continue
            q2bit = api.chromosome_2bit(query, row["chrom"])
            self._queries[query] = q2bit

    def submit(self) -> list[cli.Future[Path]]:
        return [
            cli.thread_submit(self.chr_sing_maf, query, q2bit)
            for query, q2bit in self._queries.items()
        ]

    def chr_sing_maf(self, query: str, q2bit: Path) -> Path:
        outdir = self._target_dir / query
        query_sizes = api.fasize(query)
        axt = lastz(self._t2bit, q2bit, outdir)
        chain = kent.axtChain(axt, self._t2bit, q2bit)
        net, _q_net = kent.chain_net(chain, self._target_sizes, query_sizes)
        sing_maf = net.with_name("sing.maf")
        kent.net_to_maf(net, chain, sing_maf, self._target, query)
        return sing_maf

    @property
    def target_dir(self) -> Path:
        return self._target_dir

    @property
    def target(self) -> str:
        return self._target

    @property
    def queries(self) -> list[str]:
        return list(self._queries.keys())


class PairwiseGenomeAlignment:
    def __init__(self, target: str, query: str) -> None:
        self._target = target
        self._query = query
        self._target_sizes = api.fasize(target)
        self._query_sizes = api.fasize(query)
        self._outdir = Path("pairwise") / target / query
        self._lastz_opts = _lastz_options(target, query)

    def run(self) -> list[cli.Future[Path]]:
        pool = cli.ThreadPool()
        target_chromosomes = list(api.iter_chromosome_2bit(self._target))
        query_chromosomes = list(api.iter_chromosome_2bit(self._query))
        ll: list[list[cli.Future[Path]]] = [
            [pool.submit(self.align_chr, t, q) for q in query_chromosomes]
            for t in target_chromosomes
        ]
        return [pool.submit(self.wait_integrate, futures) for futures in ll]

    def align_chr(self, target_2bit: Path, query_2bit: Path) -> Path:
        axt_gz = lastz(target_2bit, query_2bit, self._outdir, **self._lastz_opts)
        return kent.axtChain(axt_gz, target_2bit, query_2bit)

    def wait_integrate(self, futures: list[cli.Future[Path]]) -> Path:
        return self.integrate([f.result() for f in futures])

    def integrate(self, chains: list[Path]) -> Path:
        chain = kent.chainMergeSort(chains)
        net, _q_net = kent.chain_net(chain, self._target_sizes, self._query_sizes)
        if self._target.rsplit("/")[-1] == self._query.rsplit("/")[-1]:
            fs.print_if_exists(kent.net_chain_subset(net, chain))
        sing_maf = net.with_name("sing.maf")
        kent.net_to_maf(net, chain, sing_maf, self._target, self._query)
        return fs.print_if_exists(sing_maf)


def lastz(t2bit: Path, q2bit: Path, outdir: Path, **kwargs: Any) -> Path:
    target_label = t2bit.stem.rsplit("dna_sm.", 1)[1]
    query_label = q2bit.stem.rsplit("dna_sm.", 1)[1]
    subdir = outdir / target_label
    if not cli.dry_run:
        subdir.mkdir(0o755, parents=True, exist_ok=True)
    axt_gz = subdir / f"{query_label}.axt.gz"
    opts = subp.optargs(config["lastz"] | kwargs)
    args = ["lastz", t2bit, q2bit, "--format=axt", *opts]
    is_to_run = fs.is_outdated(axt_gz, [t2bit, q2bit])
    with subp.popen(args, stdout=subp.PIPE, if_=is_to_run) as p:
        return subp.gzip(p.stdout, axt_gz, if_=is_to_run)


def _lastz_options(target: str, query: str) -> dict[str, Any]:
    opts: dict[str, Any] = {}
    if query in config["lastz"].get("close", {}).get(target, []):
        opts["seed"] = "match12"
        opts["match"] = "1,5"
        opts["inner"] = int(config["lastz"].get("inner", 0) / 100)
    return opts


if __name__ == "__main__":
    main()
