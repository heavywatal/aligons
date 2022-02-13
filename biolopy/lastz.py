"""Pairwise genome alignment

src: {ensemblgenomes.prefix}/fasta/{species}/*.fa.gz
dst: ./pairwise/{target}/{query}/{chromosome}/sing.maf

https://lastz.github.io/lastz/
"""
import argparse
import concurrent.futures as confu
import gzip
import logging
import os
import shutil
from pathlib import Path
from subprocess import PIPE
from typing import Any

from . import cli, fs
from .db import ensemblgenomes, name

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("target")
    parser.add_argument("query")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    assert args.target in ensemblgenomes.list_species()
    assert args.query in ensemblgenomes.list_species()
    PairwiseAlignment(args.target, args.query, quick=args.quick, jobs=args.jobs)


class PairwiseAlignment:
    def __init__(self, target: str, query: str, quick: str, jobs: int):
        self._target = target
        self._query = query
        self._quick = quick
        self._jobs = jobs
        self._target_sizes = ensemblgenomes.get_file("fasize.chrom.sizes", target)
        self._query_sizes = ensemblgenomes.get_file("fasize.chrom.sizes", query)
        self._outdir = Path(f"pairwise/{self._target}/{self._query}")
        self.run()

    def run(self):
        if not cli.dry_run:
            self._outdir.mkdir(0o755, parents=True, exist_ok=True)
        patt = "*.chromosome.*.2bit"
        target_chromosomes = ensemblgenomes.rglob(patt, self._target)
        query_chromosomes = list(ensemblgenomes.rglob(patt, self._query))
        with confu.ThreadPoolExecutor(max_workers=self._jobs) as executor:
            nested: list[list[confu.Future[Path]]] = []
            for t in fs.sorted_naturally(target_chromosomes):
                futures = [
                    executor.submit(self.align_chr_pair, t, q)
                    for q in fs.sorted_naturally(query_chromosomes)
                ]
                nested.append(futures)
            subexe = confu.ThreadPoolExecutor(max_workers=None)
            waiters = [subexe.submit(wait_results, fts) for fts in nested]
            futures = [
                executor.submit(self.integrate, future.result())
                for future in confu.as_completed(waiters)
            ]
            for future in confu.as_completed(futures):
                product = future.result()
                if product.exists():
                    print(product)

    def align_chr_pair(self, target_2bit: Path, query_2bit: Path):
        axtgz = self.lastz(target_2bit, query_2bit)
        chain = self.axt_chain(target_2bit, query_2bit, axtgz)
        return chain

    def integrate(self, chains: list[Path]):
        pre_chain = self.merge_sort_pre(chains)
        syntenic_net = self.chain_net_syntenic(pre_chain)
        sing_maf = self.net_axt_maf(syntenic_net, pre_chain)
        return sing_maf

    def lastz(self, target_2bit: Path, query_2bit: Path):
        target_label = target_2bit.stem.rsplit("dna_sm.", 1)[1]
        query_label = query_2bit.stem.rsplit("dna_sm.", 1)[1]
        subdir = self._outdir / target_label
        if not cli.dry_run:
            subdir.mkdir(0o755, exist_ok=True)
        axtgz = subdir / f"{query_label}.axt.gz"
        args = f"lastz {target_2bit} {query_2bit} --format=axt --inner=2000 --step=7"
        if self._quick:
            args += " --notransition --nogapped"
        lastz = cli.popen_if(not axtgz.exists(), args, stdout=PIPE)
        if not axtgz.exists() and not cli.dry_run:
            assert lastz.stdout
            with gzip.open(axtgz, "wb") as fout:
                shutil.copyfileobj(lastz.stdout, fout)
        return axtgz

    def axt_chain(self, target_2bit: Path, query_2bit: Path, axtgz: Path):
        chain = axtgz.with_suffix("").with_suffix(".chain")
        cmd = "axtChain -minScore=5000 -linearGap=medium stdin"
        cmd += f" {target_2bit} {query_2bit} {chain}"
        p = cli.popen_if(not chain.exists(), cmd, stdin=PIPE)
        if not chain.exists() and not cli.dry_run:
            assert p.stdin
            with gzip.open(axtgz, "rb") as fin:
                shutil.copyfileobj(fin, p.stdin)
                p.stdin.close()
        p.communicate()
        return chain

    def merge_sort_pre(self, chains: list[Path]):
        parent = set(x.parent for x in chains)
        subdir = parent.pop()
        assert not parent, "chains are in the same directory"
        pre_chain = subdir / "pre.chain.gz"
        merge = cli.popen_if(
            not pre_chain.exists(),
            ["chainMergeSort"] + [str(x) for x in chains],
            stdout=PIPE,
        )
        pre = cli.popen_if(
            not pre_chain.exists(),
            f"chainPreNet stdin {self._target_sizes} {self._query_sizes} stdout",
            stdin=merge.stdout,
            stdout=PIPE,
        )
        if not pre_chain.exists() and not cli.dry_run:
            assert pre.stdout
            with gzip.open(pre_chain, "wb") as fout:
                shutil.copyfileobj(pre.stdout, fout)
        return pre_chain

    def chain_net_syntenic(self, pre_chain: Path):
        syntenic_net = pre_chain.parent / "syntenic.net"
        cn = cli.popen_if(
            not syntenic_net.exists(),
            f"chainNet stdin {self._target_sizes} {self._query_sizes} stdout /dev/null",
            stdin=PIPE,
            stdout=PIPE,
        )
        if not syntenic_net.exists() and not cli.dry_run:
            assert cn.stdin
            with gzip.open(pre_chain, "rb") as fout:
                shutil.copyfileobj(fout, cn.stdin)
                cn.stdin.close()
        cli.run_if(
            not syntenic_net.exists(),
            f"netSyntenic stdin {syntenic_net}",
            stdin=cn.stdout,
        )
        return syntenic_net

    def net_axt_maf(self, syntenic_net: Path, pre_chain: Path):
        sing_maf = syntenic_net.parent / "sing.maf"
        target_2bit = ensemblgenomes.get_file("*.genome.2bit", self._target)
        query_2bit = ensemblgenomes.get_file("*.genome.2bit", self._query)
        toaxt = cli.popen_if(
            not sing_maf.exists(),
            f"netToAxt {syntenic_net} stdin {target_2bit} {query_2bit} stdout",
            stdin=PIPE,
            stdout=PIPE,
        )
        if not sing_maf.exists() and not cli.dry_run:
            assert toaxt.stdin
            with gzip.open(pre_chain, "rb") as fout:
                shutil.copyfileobj(fout, toaxt.stdin)
                toaxt.stdin.close()
        sort = cli.popen_if(
            not sing_maf.exists(),
            "axtSort stdin stdout",
            stdin=toaxt.stdout,
            stdout=PIPE,
        )
        tprefix = name.shorten(self._target)
        qprefix = name.shorten(self._query)
        args = (
            f"axtToMaf -tPrefix={tprefix}. -qPrefix={qprefix}. stdin"
            f" {self._target_sizes} {self._query_sizes} {sing_maf}"
        )
        cli.run_if(not sing_maf.exists(), args, stdin=sort.stdout)
        return sing_maf


def wait_results(futures: list[confu.Future[Any]]):
    return [f.result() for f in futures]


if __name__ == "__main__":
    main()
