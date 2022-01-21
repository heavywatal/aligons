#!/usr/bin/env python3
"""Pairwise genome alignment

src: ensemblgenomes/plants/release-${VERSION}/fasta/*_*
dst: osat_****/chromosome.*/sing.maf

https://lastz.github.io/lastz/
"""
import concurrent.futures as confu
import gzip
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from subprocess import PIPE
from typing import Any, IO

from .db import ensemblgenomes, name


class PairwiseAlignment:
    def __init__(self, target: str, query: str, quick: str, jobs: int, dry_run: bool):
        self._target = target
        self._query = query
        self._quick = quick
        self._jobs = jobs
        self._dry_run = dry_run
        self._target_sizes = ensemblgenomes.single_path("fasize.chrom.sizes", target)
        self._query_sizes = ensemblgenomes.single_path("fasize.chrom.sizes", query)
        self._target_2bit = ensemblgenomes.single_path("*.genome.2bit", target)
        self._query_2bit = ensemblgenomes.single_path("*.genome.2bit", query)
        self._target_nickname = name.shorten(target)
        self._query_nickname = name.shorten(query)
        outdir = f"{self._target_nickname}_{self._query_nickname}"
        print(f"{outdir=}")
        self._outdir = Path(outdir)
        self.run()

    def run(self):
        if not self._dry_run:
            self._outdir.mkdir(0o755, exist_ok=True)
        patt = "*.chromosome.*.2bit"
        target_chromosomes = sorted(ensemblgenomes.rglob(patt, self._target))
        query_chromosomes = sorted(ensemblgenomes.rglob(patt, self._query))
        with confu.ThreadPoolExecutor(max_workers=self._jobs) as executor:
            nested: list[list[confu.Future[Path]]] = []
            for t in target_chromosomes:
                futures = [
                    executor.submit(self.align_chr_pair, t, q)
                    for q in query_chromosomes
                ]
                nested.append(futures)
            subexe = confu.ThreadPoolExecutor(max_workers=None)
            waiters = [subexe.submit(wait_results, fs) for fs in nested]
            for future in confu.as_completed(waiters):
                executor.submit(self.integrate, future.result())
            subexe.shutdown()

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
        patt = r"(?<=dna_sm\.)(.+)$"
        target_label = re_search(patt, target_2bit.stem).group(1)
        query_label = re_search(patt, query_2bit.stem).group(1)
        subdir = self._outdir / f"{target_label}"
        if not self._dry_run:
            subdir.mkdir(0o755, exist_ok=True)
        axtgz = subdir / f"{query_label}.axt.gz"
        args = shlex.split(
            f"lastz {target_2bit} {query_2bit} --format=axt --inner=2000 --step=7"
        )
        if self._quick:
            args.extend(["--notransition", "--nogapped"])
        lastz = self.popen_if(not axtgz.exists(), args, stdout=PIPE)
        if not axtgz.exists() and not self._dry_run:
            assert lastz.stdout
            with gzip.open(axtgz, "wb") as fout:
                shutil.copyfileobj(lastz.stdout, fout)
        return axtgz

    def axt_chain(self, target_2bit: Path, query_2bit: Path, axtgz: Path):
        chain = axtgz.with_suffix("").with_suffix(".chain")
        args = shlex.split(
            "axtChain -minScore=5000 -linearGap=medium"
            f" stdin {target_2bit} {query_2bit} {chain}"
        )
        p = self.popen_if(not chain.exists(), args, stdin=PIPE)
        if not chain.exists() and not self._dry_run:
            assert p.stdin
            with gzip.open(axtgz, "rb") as fout:
                shutil.copyfileobj(fout, p.stdin)
                p.stdin.close()
        p.communicate()
        return chain

    def merge_sort_pre(self, chains: list[Path]):
        parent = set(x.parent for x in chains)
        subdir = parent.pop()
        assert not parent, "chains are in the same directory"
        pre_chain = subdir / "pre.chain.gz"
        args = ["chainMergeSort"] + [str(x) for x in chains]
        merge = self.popen_if(not pre_chain.exists(), args, stdout=PIPE)
        args = shlex.split(
            f"chainPreNet stdin {self._target_sizes} {self._query_sizes} stdout"
        )
        pre = self.popen_if(
            not pre_chain.exists(), args, stdin=merge.stdout, stdout=PIPE
        )
        if not pre_chain.exists() and not self._dry_run:
            assert pre.stdout
            with gzip.open(pre_chain, "wb") as fout:
                shutil.copyfileobj(pre.stdout, fout)
        return pre_chain

    def chain_net_syntenic(self, pre_chain: Path):
        syntenic_net = pre_chain.parent / "syntenic.net"
        args = shlex.split(
            f"chainNet stdin {self._target_sizes} {self._query_sizes} stdout /dev/null"
        )
        cn = self.popen_if(not syntenic_net.exists(), args, stdin=PIPE, stdout=PIPE)
        if not syntenic_net.exists() and not self._dry_run:
            assert cn.stdin
            with gzip.open(pre_chain, "rb") as fout:
                shutil.copyfileobj(fout, cn.stdin)
                cn.stdin.close()
        args = shlex.split(f"netSyntenic stdin {syntenic_net}")
        self.run_if(not syntenic_net.exists(), args, stdin=cn.stdout)
        return syntenic_net

    def net_axt_maf(self, syntenic_net: Path, pre_chain: Path):
        sing_maf = syntenic_net.parent / "sing.maf"
        args = shlex.split(
            f"netToAxt {syntenic_net} stdin {self._target_2bit} {self._query_2bit} stdout"
        )
        toaxt: subprocess.Popen[bytes] = self.popen_if(
            not sing_maf.exists(), args, stdin=PIPE, stdout=PIPE
        )
        if not sing_maf.exists() and not self._dry_run:
            assert toaxt.stdin
            with gzip.open(pre_chain, "rb") as fout:
                shutil.copyfileobj(fout, toaxt.stdin)
                toaxt.stdin.close()
        args = shlex.split("axtSort stdin stdout")
        sort = self.popen_if(
            not sing_maf.exists(), args, stdin=toaxt.stdout, stdout=PIPE
        )
        args = shlex.split(
            "axtToMaf"
            f" -tPrefix={self._target_nickname}."
            f" -qPrefix={self._query_nickname}. stdin"
            f" {self._target_sizes} {self._query_sizes} {sing_maf}"
        )
        self.run_if(not sing_maf.exists(), args, stdin=sort.stdout)
        if not sing_maf.exists():
            print(f"## {sing_maf} does not exist.")
        return sing_maf

    def popen_if(
        self,
        cond: bool,
        args: list[str],
        stdin: IO[bytes] | int | None = None,
        stdout: IO[bytes] | int | None = None,
    ):  # kwargs hinders type inference to Popen[bytes]
        if cond and not self._dry_run:
            print(" ".join(args))
            return subprocess.Popen(args, stdin=stdin, stdout=stdout)
        print("# " + " ".join(args))
        return subprocess.Popen(["sleep", "0"], stdin=stdin, stdout=stdout)

    def run_if(self, cond: bool, args: list[str], **kwargs: Any):
        if cond and not self._dry_run:
            print(" ".join(args))
            return subprocess.run(args, **kwargs)
        print("# " + " ".join(args))
        return subprocess.run(["sleep", "0"], **kwargs)


def wait_results(futures: list[confu.Future[Any]]):
    return [f.result() for f in futures]


def re_search(pattern: str, string: str):
    mobj = re.search(pattern, string)
    assert mobj  # shut up None warning
    return mobj


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--dry_run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("target", default=os.getenv("TARGET"))
    parser.add_argument("query", default=os.getenv("QUERY"))
    args = parser.parse_args()

    targets = [x.strip() for x in args.target.split(",") if x]
    queries = [x.strip() for x in args.query.split(",") if x]
    print(f"{targets=}")
    print(f"{queries=}")
    for target in targets:
        for query in queries:
            if target == query:
                continue
            print(target, query)
            PairwiseAlignment(
                target, query, quick=args.quick, jobs=args.jobs, dry_run=args.dry_run
            )


if __name__ == "__main__":
    main()
