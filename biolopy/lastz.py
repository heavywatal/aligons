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
import sys
from pathlib import Path
from subprocess import PIPE
from typing import Any

from .db import name


def unlist(x: list[Path]):
    assert len(x) == 1
    return str(x[0])


class PairwiseAlignment:
    version = ""
    prefix_fasta = Path("")

    @classmethod
    def set_version(cls, version: str):
        cls.version = version
        x = f"~/db/ensemblgenomes/plants/release-{cls.version}"
        cls.prefix_fasta = Path(x).expanduser() / "fasta"
        print(f"Using v{cls.version}: {cls.prefix_fasta}/", file=sys.stderr)

    @classmethod
    def ls(cls):
        if cls.prefix_fasta is None:
            cls.set_version(cls.version)
        return sorted([x for x in cls.prefix_fasta.iterdir() if x.is_dir()])

    def __init__(self, target: str, query: str, quick: str, jobs: int, dry_run: bool):
        self._target = target
        self._query = query
        self._quick = quick
        self._jobs = jobs
        self._dry_run = dry_run
        pf = self.prefix_fasta
        self._target_sizes = pf / target / "dna" / "fasize.chrom.sizes"
        self._query_sizes = pf / query / "dna" / "fasize.chrom.sizes"
        self._target_2bit = unlist(self.rglob(target, "*.genome.2bit"))
        self._query_2bit = unlist(self.rglob(query, "*.genome.2bit"))
        assert self._target_sizes.exists()
        assert self._query_sizes.exists()
        assert Path(self._target_2bit).exists()
        assert Path(self._query_2bit).exists()
        self._target_nickname = name.shorten(target)
        self._query_nickname = name.shorten(query)
        outdir = f"{self._target_nickname}_{self._query_nickname}"
        print(f"{outdir=}")
        self._outdir = Path(outdir)
        self.run()

    def run(self):
        if not self._dry_run:
            self._outdir.mkdir(0o755, exist_ok=True)
        target_chromosomes = self.rglob(self._target, "*.chromosome.*.2bit")
        query_chromosomes = self.rglob(self._query, "*.chromosome.*.2bit")
        assert target_chromosomes
        assert query_chromosomes
        chains: list[str] = []
        with confu.ThreadPoolExecutor(max_workers=self._jobs) as executor:
            futures = [
                executor.submit(self.lastz_chain, t, q)
                for t in target_chromosomes
                for q in query_chromosomes
            ]
            (_done, _notdone) = confu.wait(futures)
            for future in futures:
                chains.append(future.result())
        subdirs = set(Path(x).parent for x in chains)
        with confu.ThreadPoolExecutor(max_workers=self._jobs) as executor:
            futures = [executor.submit(self.sub, x) for x in subdirs]
            (_done, _notdone) = confu.wait(futures)
            for future in futures:
                future.result()

    def rglob(self, directory: str, pattern: str):
        path = PairwiseAlignment.prefix_fasta / directory
        return sorted(list(path.rglob(pattern)))

    def sub(self, subdir: Path):
        pre_chain = subdir / "pre.chain.gz"
        syntenic_net = subdir / "syntenic.net"
        sing_maf = subdir / "sing.maf"
        if not pre_chain.exists():
            self.merge_sort_pre(subdir, pre_chain)
        if not syntenic_net.exists():
            self.chain_net_syntenic(pre_chain, syntenic_net)
        if not sing_maf.exists():
            self.net_axt_maf(syntenic_net, pre_chain, sing_maf)

    def lastz_chain(self, target_2bit: Path, query_2bit: Path, force: bool = False):
        pattern = r"(?<=dna_sm.)(.+?)\.2bit$"
        mobj = re.search(pattern, str(target_2bit))
        assert mobj
        target_label = mobj.group(1)
        mobj = re.search(pattern, str(query_2bit))
        assert mobj
        query_label = mobj.group(1)
        subdir = self._outdir / f"{target_label}"
        if not self._dry_run:
            subdir.mkdir(0o755, exist_ok=True)
        axtgz = subdir / f"{query_label}.axt.gz"
        if not axtgz.exists() or force:
            args = f"lastz {target_2bit} {query_2bit} --format=axt"
            args += " --inner=2000 --step=7"
            if self._quick:
                args += " --notransition --nogapped"
            print(args)
            lastz = self.popen(shlex.split(args), stdout=PIPE)
            assert lastz.stdout
            with gzip.open(axtgz, "wt") as fout:
                shutil.copyfileobj(lastz.stdout, fout)

        chain = re.sub(r"axt\.gz$", "chain", str(axtgz))
        if not Path(chain).exists() or force:
            args = (
                "axtChain -minScore=5000 -linearGap=medium"
                f" stdin {target_2bit} {query_2bit} {chain}"
            )
            print(args)
            if not self._dry_run:
                gz = subprocess.Popen(["gunzip", "-c", axtgz], stdout=PIPE)
                subprocess.run(shlex.split(args), stdin=gz.stdout)
        return chain

    def merge_sort_pre(self, subdir: Path, outfile: Path):
        chains = sorted([str(x) for x in subdir.glob("*.chain")])
        args = ["chainMergeSort"] + chains
        print(" ".join(args))
        assert chains
        merge = self.popen(args, stdout=PIPE)
        args = "chainPreNet stdin" f" {self._target_sizes} {self._query_sizes} stdout"
        print(args)
        pre = self.popen(shlex.split(args), stdin=merge.stdout, stdout=PIPE)
        assert pre.stdout
        with gzip.open(outfile, "wt") as fout:
            shutil.copyfileobj(pre.stdout, fout)

    def chain_net_syntenic(self, all_chain: Path, outfile: Path):
        gz = self.popen(["gunzip", "-c", str(all_chain)], stdout=PIPE)
        args = (
            f"chainNet stdin {self._target_sizes} {self._query_sizes}"
            " stdout /dev/null"
        )
        print(args)
        net = self.popen(shlex.split(args), stdin=gz.stdout, stdout=PIPE)
        print(all_chain)
        args = f"netSyntenic stdin {outfile}"
        print(args)
        if not self._dry_run:
            subprocess.run(shlex.split(args), stdin=net.stdout)

    def net_axt_maf(self, syntenic_net: Path, all_chain: Path, outfile: Path):
        gz = self.popen(["gunzip", "-c", str(all_chain)], stdout=PIPE)
        args = (
            f"netToAxt {syntenic_net} stdin"
            f" {self._target_2bit} {self._query_2bit} stdout"
        )
        print(args)
        toaxt = self.popen(shlex.split(args), stdin=gz.stdout, stdout=PIPE)
        args = "axtSort stdin stdout"
        print(args)
        sort = self.popen(shlex.split(args), stdin=toaxt.stdout, stdout=PIPE)
        args = (
            "axtToMaf"
            f" -tPrefix={self._target_nickname}."
            f" -qPrefix={self._query_nickname}. stdin"
            f" {self._target_sizes} {self._query_sizes} {outfile}"
        )
        print(args)
        if not self._dry_run:
            subprocess.run(shlex.split(args), stdin=sort.stdout)

    def popen(self, args: list[str], *pargs: Any, **kwargs: Any):
        if self._dry_run:
            return subprocess.Popen([":"], *pargs, **kwargs)
        return subprocess.Popen(args, *pargs, **kwargs)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--dry_run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("-v", "--version", default=os.environ["ENSEMBLGENOMES_VERSION"])
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("-t", "--target", default=os.getenv("TARGET"))
    parser.add_argument("-q", "--query", default=os.getenv("QUERY"))
    args = parser.parse_args()

    PairwiseAlignment.set_version(args.version)
    targets = filter(None, (x.strip() for x in args.target.split(",")))
    queries = filter(None, (x.strip() for x in args.query.split(",")))

    print(f"{targets=}")
    print(f"{queries=}")
    i = 0
    for target in targets:
        for query in queries:
            if target == query:
                continue
            print(target, query)
            PairwiseAlignment(
                target, query, quick=args.quick, jobs=args.jobs, dry_run=args.dry_run
            )
            i += 1


if __name__ == "__main__":
    main()
