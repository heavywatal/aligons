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
from typing import Any

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
        print(f"{subdirs=}")
        with confu.ThreadPoolExecutor(max_workers=self._jobs) as executor:
            futures = [executor.submit(self.sub, x) for x in subdirs]
            (_done, _notdone) = confu.wait(futures)
            for future in futures:
                future.result()

    def rglob(self, directory: str, pattern: str):
        return sorted(ensemblgenomes.rglob(pattern, directory))

    def sub(self, subdir: Path):
        pre_chain = subdir / "pre.chain.gz"
        syntenic_net = subdir / "syntenic.net"
        sing_maf = subdir / "sing.maf"
        self.merge_sort_pre(subdir, pre_chain)
        self.chain_net_syntenic(pre_chain, syntenic_net)
        self.net_axt_maf(syntenic_net, pre_chain, sing_maf)

    def lastz_chain(self, target_2bit: Path, query_2bit: Path):
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
        args = f"lastz {target_2bit} {query_2bit} --format=axt"
        args += " --inner=2000 --step=7"
        if self._quick:
            args += " --notransition --nogapped"
        lastz = self.popen_if(not axtgz.exists(), shlex.split(args), stdout=PIPE)
        assert lastz.stdout
        if not self._dry_run:
            with gzip.open(axtgz, "wt") as fout:
                shutil.copyfileobj(lastz.stdout, fout)

        chain = re.sub(r"axt\.gz$", "chain", str(axtgz))
        args = (
            "axtChain -minScore=5000 -linearGap=medium"
            f" stdin {target_2bit} {query_2bit} {chain}"
        )
        gz = self.popen_if(
            not Path(chain).exists(), ["gunzip", "-c", str(axtgz)], stdout=PIPE
        )
        self.run_if(not Path(chain).exists(), shlex.split(args), stdin=gz.stdout)
        return chain

    def merge_sort_pre(self, subdir: Path, outfile: Path):
        chains = sorted([str(x) for x in subdir.glob("*.chain")])
        args = ["chainMergeSort"] + chains
        assert chains
        merge = self.popen_if(not outfile.exists(), args, stdout=PIPE)
        args = f"chainPreNet stdin {self._target_sizes} {self._query_sizes} stdout"
        pre = self.popen_if(
            not outfile.exists(), shlex.split(args), stdin=merge.stdout, stdout=PIPE
        )
        assert pre.stdout
        if not outfile.exists() and not self._dry_run:
            with gzip.open(outfile, "wt") as fout:
                shutil.copyfileobj(pre.stdout, fout)

    def chain_net_syntenic(self, all_chain: Path, outfile: Path):
        gz = self.popen_if(
            not outfile.exists(), ["gunzip", "-c", str(all_chain)], stdout=PIPE
        )
        args = (
            f"chainNet stdin {self._target_sizes} {self._query_sizes} stdout /dev/null"
        )
        net = self.popen_if(
            not outfile.exists(), shlex.split(args), stdin=gz.stdout, stdout=PIPE
        )
        args = f"netSyntenic stdin {outfile}"
        self.run_if(not outfile.exists(), shlex.split(args), stdin=net.stdout)

    def net_axt_maf(self, syntenic_net: Path, all_chain: Path, outfile: Path):
        gz = self.popen_if(
            not outfile.exists(), ["gunzip", "-c", str(all_chain)], stdout=PIPE
        )
        args = f"netToAxt {syntenic_net} stdin {self._target_2bit} {self._query_2bit} stdout"
        toaxt = self.popen_if(
            not outfile.exists(), shlex.split(args), stdin=gz.stdout, stdout=PIPE
        )
        args = "axtSort stdin stdout"
        sort = self.popen_if(
            not outfile.exists(), shlex.split(args), stdin=toaxt.stdout, stdout=PIPE
        )
        args = (
            "axtToMaf"
            f" -tPrefix={self._target_nickname}."
            f" -qPrefix={self._query_nickname}. stdin"
            f" {self._target_sizes} {self._query_sizes} {outfile}"
        )
        self.run_if(not outfile.exists(), shlex.split(args), stdin=sort.stdout)

    def popen_if(self, cond: bool, args: list[str], *pargs: Any, **kwargs: Any):
        if cond and not self._dry_run:
            print(" ".join(args))
            return subprocess.Popen(args, *pargs, **kwargs)
        print("# " + " ".join(args))
        return subprocess.Popen([":"], shell=True, *pargs, **kwargs)

    def run_if(self, cond: bool, args: list[str], *pargs: Any, **kwargs: Any):
        if cond and not self._dry_run:
            print(" ".join(args))
            return subprocess.run(args, *pargs, **kwargs)
        print("# " + " ".join(args))
        return subprocess.run([":"], shell=True, *pargs, **kwargs)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--dry_run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("-t", "--target", default=os.getenv("TARGET"))
    parser.add_argument("-q", "--query", default=os.getenv("QUERY"))
    args = parser.parse_args()

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
