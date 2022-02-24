"""Pairwise genome alignment

src: {ensemblgenomes.prefix}/fasta/{species}/*.fa.gz
dst: ./pairwise/{target}/{query}/{chromosome}/sing.maf

https://lastz.github.io/lastz/
"""
import concurrent.futures as confu
import gzip
import logging
import os
import shutil
from pathlib import Path

from ..db import ensemblgenomes, phylo
from ..util import cli, fs, subp

_log = logging.getLogger(__name__)
_executor = confu.ThreadPoolExecutor()


def main(argv: list[str] = []):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("-c", "--clade", choices=phylo.trees.keys())
    parser.add_argument("target", choices=ensemblgenomes.species_names())
    parser.add_argument("query", nargs="*")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.clade:
        assert not args.query
        run(args.target, args.clade, args.jobs, args.quick)
    else:
        _run(args.target, args.query, args.jobs, args.quick)


def run(target: str, clade: str, jobs: int, quick: bool = False):
    tree = phylo.trees[clade]
    _run(target, phylo.extract_labels(tree), jobs, quick)
    return Path("pairwise") / target


def _run(target: str, queries: list[str], jobs: int, quick: bool):
    queries = ensemblgenomes.sanitize_queries(target, queries)
    _executor._max_workers = jobs
    futures: list[confu.Future[Path]] = []
    for query in queries:
        pa = PairwiseAlignment(target, query, quick=quick)
        futures.extend(pa.run())
    for future in confu.as_completed(futures):
        if (sing_maf := future.result()).exists():
            print(sing_maf)


class PairwiseAlignment:
    def __init__(self, target: str, query: str, quick: bool):
        self._target = target
        self._query = query
        self._quick = quick
        self._target_sizes = ensemblgenomes.get_file("fasize.chrom.sizes", target)
        self._query_sizes = ensemblgenomes.get_file("fasize.chrom.sizes", query)
        self._outdir = Path("pairwise") / target / query

    def run(self):
        if not cli.dry_run:
            self._outdir.mkdir(0o755, parents=True, exist_ok=True)
        patt = "*.chromosome.*.2bit"
        it = ensemblgenomes.rglob(patt, [self._target])
        target_chromosomes = fs.sorted_naturally(it)
        it = ensemblgenomes.rglob(patt, [self._query])
        query_chromosomes = fs.sorted_naturally(it)
        subexe = confu.ThreadPoolExecutor(max_workers=len(target_chromosomes))
        waiters: list[confu.Future[list[Path]]] = []
        for t in target_chromosomes:
            futures = [
                _executor.submit(self.align_chr, t, q) for q in query_chromosomes
            ]
            waiters.append(subexe.submit(wait_results, futures))
        return [
            _executor.submit(self.integrate, future.result())
            for future in confu.as_completed(waiters)
        ]

    def align_chr(self, target_2bit: Path, query_2bit: Path):
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
        is_to_run = fs.is_outdated(axtgz, [target_2bit, query_2bit])
        lastz = subp.run_if(is_to_run, args, stdout=subp.PIPE)
        if is_to_run and not cli.dry_run:
            with gzip.open(axtgz, "wb") as fout:
                fout.write(lastz.stdout)
        return axtgz

    def axt_chain(self, target_2bit: Path, query_2bit: Path, axtgz: Path):
        chain = axtgz.with_suffix("").with_suffix(".chain")
        cmd = "axtChain -minScore=5000 -linearGap=medium stdin"
        cmd += f" {target_2bit} {query_2bit} {chain}"
        is_to_run = fs.is_outdated(chain, axtgz)
        p = subp.popen_if(is_to_run, cmd, stdin=subp.PIPE)
        if is_to_run and not cli.dry_run:
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
        is_to_run = fs.is_outdated(pre_chain, chains)
        merge_cmd = ["chainMergeSort"] + [str(x) for x in chains]
        merge = subp.popen_if(is_to_run, merge_cmd, stdout=subp.PIPE)
        assert merge.stdout
        pre_cmd = f"chainPreNet stdin {self._target_sizes} {self._query_sizes} stdout"
        pre = subp.popen_if(is_to_run, pre_cmd, stdin=merge.stdout, stdout=subp.PIPE)
        merge.stdout.close()
        if is_to_run and not cli.dry_run:
            (stdout, _stderr) = pre.communicate()
            with gzip.open(pre_chain, "wb") as fout:
                fout.write(stdout)
        return pre_chain

    def chain_net_syntenic(self, pre_chain: Path):
        syntenic_net = pre_chain.parent / "syntenic.net"
        is_to_run = fs.is_outdated(syntenic_net, pre_chain)
        cn_cmd = (
            f"chainNet stdin {self._target_sizes} {self._query_sizes} stdout /dev/null"
        )
        cn = subp.popen_if(is_to_run, cn_cmd, stdin=subp.PIPE, stdout=subp.PIPE)
        assert cn.stdin
        assert cn.stdout
        if is_to_run and not cli.dry_run:
            with gzip.open(pre_chain, "rb") as fout:
                shutil.copyfileobj(fout, cn.stdin)
                cn.stdin.close()
        sn = subp.popen_if(
            is_to_run, f"netSyntenic stdin {syntenic_net}", stdin=cn.stdout
        )
        cn.stdout.close()
        sn.communicate()
        return syntenic_net

    def net_axt_maf(self, syntenic_net: Path, pre_chain: Path):
        sing_maf = syntenic_net.parent / "sing.maf"
        target_2bit = ensemblgenomes.get_file("*.genome.2bit", self._target)
        query_2bit = ensemblgenomes.get_file("*.genome.2bit", self._query)
        is_to_run = fs.is_outdated(sing_maf, [syntenic_net, pre_chain])
        toaxt_cmd = f"netToAxt {syntenic_net} stdin {target_2bit} {query_2bit} stdout"
        toaxt = subp.popen_if(is_to_run, toaxt_cmd, stdin=subp.PIPE, stdout=subp.PIPE)
        assert toaxt.stdin
        assert toaxt.stdout
        if is_to_run and not cli.dry_run:
            with gzip.open(pre_chain, "rb") as fout:
                shutil.copyfileobj(fout, toaxt.stdin)
                toaxt.stdin.close()
        sort = subp.popen_if(
            is_to_run, "axtSort stdin stdout", stdin=toaxt.stdout, stdout=subp.PIPE
        )
        toaxt.stdout.close()
        assert sort.stdout
        tprefix = phylo.shorten(self._target)
        qprefix = phylo.shorten(self._query)
        axttomaf_cmd = (
            f"axtToMaf -tPrefix={tprefix}. -qPrefix={qprefix}. stdin"
            f" {self._target_sizes} {self._query_sizes} {sing_maf}"
        )
        atm = subp.popen_if(is_to_run, axttomaf_cmd, stdin=sort.stdout)
        sort.stdout.close()
        atm.communicate()
        return sing_maf


def wait_results(futures: list[confu.Future[Path]]):
    return [f.result() for f in futures]


if __name__ == "__main__":
    main()
