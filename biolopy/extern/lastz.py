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
    parser.add_argument("-c", "--clade", choices=phylo.newicks.keys())
    parser.add_argument("target", choices=ensemblgenomes.species_names())
    parser.add_argument("query", nargs="*")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.clade:
        assert not args.query
        run(args.target, args.clade, args.jobs)
    else:
        _run(args.target, args.query, args.jobs)


def run(target: str, clade: str, jobs: int):
    tree = phylo.newicks[clade]
    _run(target, phylo.extract_names(tree), jobs)
    return Path("pairwise") / target


def _run(target: str, queries: list[str], jobs: int):
    queries = ensemblgenomes.sanitize_queries(target, queries)
    _executor._max_workers = jobs
    futures: list[confu.Future[Path]] = []
    for query in queries:
        pa = PairwiseAlignment(target, query)
        futures.extend(pa.run())
    for future in confu.as_completed(futures):
        if (sing_maf := future.result()).exists():
            print(sing_maf)


class PairwiseAlignment:
    def __init__(self, target: str, query: str):
        self._target = target
        self._query = query
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
        flists: list[list[confu.Future[Path]]] = []
        for t in target_chromosomes:
            flists.append(
                [_executor.submit(self.align_chr, t, q) for q in query_chromosomes]
            )
        return [_executor.submit(self.wait_integrate, futures) for futures in flists]

    def align_chr(self, t2bit: Path, q2bit: Path):
        lastz_opts: subp.Optdict = {}
        axtch_opts: subp.Optdict = {}
        axtgz = lastz(t2bit, q2bit, self._outdir, lastz_opts)
        chain = axt_chain(t2bit, q2bit, axtgz, axtch_opts)
        return chain

    def wait_integrate(self, futures: list[confu.Future[Path]]):
        return self.integrate([f.result() for f in futures])

    def integrate(self, chains: list[Path]):
        pre_chain = merge_sort_pre(chains, self._target_sizes, self._query_sizes)
        syntenic_net = chain_net_syntenic(
            pre_chain, self._target_sizes, self._query_sizes
        )
        sing_maf = net_axt_maf(syntenic_net, pre_chain, self._target, self._query)
        return sing_maf


def lastz(t2bit: Path, q2bit: Path, outdir: Path, options: subp.Optdict = {}):
    defaults: subp.Optdict = {
        "gap": None,  # 400,30 (open, extend)
        "xdrop": None,  # 910
        "hspthresh": None,  # 3000
        "ydrop": None,  # 9400 = Open + 30 * Extend
        "gappedthresh": None,  # 3000
        "inner": 2000,
    }
    options = defaults | options
    target_label = t2bit.stem.rsplit("dna_sm.", 1)[1]
    query_label = q2bit.stem.rsplit("dna_sm.", 1)[1]
    subdir = outdir / target_label
    if not cli.dry_run:
        subdir.mkdir(0o755, exist_ok=True)
    axtgz = subdir / f"{query_label}.axt.gz"
    cmd = f"lastz {t2bit} {q2bit} --format=axt"
    cmd += subp.optjoin(options)
    is_to_run = fs.is_outdated(axtgz, [t2bit, q2bit])
    lastz = subp.run_if(is_to_run, cmd, stdout=subp.PIPE)
    if is_to_run and not cli.dry_run:
        with gzip.open(axtgz, "wb") as fout:
            fout.write(lastz.stdout)
    return axtgz


def axt_chain(t2bit: Path, q2bit: Path, axtgz: Path, options: subp.Optdict = {}):
    # medium: mouse/human ~80MYA ~poales/poaceae
    # loose: chicken/human ~300MYA ~gymnosperm/monocot
    defaults: subp.Optdict = {
        "minScore": 3000,
        "linearGap": "medium",
    }
    options = defaults | options
    chain = axtgz.with_suffix("").with_suffix(".chain")
    cmd = "axtChain"
    cmd += subp.optjoin(options, "-")
    cmd += f" stdin {t2bit} {q2bit} {chain}"
    is_to_run = fs.is_outdated(chain, axtgz)
    p = subp.popen_if(is_to_run, cmd, stdin=subp.PIPE)
    if is_to_run and not cli.dry_run:
        assert p.stdin
        with gzip.open(axtgz, "rb") as fin:
            shutil.copyfileobj(fin, p.stdin)
            p.stdin.close()
    p.communicate()
    return chain


def merge_sort_pre(chains: list[Path], target_sizes: Path, query_sizes: Path):
    parent = set(x.parent for x in chains)
    subdir = parent.pop()
    assert not parent, "chains are in the same directory"
    pre_chain = subdir / "pre.chain.gz"
    is_to_run = fs.is_outdated(pre_chain, chains)
    merge_cmd = ["chainMergeSort"] + [str(x) for x in chains]
    merge = subp.popen_if(is_to_run, merge_cmd, stdout=subp.PIPE)
    assert merge.stdout
    pre_cmd = f"chainPreNet stdin {target_sizes} {query_sizes} stdout"
    pre = subp.popen_if(is_to_run, pre_cmd, stdin=merge.stdout, stdout=subp.PIPE)
    merge.stdout.close()
    if is_to_run and not cli.dry_run:
        (stdout, _stderr) = pre.communicate()
        with gzip.open(pre_chain, "wb") as fout:
            fout.write(stdout)
    return pre_chain


def chain_net_syntenic(pre_chain: Path, target_sizes: Path, query_sizes: Path):
    syntenic_net = pre_chain.parent / "syntenic.net"
    is_to_run = fs.is_outdated(syntenic_net, pre_chain)
    cn_cmd = f"chainNet stdin {target_sizes} {query_sizes} stdout /dev/null"
    ns_cmd = f"netSyntenic stdin {syntenic_net}"
    cn = subp.popen_if(is_to_run, cn_cmd, stdin=subp.PIPE, stdout=subp.PIPE)
    ns = subp.popen_if(is_to_run, ns_cmd, stdin=subp.PIPE)
    content = b""
    if is_to_run and not cli.dry_run:
        with gzip.open(pre_chain, "rb") as fout:
            content = fout.read()
    (cn_out, _) = cn.communicate(content)
    ns.communicate(cn_out)
    return syntenic_net


def net_axt_maf(syntenic_net: Path, pre_chain: Path, target: str, query: str):
    sing_maf = syntenic_net.parent / "sing.maf"
    target_2bit = ensemblgenomes.get_file("*.genome.2bit", target)
    query_2bit = ensemblgenomes.get_file("*.genome.2bit", query)
    target_sizes = ensemblgenomes.get_file("fasize.chrom.sizes", target)
    query_sizes = ensemblgenomes.get_file("fasize.chrom.sizes", query)
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
    tprefix = phylo.shorten(target)
    qprefix = phylo.shorten(query)
    axttomaf_cmd = (
        f"axtToMaf -tPrefix={tprefix}. -qPrefix={qprefix}. stdin"
        f" {target_sizes} {query_sizes} {sing_maf}"
    )
    atm = subp.popen_if(is_to_run, axttomaf_cmd, stdin=sort.stdout)
    sort.stdout.close()
    atm.communicate()
    return sing_maf


if __name__ == "__main__":
    main()
