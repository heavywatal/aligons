"""Integrate chromosome wigs into a genome-wide bigwig.

src: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz
dst: ./multiple/{target}/{clade}/phastcons.bw

https://github.com/ucscGenomeBrowser/kent
"""
import functools
import logging
import os
import re
from pathlib import Path

from aligons.db import api
from aligons.util import cli, config, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("clade", type=Path)
    args = parser.parse_args(argv or None)
    run(args.clade)


def run(clade: Path) -> None:
    if (bigwig := integrate_wigs(clade)).exists():
        print(bigwig)
        _log.info(bigWigInfo(bigwig).rstrip())


def integrate_wigs(clade: Path) -> Path:
    species = clade.parent.name
    chrom_sizes = api.fasize(species)
    name = "phastcons.wig.gz"
    wigs = [p / name for p in fs.sorted_naturally(clade.glob("chromosome.*"))]
    _log.debug(f"{[str(x) for x in wigs]}")
    outfile = clade / "phastcons.bw"
    is_to_run = not cli.dry_run and fs.is_outdated(outfile, wigs)
    args = ["wigToBigWig", "stdin", chrom_sizes, outfile]
    with subp.popen(args, stdin=subp.PIPE, if_=is_to_run) as w2bw:
        for wig in wigs:
            subp.run_zcat(wig, w2bw.stdin, if_=is_to_run)
    return outfile


def bigWigInfo(path: Path) -> bytes:  # noqa: N802
    args = ["bigWigInfo", path]
    return subp.run(args, stdout=subp.PIPE, text=True).stdout


def bedGraphToBigWig(infile: Path, chrom_sizes: Path) -> Path:  # noqa: N802
    name, cnt = re.subn(r"\.bedgraph(\.gz)?$", ".bw", infile.name, flags=re.I)
    assert cnt == 1, infile
    outfile = infile.with_name(name)
    is_to_run = fs.is_outdated(outfile, infile)
    if infile.suffix == ".gz":
        infile = _gunzip(infile, if_=is_to_run)
        # gz or stdin are not accepted
    subp.run(["bedGraphToBigWig", infile, chrom_sizes, outfile], if_=is_to_run)
    return outfile


def faToTwoBit(fa_gz: Path) -> Path:  # noqa: N802
    outfile = fa_gz.with_suffix("").with_suffix(".2bit")
    subp.run(["faToTwoBit", fa_gz, outfile], if_=fs.is_outdated(outfile, fa_gz))
    return outfile


def faSize(genome_fa_gz: Path) -> Path:  # noqa: N802
    if not str(genome_fa_gz).endswith(("genome.fa.gz", os.devnull)):
        _log.warning(f"expecting *.genome.fa.gz: {genome_fa_gz}")
    outfile = genome_fa_gz.with_name("fasize.chrom.sizes")
    is_to_run = fs.is_outdated(outfile, genome_fa_gz)
    with subp.open_(outfile, "wb", if_=is_to_run) as fout:
        subp.run(["faSize", "-detailed", genome_fa_gz], stdout=fout, if_=is_to_run)
    _log.info(f"{outfile}")
    return outfile


def axt_chain(t2bit: Path, q2bit: Path, axtgz: Path) -> Path:
    chain = axtgz.with_suffix("").with_suffix(".chain")
    opts = subp.optargs(config["axtChain"], "-")
    args = ["axtChain", *opts, "stdin", t2bit, q2bit, chain]
    is_to_run = fs.is_outdated(chain, axtgz)
    with subp.popen_zcat(axtgz, if_=is_to_run) as zcat:
        subp.run(args, stdin=zcat.stdout, if_=is_to_run)
    return chain


def merge_sort_pre(chains: list[Path], target_sizes: Path, query_sizes: Path) -> Path:
    parent = {x.parent for x in chains}
    subdir = parent.pop()
    assert not parent, "chains are in the same directory"
    pre_chain = subdir / "pre.chain.gz"
    is_to_run = fs.is_outdated(pre_chain, chains)
    merge_cmd = ["chainMergeSort"] + [str(x) for x in chains]
    pre_cmd = f"chainPreNet stdin {target_sizes} {query_sizes} stdout"
    with (
        subp.popen(merge_cmd, if_=is_to_run, stdout=subp.PIPE) as merge,
        subp.popen(pre_cmd, if_=is_to_run, stdin=merge.stdout, stdout=subp.PIPE) as pre,
    ):
        return subp.gzip(pre.stdout, pre_chain, if_=is_to_run)


def chain_net_syntenic(pre_chain: Path, target_sizes: Path, query_sizes: Path) -> Path:
    syntenic_net = pre_chain.with_name("syntenic.net")
    is_to_run = fs.is_outdated(syntenic_net, pre_chain)
    cn_args = ["chainNet", *subp.optargs(config["chainNet"], "-")]
    cn_args += ["stdin", target_sizes, query_sizes, "stdout", "/dev/null"]
    ns_args = ["netSyntenic", "stdin", syntenic_net]
    with (
        subp.popen_zcat(pre_chain, if_=is_to_run) as zcat,
        subp.popen(cn_args, stdin=zcat.stdout, stdout=subp.PIPE, if_=is_to_run) as cn,
    ):
        subp.run(ns_args, stdin=cn.stdout, if_=is_to_run)
    return syntenic_net


def net_to_maf(net: Path, chain: Path, sing_maf: Path, target: str, query: str) -> Path:
    is_to_run = fs.is_outdated(sing_maf, [net, chain])
    sort = ["axtSort", "stdin", "stdout"]
    with (
        netToAxt(net, chain, target, query, if_=is_to_run) as ptoaxt,
        subp.popen(sort, stdin=ptoaxt.stdout, stdout=subp.PIPE, if_=is_to_run) as psort,
        axtToMaf(psort.stdout, target, query, sing_maf, if_=is_to_run) as paxttomaf,
    ):
        ptoaxt.stdout.close() if ptoaxt.stdout else None
        psort.stdout.close() if psort.stdout else None
        paxttomaf.communicate()
    return sing_maf


def netToAxt(  # noqa: N802
    net: Path, chain: Path, target: str, query: str, *, if_: bool = True
) -> subp.Popen[bytes]:
    target_2bit = faToTwoBit(api.genome_fa(target))
    query_2bit = faToTwoBit(api.genome_fa(query))
    opts = subp.optargs(config["netToAxt"], "-")
    args = ["netToAxt", *opts, net, chain, target_2bit, query_2bit, "stdout"]
    return subp.popen(args, stdout=subp.PIPE, if_=if_)


def axtToMaf(  # noqa: N802
    stdin: subp.FILE, target: str, query: str, sing_maf: Path, *, if_: bool = True
) -> subp.Popen[bytes]:
    opts = [f"-tPrefix={api.shorten(target)}.", f"-qPrefix={api.shorten(query)}."]
    args = ["axtToMaf", *opts, "stdin", api.fasize(target), api.fasize(query), sing_maf]
    return subp.popen(args, stdin=stdin, if_=if_)


def _popen(args: subp.Args, stdin: subp.FILE, **kwargs: str) -> subp.Popen[bytes]:
    opts = subp.optargs(kwargs, "-")
    return subp.popen([*args, *opts], stdin=stdin, stdout=subp.PIPE)


netFilter = functools.partial(_popen, ["netFilter", "stdin"])  # noqa: N816
chainFilter = functools.partial(_popen, ["chainFilter", "stdin"])  # noqa: N816


def chain_net_filter(file: Path, **kwargs: str) -> subp.Popen[bytes]:
    options = [f"-{k}={v}".removesuffix("=True") for k, v in kwargs.items()]
    if file.name.removesuffix(".gz").endswith(".net"):
        program = "netFilter"
    else:
        program = "chainFilter"
    cmd = [program, *options, file]
    return subp.popen(cmd, stdout=subp.PIPE)


def _gunzip(infile: Path, *, if_: bool = True) -> Path:
    fs.expect_suffix(infile, ".gz")
    outfile = infile.with_suffix("")
    if_ = if_ and fs.is_outdated(outfile, infile)
    subp.run(["unzstd", "-fk", infile, "-o", outfile], if_=if_)
    if if_ and not cli.dry_run:
        outfile.touch()
    return outfile


if __name__ == "__main__":
    main()
