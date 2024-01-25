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


def axtChain(axt: Path, t2bit: Path, q2bit: Path) -> Path:  # noqa: N802
    """Chain together axt alignments."""
    chain = axt.with_suffix("").with_suffix(".chain.gz")
    opts = subp.optargs(config["axtChain"], "-")
    args = ["axtChain", *opts, axt, t2bit, q2bit, "stdout"]
    if_ = fs.is_outdated(chain, axt)
    with subp.popen(args, stdout=subp.PIPE, if_=if_) as p:
        return subp.gzip(p.stdout, chain, if_=if_)


def chainMergeSort(chains: list[Path]) -> Path:  # noqa: N802
    """Combine sorted files into larger sorted file."""
    parent = {x.parent for x in chains}
    subdir = parent.pop()
    assert not parent, "chains are in the same directory"
    outfile = subdir / "target.chain.gz"
    if_ = fs.is_outdated(outfile, chains)
    merge_cmd = ["chainMergeSort", *chains]
    with subp.popen(merge_cmd, if_=if_, stdout=subp.PIPE) as merge:
        return subp.gzip(merge.stdout, outfile, if_=if_)


def chain_net(chain: Path, target_sizes: Path, query_sizes: Path) -> tuple[Path, Path]:
    """Make alignment nets out of chains.

    chainPreNet - Remove chains that don't have a chance of being netted.
    """
    tnet = chain.with_name("target.net")
    qnet = chain.with_name("query.net")
    tsyn = tnet.with_suffix(".net.gz")
    qsyn = qnet.with_suffix(".net.gz")
    if_ = fs.is_outdated(tsyn, chain) or fs.is_outdated(qsyn, chain)
    pre_args = ["chainPreNet", chain, target_sizes, query_sizes, "stdout"]
    opts = subp.optargs(config["chainNet"], "-")
    cn_args = ["chainNet", *opts, "stdin", target_sizes, query_sizes, tnet, qnet]
    with subp.popen(pre_args, stdout=subp.PIPE, if_=if_) as p:
        try:
            subp.run(cn_args, stdin=p.stdout, if_=if_)
        except subp.CalledProcessError as e:
            # printMem() in chainNet/netSyntenic works only on Linux.
            if not e.stderr.startswith(b"Couldn't open /proc/self/stat"):
                raise
    return (netSyntenic(tnet, tsyn, if_=if_), netSyntenic(qnet, qsyn, if_=if_))


def netSyntenic(net: Path, synnet: Path, *, if_: bool) -> Path:  # noqa: N802
    """Add synteny info to net."""
    args = ["netSyntenic", net, "stdout"]
    with subp.popen(args, stdout=subp.PIPE, if_=if_) as p:
        subp.gzip(p.stdout, synnet, if_=if_)
    if net.exists() and if_ and not cli.dry_run:
        net.unlink()
    return synnet


def net_to_maf(net: Path, chain: Path, sing_maf: Path, target: str, query: str) -> Path:
    if_ = fs.is_outdated(sing_maf, [net, chain])
    if net.name.startswith("query"):
        target, query = query, target
    sort = ["axtSort", "stdin", "stdout"]
    with (
        netToAxt(net, chain, target, query, if_=if_) as ptoaxt,
        subp.popen(sort, stdin=ptoaxt.stdout, stdout=subp.PIPE, if_=if_) as psort,
        axtToMaf(psort.stdout, target, query, sing_maf, if_=if_) as paxttomaf,
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
    if net.name.startswith("query"):
        opts.append("-qChain")
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
