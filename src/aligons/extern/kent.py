"""Integrate chromosome wigs into a genome-wide bigwig.

src: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz
dst: ./multiple/{target}/{clade}/phastcons.bw

https://github.com/ucscGenomeBrowser/kent
"""
import gzip
import logging
import os
import re
import shutil
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
    p = subp.popen(args, if_=is_to_run, stdin=subp.PIPE)
    if is_to_run:
        assert p.stdin is not None
        for wig in wigs:
            with gzip.open(wig, "rb") as fin:
                p.stdin.write(fin.read())
                p.stdin.flush()
    p.communicate()
    return outfile


def bigWigInfo(path: Path) -> bytes:  # noqa: N802
    args = ["bigWigInfo", path]
    return subp.run(args, stdout=subp.PIPE, text=True).stdout


def bedGraphToBigWig(infile: Path, chrom_sizes: Path) -> Path:  # noqa: N802
    name = re.sub(r"\.bedgraph(\.gz)?$", ".bw", infile.name, flags=re.I)
    outfile = infile.with_name(name)
    if fs.is_outdated(outfile, infile):
        bedgraph = gunzip(infile)
        # gz or stdin are not accepted
        subp.run(["bedGraphToBigWig", bedgraph, chrom_sizes, outfile])
    return outfile


def faToTwoBit(fa_gz: Path) -> Path:  # noqa: N802
    outfile = fa_gz.with_suffix("").with_suffix(".2bit")
    subp.run(["faToTwoBit", fa_gz, outfile], if_=fs.is_outdated(outfile, fa_gz))
    return outfile


def faSize(genome_fa_gz: Path) -> Path:  # noqa: N802
    if not str(genome_fa_gz).endswith(("genome.fa.gz", os.devnull)):
        _log.warning(f"expecting *.genome.fa.gz: {genome_fa_gz}")
    outfile = genome_fa_gz.with_name("fasize.chrom.sizes")
    if fs.is_outdated(outfile, genome_fa_gz) and not cli.dry_run:
        with outfile.open("wb") as fout:
            subp.run(["faSize", "-detailed", genome_fa_gz], stdout=fout)
    _log.info(f"{outfile}")
    return outfile


def axt_chain(t2bit: Path, q2bit: Path, axtgz: Path) -> Path:
    chain = axtgz.with_suffix("").with_suffix(".chain")
    cmd = "axtChain"
    cmd += subp.optjoin(config["axtChain"], "-")
    cmd += f" stdin {t2bit} {q2bit} {chain}"
    is_to_run = fs.is_outdated(chain, axtgz)
    p = subp.popen(cmd, if_=is_to_run, stdin=subp.PIPE)
    if is_to_run and not cli.dry_run:
        assert p.stdin is not None
        with gzip.open(axtgz, "rb") as fin:
            shutil.copyfileobj(fin, p.stdin)
            p.stdin.close()
    p.communicate()
    return chain


def merge_sort_pre(chains: list[Path], target_sizes: Path, query_sizes: Path) -> Path:
    parent = {x.parent for x in chains}
    subdir = parent.pop()
    assert not parent, "chains are in the same directory"
    pre_chain = subdir / "pre.chain.gz"
    is_to_run = fs.is_outdated(pre_chain, chains)
    merge_cmd = ["chainMergeSort"] + [str(x) for x in chains]
    merge = subp.popen(merge_cmd, if_=is_to_run, stdout=subp.PIPE)
    assert merge.stdout is not None
    pre_cmd = f"chainPreNet stdin {target_sizes} {query_sizes} stdout"
    pre = subp.popen(pre_cmd, if_=is_to_run, stdin=merge.stdout, stdout=subp.PIPE)
    merge.stdout.close()
    if is_to_run and not cli.dry_run:
        (stdout, _stderr) = pre.communicate()
        with gzip.open(pre_chain, "wb") as fout:
            fout.write(stdout)
    return pre_chain


def chain_net_syntenic(pre_chain: Path, target_sizes: Path, query_sizes: Path) -> Path:
    syntenic_net = pre_chain.with_name("syntenic.net")
    is_to_run = fs.is_outdated(syntenic_net, pre_chain)
    cn_cmd = "chainNet"
    cn_cmd += subp.optjoin(config["chainNet"], "-")
    cn_cmd += f" stdin {target_sizes} {query_sizes} stdout /dev/null"
    ns_cmd = f"netSyntenic stdin {syntenic_net}"
    cn = subp.popen(cn_cmd, if_=is_to_run, stdin=subp.PIPE, stdout=subp.PIPE)
    ns = subp.popen(ns_cmd, if_=is_to_run, stdin=subp.PIPE)
    content = b""
    if is_to_run and not cli.dry_run:
        with gzip.open(pre_chain, "rb") as fout:
            content = fout.read()
    (cn_out, _) = cn.communicate(content)
    ns.communicate(cn_out)
    return syntenic_net


def net_to_maf(
    net: Path, chain_gz: Path, sing_maf: Path, target: str, query: str
) -> Path:
    target_2bit = api.genome_2bit(target)
    query_2bit = api.genome_2bit(query)
    target_sizes = api.fasize(target)
    query_sizes = api.fasize(query)
    tprefix = api.shorten(target)
    qprefix = api.shorten(query)
    is_to_run = fs.is_outdated(sing_maf, [net, chain_gz])
    toaxt_cmd = "netToAxt"
    toaxt_cmd += subp.optjoin(config["netToAxt"], "-")
    toaxt_cmd += f" {net} stdin {target_2bit} {query_2bit} stdout"
    toaxt = subp.popen(toaxt_cmd, if_=is_to_run, stdin=subp.PIPE, stdout=subp.PIPE)
    assert toaxt.stdin is not None
    assert toaxt.stdout is not None
    if is_to_run and not cli.dry_run:
        with gzip.open(chain_gz, "rb") as fout:
            shutil.copyfileobj(fout, toaxt.stdin)
            toaxt.stdin.close()
    sort = subp.popen(
        "axtSort stdin stdout", if_=is_to_run, stdin=toaxt.stdout, stdout=subp.PIPE
    )
    toaxt.stdout.close()
    assert sort.stdout is not None
    axttomaf_cmd = (
        f"axtToMaf -tPrefix={tprefix}. -qPrefix={qprefix}. stdin"
        f" {target_sizes} {query_sizes} {sing_maf}"
    )
    axttomaf = subp.popen(axttomaf_cmd, if_=is_to_run, stdin=sort.stdout)
    sort.stdout.close()
    axttomaf.communicate()
    return sing_maf


def chain_net_filter(file: Path, **kwargs: str) -> bytes:
    options = [f"-{k}={v}".removesuffix("=True") for k, v in kwargs.items()]
    if file.name.removesuffix(".gz").endswith(".net"):
        program = "netFilter"
    else:
        program = "chainFilter"
    cmd = [program, *options, file]
    p = subp.run(cmd, stdout=subp.PIPE)
    return p.stdout


def gunzip(infile: Path) -> Path:
    if infile.suffix != ".gz":
        _log.debug(f"non-gz {infile}")
        return infile
    outfile = infile.with_suffix("")
    if fs.is_outdated(outfile, infile):
        subp.run(["gunzip", "-fk", infile])
        subp.run(["touch", outfile])
    return outfile


if __name__ == "__main__":
    main()
