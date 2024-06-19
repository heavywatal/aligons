"""Integrate chromosome wigs into a genome-wide bigwig.

src: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz
dst: ./multiple/{target}/{clade}/phastcons.bw

https://github.com/ucscGenomeBrowser/kent
"""

import functools
import io
import logging
import os
import re
from collections.abc import Sequence
from pathlib import Path

from aligons.db import api
from aligons.util import cli, config, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("clade", type=Path)
    args = parser.parse_args(argv or None)
    _run(args.clade)


def _run(clade: Path) -> Path:
    name = "phastcons.bw"
    wigs = [p / name for p in fs.sorted_naturally(clade.glob("chromosome.*"))]
    if (bigwig := bigWigCat(clade / name, wigs)).exists():
        _log.info(bigwig)
        _log.info(bigWigInfo(bigwig).decode().rstrip())
    return bigwig


def bigWigCat(out_bw: Path, in_bws: Sequence[Path]) -> Path:  # noqa: N802
    if len(in_bws) == 1:
        args = ["cp", *in_bws, out_bw]
    else:
        args = ["bigWigCat", out_bw, *in_bws]
    subp.run(args, if_=fs.is_outdated(out_bw, in_bws))
    return out_bw


def bigWigInfo(path: Path) -> bytes:  # noqa: N802
    args = ["bigWigInfo", path]
    return subp.run(args, stdout=subp.PIPE, text=True).stdout


def wigToBigWig(wig: Path, chrom_sizes: Path, *, keep: bool = False) -> Path:  # noqa: N802
    bw = wig.with_suffix(".bw")
    args: subp.Args = ["wigToBigWig"]
    args.extend(["-fixedSummaries", "-keepAllChromosomes"])  # for bigWigCat
    args.extend([wig, chrom_sizes, bw])  # mustBeReadableAndRegularFile
    subp.run(args, if_=fs.is_outdated(bw, wig))
    if wig.exists() and bw.exists() and not keep and not cli.dry_run:
        wig.unlink()
    return bw


def bedGraphToBigWig(infile: Path, chrom_sizes: Path) -> Path:  # noqa: N802
    name, cnt = re.subn(r"\.bedgraph(\.gz)?$", ".bw", infile.name, flags=re.I)
    assert cnt == 1, infile
    outfile = infile.with_name(name)
    is_to_run = fs.is_outdated(outfile, infile)
    if infile.suffix == ".gz":
        infile = _gunzip(infile, if_=is_to_run)  # mustBeReadableAndRegularFile
    subp.run(["bedGraphToBigWig", infile, chrom_sizes, outfile], if_=is_to_run)
    return outfile


def bigWigToBed(infile: Path, min_width: int = 15) -> bytes:  # noqa: N802
    patt = re.compile(rb"^fixedStep chrom=(\S+) start=(\d+)")
    chrom = b""
    start = 0
    scores: list[float] = []
    buf = io.BytesIO()

    def write() -> None:
        if (w := len(scores)) >= min_width:
            s = int(sum(scores) / (0.001 * w))
            buf.write(chrom + f"\t{start}\t{start + w}\t.\t{s}\t.\n".encode())

    with subp.popen(["bigWigToWig", infile, "stdout"], stdout=subp.PIPE) as wig:
        assert wig.stdout is not None
        for line in wig.stdout:
            if mobj := patt.match(line):
                write()
                chrom = mobj.group(1)
                start = int(mobj.group(2)) - 1
                scores = []
            else:
                scores.append(float(line))
        write()
    return buf.getvalue()


def faToTwoBit(  # noqa: N802
    fa: Path | cli.Future[Path] | subp.FILE,
    twobit: Path | None = None,
    *,
    if_: bool = True,
) -> Path:
    if isinstance(fa, Path | cli.Future):
        return _faToTwoBit_f(fa, twobit)
    assert twobit is not None
    return _faToTwoBit_s(fa, twobit, if_=if_)


def _faToTwoBit_s(stdin: subp.FILE, twobit: Path, *, if_: bool = True) -> Path:  # noqa: N802
    fs.expect_suffix(twobit, ".2bit")
    subp.run(["faToTwoBit", "stdin", twobit], stdin=stdin, if_=if_)
    return twobit


def _faToTwoBit_f(fa: Path | cli.Future[Path], twobit: Path | None = None) -> Path:  # noqa: N802
    fa = cli.result(fa)
    if twobit is None:
        twobit = (fa.parent / fa.name.removesuffix(".gz")).with_suffix(".2bit")
    if_ = fs.is_outdated(twobit, fa)
    with subp.popen_zcat(fa, if_=if_) as zcat:
        return _faToTwoBit_s(zcat.stdout, twobit, if_=if_)


def read_fasize(genome: Path) -> dict[str, int]:
    fasize = faSize(genome)
    res: dict[str, int] = {}
    with fasize.open("rt") as fin:
        for line in fin:
            (seqid, size) = line.split(maxsplit=1)
            res[seqid] = int(size)
    return res


def faSize(genome_fa_gz: Path | cli.Future[Path]) -> Path:  # noqa: N802
    genome_fa_gz = cli.result(genome_fa_gz)
    if not str(genome_fa_gz).endswith(("genome.fa.gz", os.devnull)):
        _log.warning(f"expecting *.genome.fa.gz: {genome_fa_gz}")
    outfile = genome_fa_gz.with_name("fasize.chrom.sizes")
    if_ = fs.is_outdated(outfile, genome_fa_gz)
    with (
        subp.popen_zcat(genome_fa_gz, if_=if_) as zcat,
        subp.open_(outfile, "wb", if_=if_) as fout,
    ):
        cmd = ["faSize", "-detailed", "stdin"]
        subp.run(cmd, stdin=zcat.stdout, stdout=fout, if_=if_)
    return fs.print_if_exists(outfile)


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
    t_net = chain.with_name("target.net")
    q_net = chain.with_name("query.net")
    t_syn = t_net.with_suffix(".net.gz")
    q_syn = q_net.with_suffix(".net.gz")
    if_ = fs.is_outdated(t_syn, chain) or fs.is_outdated(q_syn, chain)
    pre_args = ["chainPreNet", chain, target_sizes, query_sizes, "stdout"]
    opts = subp.optargs(config["chainNet"], "-")
    cn_args = ["chainNet", *opts, "stdin", target_sizes, query_sizes, t_net, q_net]
    with subp.popen(pre_args, stdout=subp.PIPE, if_=if_) as p:
        try:
            subp.run(cn_args, stdin=p.stdout, if_=if_)
        except subp.CalledProcessError as e:
            # printMem() in chainNet/netSyntenic works only on Linux.
            if not e.stderr.startswith(b"Couldn't open /proc/self/stat"):
                raise
    return (netSyntenic(t_net, t_syn, if_=if_), netSyntenic(q_net, q_syn, if_=if_))


def netSyntenic(net: Path, syn_net: Path, *, if_: bool) -> Path:  # noqa: N802
    """Add synteny info to net."""
    args = ["netSyntenic", net, "stdout"]
    with subp.popen(args, stdout=subp.PIPE, if_=if_) as p:
        subp.gzip(p.stdout, syn_net, if_=if_)
    if net.exists() and if_ and not cli.dry_run:
        net.unlink()
    return syn_net


def net_to_maf(net: Path, chain: Path, sing_maf: Path, target: str, query: str) -> Path:
    if_ = fs.is_outdated(sing_maf, [net, chain])
    if net.name.startswith("query"):
        target, query = query, target
    sort = ["axtSort", "stdin", "stdout"]
    with (
        netToAxt(net, chain, target, query, if_=if_) as p_toa,
        subp.popen(sort, stdin=p_toa.stdout, stdout=subp.PIPE, if_=if_) as p_sort,
        axtToMaf(p_sort.stdout, target, query, sing_maf, if_=if_) as p_axttomaf,
    ):
        p_toa.stdout.close() if p_toa.stdout else None
        p_sort.stdout.close() if p_sort.stdout else None
        p_axttomaf.communicate()
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
