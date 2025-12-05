"""UCSC Genome Browser source.

<https://github.com/ucscGenomeBrowser/kent>
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
    """Integrate chromosome wigs into a genome-wide bigwig.

    :param out_bw: Output bigWig file.
    :param in_bws: Input bigWig files.
    :returns: The same path as `out_bw`.
    """
    if len(in_bws) == 1:
        args = ["cp", *in_bws, out_bw]
    else:
        args = ["bigWigCat", out_bw, *in_bws]
    subp.run(args, if_=fs.is_outdated(out_bw, in_bws))
    return out_bw


def bigWigInfo(path: Path) -> bytes:  # noqa: N802
    """Get bigWig file information.

    :param path: Input bigWig file.
    :returns: The stdout bytes of `bigWigInfo` command.
    """
    args = ["bigWigInfo", path]
    return subp.run(args, stdout=subp.PIPE, text=True).stdout


def wigToBigWig(wig: Path, chrom_sizes: Path, *, keep: bool = False) -> Path:  # noqa: N802
    """Convert wig to bigWig format.

    :param wig: Input wig file.
    :param chrom_sizes: `fasize.chrom.sizes`.
    :param keep: Whether to keep the intermediate wig file.
    :returns: Output bigWig file.
    """
    bw = wig.with_suffix(".bw")
    args: subp.Args = ["wigToBigWig"]
    args.extend(["-fixedSummaries", "-keepAllChromosomes"])  # for bigWigCat
    args.extend([wig, chrom_sizes, bw])  # mustBeReadableAndRegularFile
    subp.run(args, if_=fs.is_outdated(bw, wig))
    if wig.exists() and bw.exists() and not keep and not cli.dry_run:
        wig.unlink()
    return bw


def bedGraphToBigWig(infile: Path, chrom_sizes: Path) -> Path:  # noqa: N802
    """Convert bedGraph to bigWig format.

    :param infile: Input bedGraph file. Decompressed if compressed.
    :param chrom_sizes: `fasize.chrom.sizes`.
    :returns: Output bigWig file.
    """
    name, cnt = re.subn(r"\.bedgraph(\.gz)?$", ".bw", infile.name, flags=re.IGNORECASE)
    assert cnt == 1, infile
    outfile = infile.with_name(name)
    is_to_run = fs.is_outdated(outfile, infile)
    if infile.suffix == ".gz":
        infile = _gunzip(infile, if_=is_to_run)  # mustBeReadableAndRegularFile
    subp.run(["bedGraphToBigWig", infile, chrom_sizes, outfile], if_=is_to_run)
    return outfile


def bigWigToBed(infile: Path, min_width: int = 15) -> bytes:  # noqa: N802
    """Extract ranges from bigWig file and output in BED format.

    :param infile: Input bigWig file.
    :param min_width: Minimum width of regions to output.
    :returns: BED format data as bytes.
    """
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
    """Convert FASTA to 2bit format.

    :param fa: FASTA file or Popen stdout with FASTA data.
    :param twobit: Output 2bit file. Optional if `fa` is a Path.
    :returns: Output 2bit file.
    """
    if isinstance(fa, Path | cli.Future):
        return _faToTwoBit_f(fa, twobit)
    assert twobit is not None, "`twobit` name must be specified if input is a pipe"
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
    """Read `fasize.chrom.sizes` of the given genome.

    :param genome: A bgzipped genome FASTA file.
    :returns: A dictionary of chromosome sizes.
    """
    fasize = faSize(genome)
    res: dict[str, int] = {}
    with fasize.open("rt") as fin:
        for line in fin:
            (seqid, size) = line.split(maxsplit=1)
            res[seqid] = int(size)
    return res


def faSize(genome_fa_gz: Path | cli.Future[Path]) -> Path:  # noqa: N802
    """Summarize chromosome sizes from a genome FASTA.

    :param genome_fa_gz: A bgzipped genome FASTA file.
    :returns: Output file named `fasize.chrom.sizes`.
    """
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
    """Chain together axt alignments.

    :param axt: Gzipped AXT file.
    :param t2bit: Target 2bit file.
    :param q2bit: Query 2bit file.
    :returns: Gzipped chain file: `{axt.stem}.chain.gz`
    """
    chain = axt.with_suffix("").with_suffix(".chain.gz")
    opts = subp.optargs(config["axtChain"], "-")
    args = ["axtChain", *opts, axt, t2bit, q2bit, "stdout"]
    if_ = fs.is_outdated(chain, axt)
    with subp.popen(args, stdout=subp.PIPE, if_=if_) as p:
        return subp.gzip(p.stdout, chain, if_=if_)


def chainMergeSort(chains: list[Path]) -> Path:  # noqa: N802
    """Combine sorted files into larger sorted file.

    :param chains: Gzipped chain files to merge.
    :returns: Merged chain file named `target.chain.gz`
    """
    parent = {x.parent for x in chains}
    subdir = parent.pop()
    assert not parent, "chains are in the same directory"
    outfile = subdir / "target.chain.gz"
    if_ = fs.is_outdated(outfile, chains)
    merge_cmd = ["chainMergeSort", *chains]
    with subp.popen(merge_cmd, if_=if_, stdout=subp.PIPE) as merge:
        return subp.gzip(merge.stdout, outfile, if_=if_)


def chain_net(chain: Path, target_sizes: Path, query_sizes: Path) -> tuple[Path, Path]:
    """Make alignment nets with synteny information from a chain.

    - `chainPreNet`: Remove chains that don't have a chance of being netted.
    - `chainNet`: Make alignment nets out of chains.
    - `netSyntenic`: Add synteny info to net.

    :param chain: Gzipped chain file.
    :param target_sizes: Target `fasize.chrom.sizes`.
    :param query_sizes: Query `fasize.chrom.sizes`.
    :returns: Gzipped syntenic net files:
        (`target.net.gz`, `query.net.gz`)
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
    """Add synteny info to net.

    :param net: Input net file.
    :param syn_net: Output syntenic net file.
    :returns: Same path as `syn_net`.
    """
    args = ["netSyntenic", net, "stdout"]
    with subp.popen(args, stdout=subp.PIPE, if_=if_) as p:
        subp.gzip(p.stdout, syn_net, if_=if_)
    if net.exists() and if_ and not cli.dry_run:
        net.unlink()
    return syn_net


def net_to_maf(net: Path, chain: Path, sing_maf: Path, target: str, query: str) -> Path:
    """Convert net and chain to MAF format with `netToAxt | axtSort | axtToMaf`.

    :param net: Input net file.
    :param chain: Input chain file.
    :param sing_maf: Output MAF file.
    :param target: Target species name.
    :param query: Query species name.
    :returns: Same path as `sing_maf`.
    """
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
    """Convert net and chain to AXT format.

    :param net: Input net file.
    :param chain: Input chain file.
    :param target: Target species name.
    :param query: Query species name.
    :returns: Popen object with stdout of AXT data.
    """
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
    """Convert AXT to MAF format.

    :param stdin: Popen stdout with AXT data.
    :param target: Target species name.
    :param query: Query species name.
    :param sing_maf: Output MAF file.
    :returns: Popen object without stdout.
    """
    opts = [f"-tPrefix={api.shorten(target)}.", f"-qPrefix={api.shorten(query)}."]
    args = ["axtToMaf", *opts, "stdin", api.fasize(target), api.fasize(query), sing_maf]
    return subp.popen(args, stdin=stdin, if_=if_)


def net_chain_subset(net: Path, chain: Path) -> Path:
    """Filter and join chains that appear in the net.

    - `netChainSubset`: Create chain file with subset of chains that appear in the net.
    - `chainStitchId`: Join chain fragments with the same chain ID into a single chain.

    :param net: Input net file.
    :param chain: Input chain file.
    :returns: Gzipped chain file:
        `{net.stem}.over.chain.gz`
    """
    outfile = net.with_suffix("").with_suffix(".over.chain.gz")
    ncs = ["netChainSubset", net, chain, "stdout"]
    csi = ["chainStitchId", "stdin", "stdout"]
    if_ = fs.is_outdated(outfile, [net, chain])
    with (
        subp.popen(ncs, stdout=subp.PIPE, if_=if_) as p_ncs,
        subp.popen(csi, stdin=p_ncs.stdout, stdout=subp.PIPE, if_=if_) as p_csi,
    ):
        return subp.gzip(p_csi.stdout, outfile, if_=if_)


def _popen(args: subp.Args, stdin: subp.FILE, **kwargs: str) -> subp.Popen[bytes]:
    opts = subp.optargs(kwargs, "-")
    return subp.popen([*args, *opts], stdin=stdin, stdout=subp.PIPE)


netFilter = functools.partial(_popen, ["netFilter", "stdin"])  # noqa: N816
chainFilter = functools.partial(_popen, ["chainFilter", "stdin"])  # noqa: N816


def _gunzip(infile: Path, *, if_: bool = True) -> Path:
    """Decompress a gzipped file while keeping the original.

    :param infile: Input gzipped file.
    :returns: Decompressed file.
    """
    fs.expect_suffix(infile, ".gz")
    outfile = infile.with_suffix("")
    if_ = if_ and fs.is_outdated(outfile, infile)
    subp.run(["unzstd", "-fk", infile, "-o", outfile], if_=if_)
    if if_ and not cli.dry_run:
        outfile.touch()
    return outfile


if __name__ == "__main__":
    main()
