"""Phylogenetic Analysis with Space/Time Models.

src: ./multiple/{target}/{clade}/{chromosome}/multiz.maf
dst: ./conservation/{target}/{clade}/phastcons.bw

http://compgen.cshl.edu/phast/
"""

import logging
import os
import re
import sys
from pathlib import Path

from aligons.db import api, phylo
from aligons.util import cli, config, fs, gff, subp

from . import htslib, kent

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("clade", type=Path)  # multiple/{target}/{clade}
    args = parser.parse_args(argv or None)
    run(args.clade)


def run(path_clade: Path) -> tuple[Path, Path]:
    outdir = prepare_conservation(path_clade)
    (_cons_mod, noncons_mod) = estimate_models(outdir)
    target = outdir.parent.name
    chrom_sizes = api.fasize(target)
    pattern = "chromosome.*/multiz.maf"
    mafs = [x for x in outdir.glob(pattern) if not x.parent.name.endswith(".Pt")]
    futures = [
        cli.thread_submit(phastCons, maf, chrom_sizes, None, noncons_mod)
        for maf in fs.sorted_naturally(mafs)
    ]
    chrom_bws: list[Path] = []
    chrom_beds: list[Path] = []
    for ft in futures:
        res = ft.result()
        chrom_bws.append(res[0])
        chrom_beds.append(res[1])
    bigwig = kent.bigWigCat(outdir / "phastcons.bw", chrom_bws)
    mostcons = bigwig.with_name("most-cons.bed.gz")
    concat_clean_mostcons(chrom_beds, mostcons)
    return bigwig, mostcons


def phastCons(  # noqa: N802
    msa: Path, chrom_sizes: Path, cons_mod: Path | None, noncons_mod: Path
) -> tuple[Path, Path]:
    conf = config.get("phastCons", {})
    bed = msa.with_name("most-cons.bed")
    wig = msa.with_name("post-probs.wig")
    bw = wig.with_suffix(".bw")
    estimate_base = msa.with_name("estimate")
    seqname = msa.parent.name.split(".", 1)[1]  # remove "chromosome."
    opts = ["--most-conserved", bed, "--score"]
    opts.extend(["--seqname", seqname])  # column 1 in bed (default: msa.stem)
    opts.extend(["--idpref", seqname])  # column 4 prefix in bed (default: msa.stem)
    if cons_mod is None:  # option 2
        opts.extend(["--estimate-rho", estimate_base])
        mod = f"{noncons_mod}"
    else:  # option 3
        mod = f"{cons_mod},{noncons_mod}"
    args = ["phastCons", *subp.optargs(conf), *opts, msa, mod]
    deps = [msa, noncons_mod]
    if_ = fs.is_outdated(bw, deps) or fs.is_outdated(bed, deps)
    with subp.open_(wig, "wb", if_=if_) as fout:
        subp.run(args, stdout=fout, if_=if_)
    fs.print_if_exists(bed)
    kent.wigToBigWig(wig, chrom_sizes)
    fs.print_if_exists(bw)
    e_mods = [estimate_base.with_suffix(x) for x in (".cons.mod", ".noncons.mod")]
    if e_mods[0].exists() and e_mods[1].exists():
        consEntropy(conf["target-coverage"], conf["expected-length"], *e_mods)
    return (bw, bed)


def prepare_conservation(multiple_target_clade: Path) -> Path:
    clade = multiple_target_clade.name
    target = multiple_target_clade.parent.name
    cons_dir = multiple_target_clade.parent.parent.with_name("conservation")
    outdir = cons_dir / target / clade
    if not cli.dry_run:
        outdir.mkdir(parents=True, exist_ok=True)
    for maf in multiple_target_clade.glob("chromosome*/multiz.maf"):
        ln = outdir / maf.parent.name / maf.name
        fs.symlink(maf, ln, relative=True)
    return outdir


def estimate_models(cons_target_clade: Path) -> tuple[Path, Path]:
    target = cons_target_clade.parent.name
    if not cli.dry_run:
        assert cds_gff3(target).exists()
    cons_mod = cons_target_clade / "phylofit.cons.mod"
    noncons_mod = cons_target_clade / "phylofit.noncons.mod"
    tree = phylo.get_subtree(cons_target_clade.name.split("-"), phylo.shorten_names)
    c_futures: list[cli.Future[Path]] = []
    n_futures: list[cli.Future[Path]] = []
    pool = cli.ThreadPool()
    for maf in cons_target_clade.glob("chromosome*/multiz.maf"):
        seqid = maf.parent.name.removeprefix("chromosome.")
        c_futures.append(pool.submit(make_cons_mod, maf, target, seqid, tree))
        n_futures.append(pool.submit(make_noncons_mod, maf, target, seqid, tree))
    phyloBoot([f.result() for f in c_futures], cons_mod)
    phyloBoot([f.result() for f in n_futures], noncons_mod)
    return (cons_mod, noncons_mod)


def make_cons_mod(maf: Path, species: str, seqid: str, tree: str) -> Path:
    codons_ss = msa_view_features(maf, species, seqid, conserved=True)
    codons_mods = phyloFit(codons_ss, tree, conserved=True)
    return most_conserved_mod(codons_mods)


def make_noncons_mod(maf: Path, species: str, seqid: str, tree: str) -> Path:
    _4d_codons_ss = msa_view_features(maf, species, seqid, conserved=False)
    _4d_sites_ss = msa_view_ss(_4d_codons_ss)
    return phyloFit(_4d_sites_ss, tree, conserved=False)[0]


def msa_view_features(maf: Path, species: str, seqid: str, *, conserved: bool) -> Path:
    """Calculate sufficient statistics with a custom GFF.

    - Seqid must be labeled with species names as in maf, e.g., osat.1, slyc.2.
    - `--4d` requires all features are on the same chromosome.
    - BED is not recognized as opposed to `--help` description.
    """
    cmd = ["msa_view", maf, "--in-format", "MAF", "--features", "/dev/stdin"]
    if conserved:
        outfile = maf.with_name("codons.ss")
        cmd.extend(["--catmap", "NCATS = 3; CDS 1-3"])
        cmd.extend(["--out-format", "SS", "--unordered-ss"])
        cmd.extend(["--reverse-groups", "transcript_id"])
    else:
        outfile = maf.with_name("4d-codons.ss")
        cmd.append("--4d")
    if_ = fs.is_outdated(outfile, maf)
    shortname = phylo.shorten(species)
    grep_cmd = ["zstdgrep", f"^{seqid}\t", cds_gff3(species)]
    with (
        subp.popen(grep_cmd, stdout=subp.PIPE, if_=if_) as grep,
        subp.popen_sd(r"^(\w)", f"{shortname}.$1", stdin=grep.stdout, if_=if_) as sd,
        subp.open_(outfile, "wb", if_=if_) as fout,
    ):
        subp.run(cmd, stdin=sd.stdout, stdout=fout, if_=if_)
    return outfile


def msa_view_ss(codons_ss: Path) -> Path:
    outfile = codons_ss.with_name("4d-sites.ss")
    s = f"msa_view {codons_ss!s} --in-format SS --out-format SS --tuple-size 1"
    is_to_run = fs.is_outdated(outfile, codons_ss)
    with subp.open_(outfile, "wb", if_=is_to_run) as fout:
        subp.run(s, if_=is_to_run, stdout=fout)
    return outfile


def phyloFit(ss: Path, tree: str, *, conserved: bool) -> list[Path]:  # noqa: N802
    if conserved:
        out_root = str(ss.with_name("codons"))
        outfiles = [Path(f"{out_root}.{i}.mod") for i in range(1, 4)]
        option = "--do-cats 1,2,3"
    else:
        out_root = str(ss.with_name("4d-sites"))
        outfiles = [Path(f"{out_root}.mod")]
        option = ""
    cmd = (
        f"phyloFit --tree {tree} --msa-format SS {option} --out-root {out_root} {ss!s}"
    )
    subp.run(cmd, if_=fs.is_outdated(outfiles[0], ss))
    return outfiles


def phyloBoot(mods: list[Path], outfile: Path) -> None:  # noqa: N802
    read_mods = ",".join(str(x) for x in mods)
    subp.run(
        f"phyloBoot --read-mods {read_mods} --output-average {outfile}",
        if_=fs.is_outdated(outfile, mods),
    )


def consEntropy(  # noqa: N802
    target_coverage: float, expected_length: int, cons_mod: Path, noncons_mod: Path
) -> Path:
    with cons_mod.open("rt") as fin:
        tree = extract_tree(fin.read())
    tree_size = len(phylo.extract_tip_names(tree))
    if tree_size > 12:  # noqa: PLR2004
        _log.debug(f"{tree_size=} is too large for consEntropy: {tree}")
        return Path(os.devnull)
    cov = str(target_coverage)
    exp_len = str(expected_length)
    args = ["consEntropy", cov, exp_len, cons_mod, noncons_mod]
    outfile = cons_mod.with_name("entropy.txt")
    if_ = fs.is_outdated(outfile, [cons_mod, noncons_mod])
    with subp.open_(outfile, "wb", if_=if_) as fout:
        subp.run(args, stdout=fout, if_=if_)
    return outfile


def most_conserved_mod(mods: list[Path]) -> Path:
    outfile = mods[0].with_name("phylofit.cons.mod")
    if not fs.is_outdated(outfile, mods) or cli.dry_run:
        return outfile
    shortest_length = sys.float_info.max
    conserved = ""
    for mod in mods:
        with mod.open() as fin:
            content = fin.read()
        lengths = phylo.extract_lengths(extract_tree(content))
        total = sum(lengths)
        _log.debug(f"{mod}: {total}")
        if total < shortest_length:
            shortest_length = total
            conserved = content
    _log.debug(f"{shortest_length=}")
    with outfile.open("w") as fout:
        fout.write(conserved)
    return outfile


def extract_tree(mod: str) -> str:
    m = re.search(r"TREE: (.+;)", mod)
    assert m, mod
    return m.group(1)


def cds_gff3(species: str) -> Path:
    genome_gff3 = api.genome_gff3(species)
    name = genome_gff3.name.removesuffix(".genome.gff3.gz") + ".cds.gff3.gz"
    cds = genome_gff3.with_name(name)
    if fs.is_outdated(cds, genome_gff3) and not cli.dry_run:
        x = gff.GFF(genome_gff3)
        x.body = x.body.filter(type="CDS")
        with htslib.popen_bgzip(cds) as bgzip:
            assert bgzip.stdin, cds
            x.write(bgzip.stdin)
    return fs.print_if_exists(cds)


def concat_clean_mostcons(beds: list[Path], outfile: Path) -> Path:
    if_ = fs.is_outdated(outfile, beds)
    with htslib.popen_bgzip(outfile, if_=if_) as bgzip:
        for infile in beds:
            with subp.popen_zcat(infile, if_=if_) as zcat:
                subp.run_sd(r"\+$", ".", stdin=zcat.stdout, stdout=bgzip.stdin, if_=if_)
    htslib.tabix(outfile)
    return fs.print_if_exists(outfile)


if __name__ == "__main__":
    main()
