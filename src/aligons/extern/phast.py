"""Phylogenetic Analysis with Space/Time Models.

<http://compgen.cshl.edu/phast/>
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
    parser.add_argument("indir", type=Path, help="multiple/{target}/{clade}")
    args = parser.parse_args(argv or None)
    run(args.indir)


def run(input_dir: Path) -> tuple[Path, Path]:
    """Perform phastCons for each chromosome and concatenate the results.

    :param input_dir: Path to the directory with chromosome MAFs:
        `./multiple/{target}/{clade}`
    :returns: Paths to bigWig and BED files:
        `./conservation/{target}/{clade}/{phastcons.bw,most-cons.bed.gz}`
    """
    outdir = prepare_conservation(input_dir)
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
    """Calculate conservation scores and predict conserved elements.

    :param msa: Input MAF file.
    :param chrom_sizes: Path to `fasize.chrom.sizes` of the target species.
    :param cons_mod: Path to `phylofit.cons.mod` or None.
    :param noncons_mod: Path to `4d-sites.mod`.
    :returns: Paths to bigWig and BED files: (`post-probs.bw`, `most-cons.bed`).
    """
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


def prepare_conservation(input_dir: Path) -> Path:
    """Prepare output directory and input MAFs.

    :param input_dir: Path to the directory with chromosome MAFs:
        `./multiple/{target}/{clade}`
    :returns: Path to the output directory:
        `./conservation/{target}/{clade}`
    """
    clade = input_dir.name
    target = input_dir.parent.name
    cons_dir = input_dir.parent.parent.with_name("conservation")
    outdir = cons_dir / target / clade
    if not cli.dry_run:
        outdir.mkdir(parents=True, exist_ok=True)
    for maf in input_dir.glob("chromosome*/multiz.maf"):
        ln = outdir / maf.parent.name / maf.name
        fs.symlink(maf, ln, relative=True)
    return outdir


def estimate_models(output_dir: Path) -> tuple[Path, Path]:
    """Estimate models for conserved and non-conserved sites.

    :param output_dir: Path to the output directory:
        `./conservation/{target}/{clade}`
    :returns: Paths to models: (`phylofit.cons.mod`, `phylofit.noncons.mod`)
    """
    target = output_dir.parent.name
    if not cli.dry_run:
        assert cds_gff3(target).exists()
    cons_mod = output_dir / "phylofit.cons.mod"
    noncons_mod = output_dir / "phylofit.noncons.mod"
    tree = phylo.get_subtree(output_dir.name.split("-"), phylo.shorten_names)
    c_futures: list[cli.Future[Path]] = []
    n_futures: list[cli.Future[Path]] = []
    pool = cli.ThreadPool()
    for maf in output_dir.glob("chromosome*/multiz.maf"):
        seqid = maf.parent.name.removeprefix("chromosome.")
        c_futures.append(pool.submit(make_cons_mod, maf, target, seqid, tree))
        n_futures.append(pool.submit(make_noncons_mod, maf, target, seqid, tree))
    phyloBoot([f.result() for f in c_futures], cons_mod)
    phyloBoot([f.result() for f in n_futures], noncons_mod)
    return (cons_mod, noncons_mod)


def make_cons_mod(maf: Path, species: str, seqid: str, tree: str) -> Path:
    """Estimate phylogenetic models for conserved sites.

    :param maf: Input MAF file.
    :param species: Species name to get its genome GFF3.
    :param seqid: Chromosome name.
    :param tree: Newick tree string.
    :returns: Path to the generated model file: `phylofit.cons.mod`.
    """
    codons_ss = msa_view_features(maf, species, seqid, conserved=True)
    codons_mods = phyloFit(codons_ss, tree, conserved=True)
    return most_conserved_mod(codons_mods)


def make_noncons_mod(maf: Path, species: str, seqid: str, tree: str) -> Path:
    """Estimate phylogenetic models for non-conserved 4d sites.

    See <http://compgen.cshl.edu/phast/phyloFit-tutorial.php>

    :param maf: Input MAF file.
    :param species: Species name to get its genome GFF3.
    :param seqid: Chromosome name.
    :param tree: Newick tree string.
    :returns: Path to the generated model file: `4d-sites.mod`.
    """
    _4d_codons_ss = msa_view_features(maf, species, seqid, conserved=False)
    _4d_sites_ss = msa_view_ss(_4d_codons_ss)
    return phyloFit(_4d_sites_ss, tree, conserved=False)[0]


def msa_view_features(maf: Path, species: str, seqid: str, *, conserved: bool) -> Path:
    """Extract sufficient statistics (SS) for codon sites separately.

    - Seqid must be labeled with species names as in maf, e.g., osat.1, slyc.2.
    - `--4d` requires all features are on the same chromosome.
    - BED is not recognized as opposed to `--help` description.

    See <http://compgen.cshl.edu/phast/help-pages/msa_view.txt>

    :param maf: Input MAF file.
    :param species: Species name to get its genome GFF3.
    :param seqid: Chromosome name.
    :param conserved: Whether to extract SS for conserved codon sites.
    :returns: Path to the generated SS file: `codons.ss` or `4d-codons.ss`.
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
    """Extract 4d site in SS.

    :param codons_ss: Input codon SS file from `msa_view`: `4d-codons.ss`.
    :returns: Path to the generated SS file: `4d-sites.ss`.
    """
    outfile = codons_ss.with_name("4d-sites.ss")
    s = f"msa_view {codons_ss!s} --in-format SS --out-format SS --tuple-size 1"
    is_to_run = fs.is_outdated(outfile, codons_ss)
    with subp.open_(outfile, "wb", if_=is_to_run) as fout:
        subp.run(s, if_=is_to_run, stdout=fout)
    return outfile


def phyloFit(ss: Path, tree: str, *, conserved: bool) -> list[Path]:  # noqa: N802
    """Estimate phylogenetic models from SS.

    See <http://compgen.cshl.edu/phast/phyloFit-tutorial.php>

    :param ss: Input sufficient statistics file from `msa_view`.
    :param tree: Newick tree string.
    :param conserved: Whether to estimate models for conserved codon sites.
    :returns: List of generated model files: `codons.{1,2,3}.mod` or `4d-sites.mod`.
    """
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
    """Calculate the average of all input models.

    See <http://compgen.cshl.edu/phast/help-pages/phyloBoot.txt>

    :param mods: List of model files to average.
    :param outfile: Output model file.
    """
    read_mods = ",".join(str(x) for x in mods)
    subp.run(
        f"phyloBoot --read-mods {read_mods} --output-average {outfile}",
        if_=fs.is_outdated(outfile, mods),
    )


def consEntropy(  # noqa: N802
    target_coverage: float, expected_length: int, cons_mod: Path, noncons_mod: Path
) -> Path:
    """Compute the relative entropy and some statistics of the phylogenetic models.

    See <http://compgen.cshl.edu/phast/help-pages/consEntropy.txt>

    :param target_coverage:
    :param expected_length:
    :param cons_mod: Path to `entropy.txt`.
    """
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
    """Select and copy the most conserved model among codon models.

    :param mods: List of codon model files: `codons.{1,2,3}.mod`.
    :returns: Path to the copy of the selected model: `phylofit.cons.mod`.
    """
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
    """Extract Newick tree string from model file content.

    :param mod: Content of a model file.
    :returns: Newick tree string.
    """
    m = re.search(r"TREE: (.+;)", mod)
    assert m, mod
    return m.group(1)


def cds_gff3(species: str) -> Path:
    """Extract CDS features from genome GFF3.

    :param species: Species name to get genome GFF3.
    :returns: Path to the generated GFF3 file with only CDS features.
    """
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
    """Concatenate BED files and clean up the score column.

    :param beds: List of BED files.
    :param outfile: Concatenated output BED file.
    :returns: The same path as `outfile`.
    """
    if_ = fs.is_outdated(outfile, beds)
    with htslib.popen_bgzip(outfile, if_=if_) as bgzip:
        for infile in beds:
            with subp.popen_zcat(infile, if_=if_) as zcat:
                subp.run_sd(r"\+$", ".", stdin=zcat.stdout, stdout=bgzip.stdin, if_=if_)
    htslib.tabix(outfile)
    return fs.print_if_exists(outfile)


if __name__ == "__main__":
    main()
