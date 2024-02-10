"""Phylogenetic Analysis with Space/Time Models.

src: ./multiple/{target}/{clade}/{chromosome}/multiz.maf
dst: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz

http://compgen.cshl.edu/phast/
"""
import csv
import gzip
import itertools
import logging
import re
import sys
from pathlib import Path

from aligons.db import api, phylo
from aligons.util import cli, config, fs, subp

from . import kent

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("clade", type=Path)  # multiple/{target}/{clade}
    args = parser.parse_args(argv or None)
    if args.clean:
        clean(args.clade)
        return
    run(args.clade)


def run(path_clade: Path) -> tuple[list[Path], list[Path]]:
    (_cons_mod, noncons_mod) = prepare_mods(path_clade)
    target = path_clade.parent.name
    chrom_sizes = api.fasize(target)
    futures = [
        cli.thread_submit(phastCons, maf, chrom_sizes, None, noncons_mod)
        for maf in fs.sorted_naturally(path_clade.glob("chromosome.*/multiz.maf"))
    ]
    bigwigs: list[Path] = []
    beds: list[Path] = []
    for ft in futures:
        res = ft.result()
        bigwigs.append(res[0])
        beds.append(res[1])
    return (bigwigs, beds)


def phastCons(  # noqa: N802
    msa: Path, chrom_sizes: Path, cons_mod: Path | None, noncons_mod: Path
) -> tuple[Path, Path]:
    conf = config.get("phastCons", {})
    bed = msa.with_name("most-cons.bed")
    wig = msa.with_name("post-probs.wig")
    bw = wig.with_suffix(".bw")
    estimate_base = msa.with_name("estimate")
    seqname = msa.parent.name.split(".", 1)[1]  # remove "chromosome."
    opts = ["--most-conserved", bed, "--score", "--seqname", seqname]
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
    emods = [estimate_base.with_suffix(x) for x in (".cons.mod", ".noncons.mod")]
    if emods[0].exists() and emods[1].exists():
        consEntropy(conf["target-coverage"], conf["expected-length"], *emods)
    kent.wigToBigWig(wig, chrom_sizes)
    if bw.exists():
        print(bw)
    return (bw, bed)


def prepare_mods(clade: Path) -> tuple[Path, Path]:
    target = clade.parent.name
    prepare_labeled_gff3(target)
    cons_mod = clade / "cons.mod"
    noncons_mod = clade / "noncons.mod"
    tree = phylo.get_subtree(clade.name.split("-"), phylo.shorten_names)
    cfutures: list[cli.FuturePath] = []
    nfutures: list[cli.FuturePath] = []
    pool = cli.ThreadPool()
    for chromosome in clade.glob("chromosome*"):
        maf = chromosome / "multiz.maf"
        gff = path_labeled_gff3(target, chromosome.name)
        cfutures.append(pool.submit(make_cons_mod, maf, gff, tree))
        nfutures.append(pool.submit(make_noncons_mod, maf, gff, tree))
    phyloBoot([f.result() for f in cfutures], cons_mod)
    phyloBoot([f.result() for f in nfutures], noncons_mod)
    return (cons_mod, noncons_mod)


def make_cons_mod(maf: Path, gff: Path, tree: str) -> Path:
    codons_ss = msa_view_features(maf, gff, conserved=True)
    codons_mods = phyloFit(codons_ss, tree, conserved=True)
    return most_conserved_mod(codons_mods)


def make_noncons_mod(maf: Path, gff: Path, tree: str) -> Path:
    _4d_codons_ss = msa_view_features(maf, gff, conserved=False)
    _4d_sites_ss = msa_view_ss(_4d_codons_ss)
    return phyloFit(_4d_sites_ss, tree, conserved=False)[0]


def msa_view_features(maf: Path, gff: Path, *, conserved: bool) -> Path:
    cmd = f"msa_view {maf!s} --in-format MAF --features <(gunzip -c {gff!s})"
    if conserved:
        outfile = maf.with_name("codons.ss")
        cmd += " --catmap 'NCATS = 3; CDS 1-3' --out-format SS --unordered-ss"
        cmd += " --reverse-groups transcript_id"
    else:
        outfile = maf.with_name("4d-codons.ss")
        cmd += " --4d"
    is_to_run = fs.is_outdated(outfile, maf)
    with subp.open_(outfile, "wb", if_=is_to_run) as fout:
        subp.run(
            cmd,
            if_=is_to_run,
            stdout=fout,
            shell=True,  # noqa: S604,
            executable="/bin/bash",
        )
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
        f"phyloFit --tree {tree} --msa-format SS {option}"
        f" --out-root {out_root} {ss!s}"
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
    cov = str(target_coverage)
    exp_len = str(expected_length)
    args = ["consEntropy", cov, exp_len, cons_mod, noncons_mod]
    outfile = cons_mod.with_name("entropy.txt")
    if_ = fs.is_outdated(outfile, [cons_mod, noncons_mod])
    with subp.open_(outfile, "wb", if_=if_) as fout:
        subp.run(args, stdout=fout, if_=if_)
    return outfile


def most_conserved_mod(mods: list[Path]) -> Path:
    outfile = mods[0].with_name("cons.mod")
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


def path_labeled_gff3(species: str, chromosome: str) -> Path:
    return Path("gff3") / species / f"labeled-{chromosome}.gff3.gz"


def prepare_labeled_gff3(species: str) -> None:
    """Deploy labeled copies of GFF3.

    src: {db.api.prefix}/gff3/{species}/*.{chromosome}.gff3.gz
    dst: ./gff3/{species}/labeled-{chromosome}.gff3.gz
    """
    shortname = phylo.shorten(species)
    for infile in api.list_chromosome_gff3(species):
        mobj = re.search(r"(chromosome.+)\.gff3\.gz$", infile.name)
        assert mobj, infile
        outfile = path_labeled_gff3(species, mobj.group(1))
        if fs.is_outdated(outfile, infile) and not cli.dry_run:
            _log.info(f"{outfile}")
            add_label_to_chr(infile, outfile, shortname + ".")


def add_label_to_chr(infile: Path, outfile: Path, label: str) -> None:
    """Modify GFF3 for msa_view.

    - Add species name to chromosome name, e.g., osat.1, zmay.2
    - Extract CDS
    """
    outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
    with gzip.open(infile, "rt") as fin, gzip.open(outfile, "wt") as fout:
        for line in fin:
            if not line.startswith("##"):
                break
            fout.write(line)
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(
            fout, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE
        )
        for row in reader:
            if len(row) < 8 or row[0].startswith("#"):  # noqa: PLR2004
                continue
            if row[2] == "CDS":
                row[0] = label + row[0]
                writer.writerow(row)


def clean(path: Path) -> None:
    it = itertools.chain(
        path.glob("*.mod"),
        path.glob("*.ss"),
        path.glob("chromosome*/*.mod"),
        path.glob("chromosome*/*.ss"),
        path.glob("chromosome*/phastcons.wig.gz"),
    )
    for file in it:
        print(file)
        if not cli.dry_run:
            file.unlink()


if __name__ == "__main__":
    main()
