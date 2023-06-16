"""Phylogenetic Analysis with Space/Time Models.

src: ./multiple/{target}/{clade}/{chromosome}/multiz.maf
dst: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz

http://compgen.cshl.edu/phast/
"""
import csv
import gzip
import itertools
import logging
import os
import re
import sys
from pathlib import Path

from aligons.db import ensemblgenomes, phylo
from aligons.util import ConfDict, cli, config, empty_options, fs, read_config, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("clade", type=Path)  # multiple/{target}/{clade}
    args = parser.parse_args(argv or None)
    if args.config:
        read_config(args.config)
    if args.clean:
        clean(args.clade)
        return
    run(args.clade)


def run(path_clade: Path):
    (cons_mod, noncons_mod) = prepare_mods(path_clade)
    opts = config.get("phastCons", empty_options)
    pool = cli.ThreadPool()
    chrs = path_clade.glob("chromosome*")
    fts = [pool.submit(phastCons, d, cons_mod, noncons_mod, opts) for d in chrs]
    cli.wait_raise(fts)


def phastCons(  # noqa: N802
    path: Path, cons_mod: Path, noncons_mod: Path, options: ConfDict = empty_options
):
    maf = str(path / "multiz.maf")
    seqname = path.name.split(".", 1)[1]  # remove "chromosome."
    cmd = "phastCons"
    cmd += subp.optjoin(options)
    cmd += f" --seqname {seqname} --msa-format MAF {maf} {cons_mod},{noncons_mod}"
    wig = path / "phastcons.wig.gz"
    is_to_run = fs.is_outdated(wig, [cons_mod, noncons_mod])
    p = subp.run(cmd, if_=is_to_run, stdout=subp.PIPE)
    with gzip.open(devnull_if(not is_to_run, wig), "wb") as fout:
        fout.write(p.stdout)
    if wig.exists():
        print(wig)
    return wig


def prepare_mods(clade: Path):
    target = clade.parent.name
    prepare_labeled_gff3(target)
    cons_mod = clade / "cons.mod"
    noncons_mod = clade / "noncons.mod"
    tree = phylo.get_newick(clade.name.split("-"), phylo.shorten_names)
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


def make_cons_mod(maf: Path, gff: Path, tree: str):
    codons_ss = msa_view_features(maf, gff, conserved=True)
    codons_mods = phyloFit(codons_ss, tree, conserved=True)
    return most_conserved_mod(codons_mods)


def make_noncons_mod(maf: Path, gff: Path, tree: str):
    _4d_codons_ss = msa_view_features(maf, gff, conserved=False)
    _4d_sites_ss = msa_view_ss(_4d_codons_ss)
    return phyloFit(_4d_sites_ss, tree, conserved=False)[0]


def msa_view_features(maf: Path, gff: Path, *, conserved: bool):
    cmd = f"msa_view {maf!s} --in-format MAF --features <(gunzip -c {gff!s})"
    if conserved:
        outfile = maf.parent / "codons.ss"
        cmd += " --catmap 'NCATS = 3; CDS 1-3' --out-format SS --unordered-ss"
        cmd += " --reverse-groups transcript_id"
    else:
        outfile = maf.parent / "4d-codons.ss"
        cmd += " --4d"
    is_to_run = fs.is_outdated(outfile, maf)
    p = subp.run(
        cmd, if_=is_to_run, stdout=subp.PIPE, shell=True, executable="/bin/bash"
    )
    with devnull_if(not is_to_run, outfile).open("wb") as fout:
        fout.write(p.stdout)
    return outfile


def msa_view_ss(codons_ss: Path):
    outfile = codons_ss.parent / "4d-sites.ss"
    s = f"msa_view {codons_ss!s} --in-format SS --out-format SS --tuple-size 1"
    is_to_run = fs.is_outdated(outfile, codons_ss)
    p = subp.run(s, if_=is_to_run, stdout=subp.PIPE)
    with devnull_if(not is_to_run, outfile).open("wb") as fout:
        fout.write(p.stdout)
    return outfile


def phyloFit(ss: Path, tree: str, *, conserved: bool):  # noqa: N802
    if conserved:
        out_root = str(ss.parent / "codons")
        outfiles = [Path(f"{out_root}.{i}.mod") for i in range(1, 4)]
        option = "--do-cats 1,2,3"
    else:
        out_root = str(ss.parent / "4d-sites")
        outfiles = [Path(f"{out_root}.mod")]
        option = ""
    cmd = (
        f"phyloFit --tree {tree} --msa-format SS {option}"
        f" --out-root {out_root} {ss!s}"
    )
    subp.run(cmd, if_=fs.is_outdated(outfiles[0], ss))
    return outfiles


def phyloBoot(mods: list[Path], outfile: Path):  # noqa: N802
    read_mods = ",".join(str(x) for x in mods)
    subp.run(
        f"phyloBoot --read-mods {read_mods} --output-average {outfile}",
        if_=fs.is_outdated(outfile, mods),
    )


def most_conserved_mod(mods: list[Path]):
    outfile = mods[0].parent / "cons.mod"
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


def extract_tree(mod: str):
    assert (m := re.search(r"TREE: (.+;)", mod)), mod
    return m.group(1)


def path_labeled_gff3(species: str, chromosome: str):
    return Path("gff3") / species / f"labeled-{chromosome}.gff3.gz"


def prepare_labeled_gff3(species: str):
    """Deploy labeled copies of GFF3.

    src: {ensemblgenomes.prefix}/gff3/{species}/*.{chromosome}.gff3.gz
    dst: ./gff3/{species}/labeled-{chromosome}.gff3.gz
    """
    shortname = phylo.shorten(species)
    for infile in ensemblgenomes.glob("*.chromosome*.gff3.gz", [species]):
        mobj = re.search(r"(chromosome.+)\.gff3\.gz$", infile.name)
        assert mobj
        outfile = path_labeled_gff3(species, mobj.group(1))
        if fs.is_outdated(outfile, infile) and not cli.dry_run:
            _log.info(f"{outfile}")
            add_label_to_chr(infile, outfile, shortname + ".")


def add_label_to_chr(infile: Path, outfile: Path, label: str):
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


def clean(path: Path):
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


def devnull_if(cond: bool, file: Path):  # noqa: FBT001
    if cond:
        return Path(os.devnull)
    return file


if __name__ == "__main__":
    main()
