"""Phylogenetic Analysis with Space/Time Models

src: ./multiple/{target}/{clade}/{chromosome}/multiz.maf
dst: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz

http://compgen.cshl.edu/phast/
"""
import argparse
import concurrent.futures as confu
import csv
import gzip
import itertools
import logging
import os
import re
import sys
from pathlib import Path
from subprocess import PIPE
from typing import IO, AnyStr, cast

from . import cli
from .db import ensemblgenomes, name

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("clade", type=Path)
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.clean:
        clean(args.clade)
        return
    (cons_mod, noncons_mod) = prepare_mods(args.clade, args.jobs)
    with confu.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        chrs = args.clade.glob("chromosome*")
        futures = [pool.submit(phastCons, d, cons_mod, noncons_mod) for d in chrs]
        for future in confu.as_completed(futures):
            print(future.result())


def phastCons(path: Path, cons_mod: Path, noncons_mod: Path):
    maf = str(path / "multiz.maf")
    seqname = path.name.split(".", 1)[1]  # remove "chromosome."
    cmd = (
        "phastCons --target-coverage 0.25 --expected-length 12"
        f" --seqname {seqname} --msa-format MAF {maf} {cons_mod},{noncons_mod}"
    )
    wig = path / "phastcons.wig.gz"
    p = cli.run(cmd, stdout=PIPE)
    with open_if_not_dry_run(wig, "wb") as fout:
        fout.write(p.stdout)
    return wig


def prepare_mods(clade: Path, jobs: int):
    cons_mod = clade / "cons.mod"
    noncons_mod = clade / "noncons.mod"
    if cons_mod.exists() and noncons_mod.exists():
        return (cons_mod, noncons_mod)
    treefile = clade / "tree.nh"
    target = clade.parent.name
    prepare_labeled_gff3(target)
    cfutures: list[confu.Future[Path]] = []
    nfutures: list[confu.Future[Path]] = []
    with confu.ThreadPoolExecutor(max_workers=jobs) as pool:
        for chromosome in clade.glob("chromosome*"):
            maf = chromosome / "multiz.maf"
            gff = labeled_gff3(target, chromosome.name)
            cfutures.append(pool.submit(make_cons_mod, maf, gff, treefile))
            nfutures.append(pool.submit(make_noncons_mod, maf, gff, treefile))
    phyloBoot([str(f.result()) for f in cfutures], cons_mod)
    phyloBoot([str(f.result()) for f in nfutures], noncons_mod)
    return (cons_mod, noncons_mod)


def make_cons_mod(maf: Path, gff: Path, treefile: Path):
    codons_ss = msa_view_features(maf, gff, True)
    codons_mods = phyloFit(codons_ss, treefile, True)
    return most_conserved_mod(codons_mods)


def make_noncons_mod(maf: Path, gff: Path, treefile: Path):
    _4d_codons_ss = msa_view_features(maf, gff, False)
    _4d_sites_ss = msa_view_ss(_4d_codons_ss)
    return phyloFit(_4d_sites_ss, treefile, False)[0]


def msa_view_features(maf: Path, gff: Path, conserved: bool):
    cmd = f"msa_view {str(maf)} --in-format MAF --features <(gunzip -c {str(gff)})"
    if conserved:
        outfile = maf.parent / "codons.ss"
        cmd += " --catmap 'NCATS = 3; CDS 1-3' --out-format SS --unordered-ss"
        cmd += " --reverse-groups transcript_id"
    else:
        outfile = maf.parent / "4d-codons.ss"
        cmd += " --4d"
    with open_if_not_dry_run(outfile, "w") as fout:
        cli.run(cmd, stdout=fout, shell=True, executable="/bin/bash")
    return outfile


def msa_view_ss(codons_ss: Path):
    outfile = codons_ss.parent / "4d-sites.ss"
    s = f"msa_view {str(codons_ss)} --in-format SS --out-format SS --tuple-size 1"
    p = cli.run(s, stdout=PIPE)
    with open_if_not_dry_run(outfile, "wb") as fout:
        fout.write(p.stdout)
    return outfile


def phyloFit(ss: Path, treefile: Path, conserved: bool):
    if conserved:
        out_root = str(ss.parent / "codons")
        outfiles = [Path(f"{out_root}.{i}.mod") for i in range(1, 4)]
        option = "--do-cats 1,2,3"
    else:
        out_root = str(ss.parent / "4d-sites")
        outfiles = [Path(f"{out_root}.mod")]
        option = ""
    cmd = (
        f"phyloFit --tree {str(treefile)} --msa-format SS {option}"
        f" --out-root {str(out_root)} {str(ss)}"
    )
    cli.run(cmd)
    return outfiles


def phyloBoot(mods: list[str], outfile: Path):
    read_mods = ",".join(mods)
    cli.run(f"phyloBoot --read-mods {read_mods} --output-average {outfile}")


def most_conserved_mod(mod_files: list[Path]):
    outfile = mod_files[0].parent / "cons.mod"
    if cli.dry_run:
        return outfile
    shortest_length = sys.float_info.max
    conserved = ""
    for f in mod_files:
        with open(f, "r") as fin:
            content = fin.read()
            tree = extract_tree(content)
            assert tree
            lengths = branch_lengths(tree)
            assert lengths
            total = sum(lengths)
            if total < shortest_length:
                shortest_length = total
                conserved = content
    with open(outfile, "w") as fout:
        fout.write(conserved)
    return outfile


def extract_tree(content: str):
    if m := re.search(r"TREE: (.+);", content):
        return m.group(1)


def branch_lengths(newick: str):
    return [float(m.group(0)) for m in re.finditer(r"[\d.]+", newick)]


def labeled_gff3(species: str, chromosome: str):
    return Path("gff3") / species / f"labeled-{chromosome}.gff3.gz"


def prepare_labeled_gff3(species: str):
    """
    src: {ensemblgenomes.prefix}/gff3/{species}/*.{chromosome}.gff3.gz
    dst: ./gff3/{species}/labeled-{chromosome}.gff3.gz
    """
    shortname = name.shorten(species)
    for infile in ensemblgenomes.rglob("*.chromosome*.gff3.gz", species, "gff3"):
        mobj = re.search(r"(chromosome.+)\.gff3\.gz$", infile.name)
        assert mobj
        outfile = labeled_gff3(species, mobj.group(1))
        if not outfile.exists() and not cli.dry_run:
            _log.info(f"{outfile}")
            add_label_to_chr(infile, outfile, shortname + ".")


def add_label_to_chr(infile: Path, outfile: Path, label: str):
    """Modify GFF3 for msa_view

    - Add species name to chromosome name, e.g., osat.1, zmay.2
    - Extract CDS
    """
    if not cli.dry_run:
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
    with gzip.open(infile, "rt") as fin, gzip.open(outfile, "wt") as fout:
        reader = csv.reader(fin, delimiter="\t")
        for row in reader:
            if len(row) < 8 or row[0].startswith("#"):
                continue
            if row[2] == "CDS":
                fout.write(label + "\t".join(row) + "\n")


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


def open_if_not_dry_run(file: Path, mode: str = "r") -> IO[AnyStr]:
    suffix = file.suffix
    if cli.dry_run:
        file = Path("/dev/null")
    if suffix == ".gz":
        f = cast(IO[AnyStr], gzip.open(file, mode))
    else:
        f = open(file, mode)
    return f


if __name__ == "__main__":
    main()
