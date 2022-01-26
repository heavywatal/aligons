"""Phylogenetic Analysis with Space/Time Models

src: ./multiple/{target}/{clade}/{chromosome}/multiz.maf
dst: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz

http://compgen.cshl.edu/phast/
"""

import logging
import gzip
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, IO

from . import cli
from .db import ensemblgenomes

_log = logging.getLogger(__name__)
_dry_run = False


def main(argv: list[str] = []):
    import argparse

    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("chromosome", nargs="+", type=Path)
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    global _dry_run
    _dry_run = args.dry_run
    _log.info("## msa_view, phyloFit, phyloBoot")
    (cons_mod, noncons_mod) = prepare(args.chromosome)
    _log.info("## phastCons")
    for dir in args.chromosome:
        phastCons(dir, cons_mod, noncons_mod)
    _log.info("## done")


def phastCons(path: Path, cons_mod: Path, noncons_mod: Path):
    maf = str(path / "multiz.maf")
    cmd = (
        "phastCons --target-coverage 0.25 --expected-length 12"
        f"--seqname {path.name} --msa-format MAF {maf} {cons_mod},{noncons_mod}"
    )
    wig = path / "phastcons.wig.gz"
    with gzip.open(wig, "w") as fout:
        run(cmd, stdout=fout)  # type: ignore
    return wig


def prepare(paths: list[Path]):
    cons_mods: list[str] = []
    _4d_sites_mods: list[str] = []
    for dir in paths:
        (cons, noncons) = make_mods(dir)
        cons_mods.append(str(cons))
        _4d_sites_mods.append(str(noncons))
    outdir = paths[0].parent
    cons_mod = outdir / "cons.mod"
    noncons_mod = outdir / "noncons.mod"
    phyloBoot(cons_mods, cons_mod)
    phyloBoot(_4d_sites_mods, noncons_mod)
    return (cons_mod, noncons_mod)


def make_mods(chromosome: Path):
    treefile = chromosome.parent / "tree.nh"
    target = chromosome.parent.parent.name
    gff = ensemblgenomes.get_file(f"*.{chromosome.name}.gff3.gz", target, "gff3")
    maf = chromosome / "multiz.maf"
    codons_ss = msa_view_features(maf, gff, True)
    _4d_codons_ss = msa_view_features(maf, gff, False)
    _4d_sites_ss = msa_view_ss(_4d_codons_ss)
    codons_mods = phyloFit(codons_ss, treefile, True)
    _4d_sites_mod = phyloFit(_4d_sites_ss, treefile, False)[0]
    cons_mod = most_conserved_mod(codons_mods)
    return (cons_mod, _4d_sites_mod)


def msa_view_features(maf: Path, gff: Path, conserved: bool):
    cmd = f"msa_view {str(maf)} --in-format MAF --features <(gunzip -c {str(gff)})"
    if conserved:
        outfile = maf.parent / "codons.ss"
        cmd += " --catmap 'NCATS = 3; CDS 1-3' --out-format SS --unordered-ss"
        cmd += " --reverse-groups transcript_id"
    else:
        outfile = maf.parent / "4d-codons.ss"
        cmd += " --4d"
    with open(outfile, "w") as fout:
        bash(cmd, stdout=fout)
    return outfile


def msa_view_ss(codons_ss: Path):
    outfile = codons_ss.parent / "4d-sites.ss"
    s = f"msa_view {str(codons_ss)} --in-format SS --out-format SS --tuple-size 1"
    with open(outfile, "w") as fout:
        run(s, stdout=fout)
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
    run(cmd)
    return outfiles


def phyloBoot(mods: list[str], outfile: Path):
    read_mods = ",".join(mods)
    run(f"phyloBoot --read-mods {read_mods} --output-average {outfile}")


def most_conserved_mod(mod_files: list[Path]):
    outfile = mod_files[0].parent / "cons.mod"
    if _dry_run:
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


def bash(cmd: str, stdout: IO[Any] | int | None = None):
    (_args, cmd) = cli.prepare_args(cmd, _dry_run)
    _log.info(cmd)
    return subprocess.run(cmd, stdout=stdout, shell=True, executable="/bin/bash")


def run(
    args: list[str] | str,
    stdin: IO[Any] | int | None = None,
    stdout: IO[Any] | int | None = None,
    shell: bool = False,
):  # kwargs hinders type inference to Popen[bytes]
    (args, cmd) = cli.prepare_args(args, _dry_run)
    _log.info(cmd)
    return subprocess.run(args, stdin=stdin, stdout=stdout, shell=shell)


if __name__ == "__main__":
    main()
