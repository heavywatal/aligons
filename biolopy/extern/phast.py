"""Phylogenetic Analysis with Space/Time Models

src: ./multiple/{target}/{clade}/{chromosome}/multiz.maf
dst: ./multiple/{target}/{clade}/{chromosome}/phastcons.wig.gz

http://compgen.cshl.edu/phast/
"""
import concurrent.futures as confu
import csv
import gzip
import itertools
import logging
import os
import re
import sys
from pathlib import Path
from typing import IO, AnyStr, cast

from ..db import ensemblgenomes, phylo
from ..util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = cli.logging_argparser()
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("clade", type=Path)  # multiple/{target}/{clade}
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.clean:
        clean(args.clade)
        return
    run(args.clade, args.jobs)


def run(path_clade: Path, jobs: int):
    (cons_mod, noncons_mod) = prepare_mods(path_clade, jobs)
    with confu.ThreadPoolExecutor(max_workers=jobs) as pool:
        chrs = path_clade.glob("chromosome*")
        futures = [pool.submit(phastCons, d, cons_mod, noncons_mod) for d in chrs]
        for future in confu.as_completed(futures):
            if (wig := future.result()).exists():
                print(wig)


def phastCons(path: Path, cons_mod: Path, noncons_mod: Path):
    maf = str(path / "multiz.maf")
    seqname = path.name.split(".", 1)[1]  # remove "chromosome."
    cmd = (
        "phastCons --target-coverage 0.25 --expected-length 12"
        f" --seqname {seqname} --msa-format MAF {maf} {cons_mod},{noncons_mod}"
    )
    wig = path / "phastcons.wig.gz"
    is_to_run = fs.is_outdated(wig, [cons_mod, noncons_mod])
    p = subp.run_if(is_to_run, cmd, stdout=subp.PIPE)
    with open_if(is_to_run, wig, "wb") as fout:
        fout.write(p.stdout)
    return wig


def prepare_mods(clade: Path, jobs: int):
    target = clade.parent.name
    prepare_labeled_gff3(target)
    cons_mod = clade / "cons.mod"
    noncons_mod = clade / "noncons.mod"
    tree = phylo.shorten_labels(phylo.trees[clade.name])
    cfutures: list[confu.Future[Path]] = []
    nfutures: list[confu.Future[Path]] = []
    with confu.ThreadPoolExecutor(max_workers=jobs) as pool:
        for chromosome in clade.glob("chromosome*"):
            maf = chromosome / "multiz.maf"
            gff = path_labeled_gff3(target, chromosome.name)
            cfutures.append(pool.submit(make_cons_mod, maf, gff, tree))
            nfutures.append(pool.submit(make_noncons_mod, maf, gff, tree))
    phyloBoot([f.result() for f in cfutures], cons_mod)
    phyloBoot([f.result() for f in nfutures], noncons_mod)
    return (cons_mod, noncons_mod)


def make_cons_mod(maf: Path, gff: Path, tree: str):
    codons_ss = msa_view_features(maf, gff, True)
    codons_mods = phyloFit(codons_ss, tree, True)
    return most_conserved_mod(codons_mods)


def make_noncons_mod(maf: Path, gff: Path, tree: str):
    _4d_codons_ss = msa_view_features(maf, gff, False)
    _4d_sites_ss = msa_view_ss(_4d_codons_ss)
    return phyloFit(_4d_sites_ss, tree, False)[0]


def msa_view_features(maf: Path, gff: Path, conserved: bool):
    cmd = f"msa_view {str(maf)} --in-format MAF --features <(gunzip -c {str(gff)})"
    if conserved:
        outfile = maf.parent / "codons.ss"
        cmd += " --catmap 'NCATS = 3; CDS 1-3' --out-format SS --unordered-ss"
        cmd += " --reverse-groups transcript_id"
    else:
        outfile = maf.parent / "4d-codons.ss"
        cmd += " --4d"
    is_to_run = fs.is_outdated(outfile, maf)
    p = subp.run_if(
        is_to_run, cmd, stdout=subp.PIPE, shell=True, executable="/bin/bash"
    )
    with open_if(is_to_run, outfile, "wb") as fout:
        fout.write(p.stdout)
    return outfile


def msa_view_ss(codons_ss: Path):
    outfile = codons_ss.parent / "4d-sites.ss"
    s = f"msa_view {str(codons_ss)} --in-format SS --out-format SS --tuple-size 1"
    is_to_run = fs.is_outdated(outfile, codons_ss)
    p = subp.run_if(is_to_run, s, stdout=subp.PIPE)
    with open_if(is_to_run, outfile, "wb") as fout:
        fout.write(p.stdout)
    return outfile


def phyloFit(ss: Path, tree: str, conserved: bool):
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
        f" --out-root {out_root} {str(ss)}"
    )
    subp.run_if(fs.is_outdated(outfiles[0], ss), cmd)
    return outfiles


def phyloBoot(mods: list[Path], outfile: Path):
    read_mods = ",".join(str(x) for x in mods)
    subp.run_if(
        fs.is_outdated(outfile, mods),
        f"phyloBoot --read-mods {read_mods} --output-average {outfile}",
    )


def most_conserved_mod(mods: list[Path]):
    outfile = mods[0].parent / "cons.mod"
    if not fs.is_outdated(outfile, mods) or cli.dry_run:
        return outfile
    shortest_length = sys.float_info.max
    conserved = ""
    for mod in mods:
        with open(mod, "r") as fin:
            content = fin.read()
            lengths = phylo.extract_lengths(extract_tree(content))
            total = sum(lengths)
            _log.debug(f"{mod}: {total}")
            if total < shortest_length:
                shortest_length = total
                conserved = content
    _log.debug(f"{shortest_length=}")
    with open(outfile, "w") as fout:
        fout.write(conserved)
    return outfile


def extract_tree(mod: str):
    assert (m := re.search(r"TREE: (.+;)", mod)), mod
    return m.group(1)


def path_labeled_gff3(species: str, chromosome: str):
    return Path("gff3") / species / f"labeled-{chromosome}.gff3.gz"


def prepare_labeled_gff3(species: str):
    """
    src: {ensemblgenomes.prefix}/gff3/{species}/*.{chromosome}.gff3.gz
    dst: ./gff3/{species}/labeled-{chromosome}.gff3.gz
    """
    shortname = phylo.shorten(species)
    for infile in ensemblgenomes.rglob("*.chromosome*.gff3.gz", [species]):
        mobj = re.search(r"(chromosome.+)\.gff3\.gz$", infile.name)
        assert mobj
        outfile = path_labeled_gff3(species, mobj.group(1))
        if fs.is_outdated(outfile, infile) and not cli.dry_run:
            _log.info(f"{outfile}")
            add_label_to_chr(infile, outfile, shortname + ".")


def add_label_to_chr(infile: Path, outfile: Path, label: str):
    """Modify GFF3 for msa_view

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
            if len(row) < 8 or row[0].startswith("#"):
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


def open_if(cond: bool, file: Path, mode: str = "r") -> IO[AnyStr]:
    suffix = file.suffix
    if cli.dry_run or not cond:
        file = Path(os.devnull)
    if suffix == ".gz":
        f = cast(IO[AnyStr], gzip.open(file, mode))
    else:
        f = open(file, mode)
    return f


if __name__ == "__main__":
    main()
