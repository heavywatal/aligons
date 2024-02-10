import itertools
import logging
from pathlib import Path

from .db import api, phylo
from .extern import bedtools, htslib, kent, lastz, mafs2cram, multiz, phast
from .util import cli, fs, gff, log_config

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    tree = phylo.get_tree()
    parser = cli.ArgumentParser()
    parser.add_argument("-N", "--check-args", action="store_true")
    parser.add_argument("-t", "--tips", type=int, default=0)
    parser.add_argument("-g", "--max-bp", type=float, default=float("inf"))
    parser.add_argument("--compara", action="store_true")
    parser.add_argument("target", choices=phylo.extract_tip_names(tree))
    parser.add_argument("clade", choices=phylo.extract_inner_names(tree))
    args = parser.parse_args(argv or None)
    log_config()
    if args.check_args:
        return
    phastcons(args.target, args.clade, args.tips, args.max_bp, compara=args.compara)


def phastcons_block(bed: Path, clade: str) -> None:
    lst_species = phylo.list_species(clade)
    aln = lastz.PairwiseChromosomeAlignment(bed)
    not_in_aln = set(lst_species) - {aln.target, *aln.queries}
    assert not not_in_aln, not_in_aln
    ref_fa = api.genome_fa(aln.target)
    fts = [cli.thread_submit(mafs2cram.maf2cram, ft, ref_fa) for ft in aln.submit()]
    for future in fts:
        htslib.index(future.result())
    multiple = multiz.run(aln.target_dir, lst_species)
    if not (clade_alias := multiple.with_name(clade)).exists():
        clade_alias.symlink_to(multiple.name)
    (chrom_bws, chrom_beds) = phast.run(multiple)
    bigwig = kent.bigWigCat(multiple / "phastcons.bw", chrom_bws)
    mostcons_bed = bigwig.with_name("mostcons.bed.gz")
    htslib.concat_bgzip(chrom_beds, mostcons_bed)
    htslib.tabix(mostcons_bed)
    bed = bigwig.with_suffix(".bed.gz")
    if fs.is_outdated(bed, bigwig):
        htslib.bgzip(kent.bigWigToBed(bigwig), bed)
    cns_bed = bed.with_name("cns.bed.gz")
    if fs.is_outdated(cns_bed, bed):
        cds_df = gff.extract_cds_bed(api.genome_gff3(aln.target))
        cds = cds_df.write_csv(separator="\t", include_header=False).encode()
        cns = bedtools.subtract(bed, cds)
        cns = bedtools.remove_short(cns, 15)
        htslib.tabix(htslib.bgzip(cns, cns_bed))


def phastcons(
    target: str, clade: str, tips: int, max_bp: float, *, compara: bool
) -> None:
    lst_species = phylo.list_species(clade)
    lst_species = list(filter(lambda x: test_fasize(x, max_bp), lst_species))
    if compara:  # noqa: SIM108
        pairwise = Path("compara") / target
    else:
        pairwise = lastz.run(target, lst_species)
    fts = mafs2cram.run(pairwise, lst_species)
    n = tips or len(lst_species)
    for species in itertools.combinations(lst_species, n):
        if target not in species:
            continue
        _log.info(f"{species = }")
        multiple = multiz.run(pairwise, species)
        if n == len(lst_species):
            multiple.with_name(clade).symlink_to(multiple.name)
        chrom_bws = phast.run(multiple)
        kent.bigWigCat(multiple / "phastcons.bw", chrom_bws)
    cli.wait_raise(fts)


def test_fasize(species: str, max_bp: float) -> bool:
    bp = sum(x for x in api.chrom_sizes(species).values())
    ret = bp < max_bp
    _log.info(f"{species:30}{round(bp / 1e6):>5} Mbp {ret}")
    return ret


if __name__ == "__main__":
    main()
