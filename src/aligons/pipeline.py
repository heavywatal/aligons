import itertools
import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

from .db import api, phylo
from .extern import bedtools, htslib, kent, lastz, mafs2cram, multiz, phast
from .util import cli, fs, log_config

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    tree = phylo.get_tree()
    parser = cli.ArgumentParser()
    parser.add_argument("-N", "--check-args", action="store_true")
    parser.add_argument("-t", "--tips", type=int, default=0)
    parser.add_argument("-g", "--max-bp", type=float, default=float("inf"))
    parser.add_argument("--compara", action="store_true")
    parser.add_argument("--bed", type=Path)
    parser.add_argument("target", choices=phylo.extract_tip_names(tree))
    parser.add_argument("clade", choices=phylo.extract_inner_names(tree))
    args = parser.parse_args(argv or None)
    log_config()
    if args.check_args:
        return
    if args.bed:
        one_block(args.bed, args.clade)
        return
    genome_wide(args.target, args.clade, args.tips, args.max_bp, compara=args.compara)


def one_block(bed: Path, clade: str) -> None:
    lst_species = phylo.list_species(clade)
    aln = lastz.PairwiseChromosomeAlignment(bed, lst_species)
    not_in_aln = set(lst_species) - {aln.target, *aln.queries}
    assert not not_in_aln, not_in_aln
    ref_fa = api.genome_fa(aln.target)
    fts = [cli.thread_submit(mafs2cram.maf2cram, ft, ref_fa) for ft in aln.submit()]
    for future in fts:
        htslib.index(future.result())
    multiz_phast(aln.target_dir, lst_species, clade)


def genome_wide(
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
        multiz_phast(pairwise, species, clade)
    cli.wait_raise(fts)


def multiz_phast(pairwise: Path, lst_species: Sequence[str], clade: str) -> Path:
    target = pairwise.name
    multiple = multiz.run(pairwise, lst_species)
    is_clade = set(lst_species) == set(phylo.list_species(clade))
    if is_clade:
        clade_alias = multiple.with_name(clade)
        fs.symlink(multiple, clade_alias, relative=True)
    (bigwig, mostcons) = phast.run(multiple)
    cds = phast.cds_gff3(target)
    cns_bed = mostcons.with_name("cns.bed.gz")
    _subtract_cds(mostcons, cds, cns_bed)
    cns0_bed = bigwig.with_name("cns0.bed.gz")
    _subtract_cds(bigwig, cds, cns0_bed)
    if is_clade:
        clade_alias = mostcons.parent.with_name(clade)
        fs.symlink(mostcons.parent, clade_alias, relative=True)
    return cns_bed


def _subtract_cds(bed: Path, cds: Path, outfile: Path) -> Path:
    if bed.suffix == ".bw":
        bigwig = bed
        bed = bigwig.with_suffix(".bed.gz")
        if fs.is_outdated(bed, bigwig):
            htslib.bgzip(kent.bigWigToBed(bigwig), bed)
    if fs.is_outdated(outfile, bed):
        cns = bedtools.subtract(bed, cds)
        cns = bedtools.remove_short(cns, 15)
        htslib.tabix(htslib.bgzip(cns, outfile))
    return outfile


def test_fasize(species: str, max_bp: float) -> bool:
    bp = sum(x for x in api.chrom_sizes(species).values())
    ret = bp < max_bp
    _log.info(f"{species:30}{round(bp / 1e6):>5} Mbp {ret}")
    return ret


if __name__ == "__main__":
    main()
