import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

from aligons.util import cli

from . import _rsrc, api, ensemblgenomes, phylo, tools

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("--compara", choices=phylo.list_species())
    parser.add_argument("-N", "--naro", action="store_true")
    parser.add_argument("-C", "--clade", default="bep")
    args = parser.parse_args(argv or None)
    if args.compara:
        ensemblgenomes.download_compara(args.compara)
        return
    if args.naro:
        cli.wait_raise(shumari_v2())
        return
    species = phylo.list_species(args.clade)
    ensemblgenomes.download_via_ftp(species)


def shumari_v2() -> tuple[cli.Future[Path], cli.Future[Path]]:
    prefix = api.prefix("naro")
    src_root = _rsrc.db_root("naro.go.jp")
    src = src_root / "shumari_v2"
    in_fa = src / "Vangularis_v2.fa"
    in_gff3 = src / "Vangularis_v2.a1.genes.gff"
    species = "vigna_angularis"
    version = "naito_v2"
    stem = f"{species}_{version}"
    (prefix / species).mkdir(0o755, parents=True, exist_ok=True)
    out_fa = prefix / species / f"{stem}.dna.genome.fa.gz"
    out_gff3 = prefix / species / f"{stem}.genome.gff3.gz"
    ft_fa = cli.thread_submit(tools.index_bgzip, in_fa, out_fa)
    ft_gff3 = cli.thread_submit(tools.index_bgzip, in_gff3, out_gff3)
    ft_masked_fa = tools.softmask(ft_fa.result(), "fabaceae")
    cli.wait_raise(tools.genome_to_twobits(ft_masked_fa))
    return ft_masked_fa, ft_gff3


if __name__ == "__main__":
    main()
