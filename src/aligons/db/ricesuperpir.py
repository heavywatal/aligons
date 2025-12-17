"""Rice Super Pan-genome Information Resource Database.

<http://www.ricesuperpir.com/>
"""

import logging
from pathlib import Path

from aligons.extern import lastz
from aligons.util import cli, fs, subp

from . import _rsrc, api, tools

_log = logging.getLogger(__name__)
species = "oryza_sativa"


def main(argv: list[str] | None = None) -> None:
    """CLI for downloading and processing NIP-T2T datasets."""
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-A", "--alignment", action="store_true")
    parser.add_argument("-Q", "--old-query", action="store_true")
    args = parser.parse_args(argv or None)
    if args.download:
        download()
    if args.alignment:
        align_dir = alignment(old_query=args.old_query)
        cat_chains(align_dir, old_query=args.old_query)
        return
    api.print_existing(species_label())


def prefix() -> Path:
    """Directory of preprocessed Rice Super PIR datasets."""
    return api.prefix("ricesuperpir")


def species_label() -> str:
    """Get species label for Rice Super PIR datasets."""
    return f"{prefix().name}/{species}"


def alignment(*, old_query: bool = False) -> Path:
    """Perform pairwise genome alignment between NIP-T2T and IRGSP-1.0.

    :param old_query: Use IRGSP-1.0 as query and NIP-T2T as target if True.
    :returns: Output directory for the pairwise alignment:
        `./pairwise/{target}/{query}/`
    """
    target = species
    query = species_label()
    if old_query:
        target, query = query, target
    return lastz.run(target, [query]) / query


def cat_chains(align_dir: Path, *, old_query: bool = False) -> Path:
    """Concatenate chain files from pairwise genome alignment.

    :param align_dir: Directory containing chain files.
    :param old_query: Use IRGSP-1.0 as query and NIP-T2T as target if True.
    :returns: Concatenated chain file path.
    """
    chains = list(fs.sorted_naturally(align_dir.glob("chr*/target.over.chain.gz")))
    target = "IRGSP-1.0"
    query = "NIP-T2T"
    if old_query:
        target, query = query, target
    stem = f"{target}To{query}"
    outfile = Path(f"{stem}.over.chain.gz")
    if_ = fs.is_outdated(outfile, chains)
    with subp.open_(outfile, "wb", if_=if_) as fout:
        subp.run(["cat", *chains], stdout=fout, if_=if_)
    return outfile


def download() -> None:
    """Download and preprocess NIP-T2T genome and annotation."""
    dataset: _rsrc.DataSet = {
        "url_prefix": "http://www.ricesuperpir.com/uploads/common",
        "species": "oryza_sativa",
        "version": "T2T-NIP",
        "draft": False,
        "label": "osat-t2t",
        "clade": "",
        "sequences": ["/genome_sequence/NIP-T2T.fa.gz"],
        "annotation": "/gene_annotation/NIP-T2T.gff3.gz",
    }
    fts: list[cli.Future[Path]] = []
    ft_fa, ft_gff = tools.fetch_and_bgzip(dataset, prefix())
    fts.append(ft_gff)
    masked = tools.softmask(ft_fa.result())
    fts.extend(tools.genome_to_twobits(masked))
    cli.wait_raise(fts)


if __name__ == "__main__":
    main()
