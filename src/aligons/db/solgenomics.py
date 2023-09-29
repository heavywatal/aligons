"""https://solgenomics.net/ftp/genomes/."""
import logging

from aligons import db
from aligons.db import api
from aligons.util import cli

from . import phylo, tools

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-M", "--mask", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    args = parser.parse_args(argv or None)
    _test_newick()
    if args.download:
        fts = download()
        cli.wait_raise(fts)
    if args.mask:
        fts = split_mask_index()
        cli.wait_raise(fts)


def download():
    fts: list[cli.FuturePath] = []
    for entry in tools.iter_dataset("solgenomics.toml"):
        entry["url_prefix"] = "https://solgenomics.net/ftp/genomes/"
        fts.extend(tools.retrieve(entry, db_prefix()))
    return fts


def split_mask_index():
    fts: list[cli.FuturePath] = []
    for entry in tools.iter_dataset("solgenomics.toml"):
        species = entry["species"]
        fts.append(tools.prepare_fasta(species))
        gff3_gz = api.get_file_nolabel("*.gff3.gz", species)
        fts.append(cli.thread_submit(tools.index_gff3, [gff3_gz]))
    return fts


def db_prefix():
    return db.path("solgenomics")


def _test_newick():
    for entry in tools.iter_dataset("solgenomics.toml"):
        species = entry["species"]
        if species not in phylo.list_species():
            _log.warning(f"{species} not found in phylo")


if __name__ == "__main__":
    main()
