import concurrent.futures as confu
import logging
import os
import re
from pathlib import Path

from ..extern import htslib, kent
from ..util import cli, fs
from . import ensemblgenomes, phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-C", "--compara")
    parser.add_argument("-c", "--clade", default="bep")
    args = parser.parse_args(argv or None)
    cli.dry_run = args.dry_run
    cli.logging_config(args.loglevel)
    if args.compara:
        with ensemblgenomes.FTPensemblgenomes() as ftp:
            dirs = ftp.download_maf(args.compara)
        for dir in dirs:
            if dir.is_dir():
                print(dir)
                use_compara(dir)
        return
    tree = phylo.newicks[args.clade]
    species = phylo.extract_names(tree)
    if args.download:
        download(species, jobs=args.jobs)
    else:
        index(species, jobs=args.jobs)


def use_compara(dir: Path):
    mobj = re.search(r"([^_]+)_.+?\.v\.([^_]+)", dir.name)
    assert mobj
    target_short = mobj.group(1)
    target = list(ensemblgenomes.expand_shortnames([target_short]))[0]
    query_short = mobj.group(2)
    query = list(ensemblgenomes.expand_shortnames([query_short]))[0]
    pat = re.compile(r"lastz_net\.([^_]+)_(\d+)\.maf$")
    up_to_date: set[Path] = set()
    for maf in fs.sorted_naturally(dir.glob("*_*.maf")):
        mobj = pat.search(maf.name)
        assert mobj
        seq = mobj.group(1)
        serial = mobj.group(2)
        outdir = Path(f"compara/{target}/{query}/chromosome.{seq}")
        outfile = outdir / "sing.maf"
        lines: list[str] = []
        if serial == "1":
            _log.info(str(outfile))
            if not fs.is_outdated(outfile, maf):
                up_to_date.add(outfile)
            mode = "w"
            lines.append("##maf version=1 scoring=LASTZ_NET\n")
            outdir.mkdir(0o755, parents=True, exist_ok=True)
        else:
            mode = "a"
        if outfile in up_to_date:
            continue
        with open(maf, "r") as fin:
            for line in fin:
                if line.startswith("#") or line.startswith("a#"):
                    continue
                elif line.startswith(" score"):
                    line = "a" + line
                elif line.startswith("s"):
                    line = line.replace(target, target_short)
                    line = line.replace(query, query_short)
                    # TODO: multiz inconsistent row size
                lines.append(line)
        with open(outfile, mode) as fout:
            fout.writelines(lines)


def download(species: list[str], jobs: int | None = os.cpu_count()):
    assert species
    assert not (d := set(species) - set(ensemblgenomes.species_names_all())), d
    with (
        confu.ThreadPoolExecutor(max_workers=jobs) as pool,
        ensemblgenomes.FTPensemblgenomes() as ftp,
    ):
        for sp in species:
            dir = ftp.download_fasta(sp)
            pool.submit(index_fasta, dir)
            dir = ftp.download_gff3(sp)
            pool.submit(index_gff3, dir)
    # for sp in species:
    #     options = "--include *_sm.chromosome.*.fa.gz --exclude *.gz"
    #     ensemblgenomes.rsync(f"fasta/{sp}/dna", options)
    #     options = "--include *.chromosome.*.gff3.gz --exclude *.gz"
    #     ensemblgenomes.rsync(f"gff3/{sp}", options)


def index(species: list[str] = [], jobs: int | None = os.cpu_count()):
    with confu.ThreadPoolExecutor(max_workers=jobs) as pool:
        pool.map(index_fasta, ensemblgenomes.species_dirs("fasta", species))
        pool.map(index_gff3, ensemblgenomes.species_dirs("gff3", species))


def index_fasta(path: Path):  # fasta/{species}
    """Create bgzipped and indexed genome.fa."""
    if path.name != "dna":
        path /= "dna"
    for chromosome in fs.sorted_naturally(path.glob(r"*.chromosome.*.fa.gz")):
        _log.info(str(kent.faToTwoBit(chromosome)))
    genome = htslib.create_genome_bgzip(path)
    _log.info(str(genome))
    _log.info(str(htslib.faidx(genome)))
    _log.info(str(kent.faToTwoBit(genome)))
    _log.info(str(kent.faSize(genome)))
    print(path)
    return genome


def index_gff3(path: Path):  # gff3/{species}
    """Create bgzipped and indexed genome.gff3."""
    genome = htslib.create_genome_bgzip(path)
    _log.info(str(genome))
    _log.info(str(htslib.tabix(genome)))
    print(path)
    return genome


if __name__ == "__main__":
    main()
