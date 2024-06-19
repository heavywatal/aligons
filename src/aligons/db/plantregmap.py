"""http://plantregmap.gao-lab.org/."""

import logging
import re
from collections.abc import Iterator
from pathlib import Path

from aligons.extern import htslib, kent, mafs2cram
from aligons.util import cli, dl, fs, subp, tomli_w, tomllib

from . import _rsrc, api, tools

_log = logging.getLogger(__name__)
_HOST = "plantregmap.gao-lab.org"
_FTP_HOST = "ftp.cbi.pku.edu.cn"


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-F", "--ftp", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-G", "--genome", action="store_true")
    parser.add_argument("-C", "--to-cram", type=Path)
    parser.add_argument("pattern", nargs="?", default="*")
    args = parser.parse_args(argv or None)
    if args.ftp:
        fts = download_via_ftp()
        cli.wait_raise(fts)
    if args.download:
        fts = [retrieve_deploy(q) for q in iter_download_queries()]
        cli.wait_raise(fts)
    if args.genome:
        fts_fa: list[cli.Future[Path]] = []
        fts_gff: list[cli.Future[Path]] = []
        for ft_fa, ft_gff in iter_fetch_and_bgzip():
            fts_fa.append(ft_fa)
            fts_gff.append(ft_gff)
        fts: list[cli.Future[Path]] = []
        for ft in cli.as_completed(fts_fa):
            masked = tools.softmask(ft.result())
            fts.extend(tools.genome_to_twobits(masked))
        cli.wait_raise(fts)
        cli.wait_raise(fts_gff)
    if args.to_cram:
        maf = net_to_maf(args.to_cram)
        to_cram(maf)


def iter_fetch_and_bgzip() -> Iterator[tuple[cli.Future[Path], cli.Future[Path]]]:
    for entry in _rsrc.iter_builtin_dataset("plantregmap.toml"):
        yield tools.fetch_and_bgzip(entry, db_prefix())


def retrieve_deploy(query: str) -> cli.Future[Path]:
    url = f"http://{_HOST}/download_ftp.php?{query}"
    relpath = query.split("=", 1)[1]
    rawfile = _prefix_mirror() / relpath
    outfile = db_prefix() / relpath.removeprefix("08-download/")
    if outfile.suffix in (".bed", ".gff"):
        outfile = outfile.with_suffix(outfile.suffix + ".gz")
    elif outfile.name.endswith(".gtf.gz"):
        outfile = outfile.with_suffix("").with_suffix(".gff.gz")
    response = dl.fetch(url, rawfile)
    return cli.thread_submit(tools.bgzip_or_symlink, response, outfile)


def iter_download_queries() -> Iterator[str]:
    targets = [
        "Oryza_sativa_Japonica_Group",
        "Arabidopsis_thaliana",
        "Solanum_lycopersicum",
    ]
    pattern = "|".join(targets)
    for query in iter_download_queries_all():
        if re.search(pattern, query):
            yield query


def iter_download_queries_all() -> Iterator[str]:
    content = download_php()
    for mobj in re.finditer(r"download_ftp\.php\?([^\"']+)", content):
        yield mobj[1]


def download_php() -> str:
    url = f"http://{_HOST}/download.php"
    cache = _prefix_mirror() / "download.php.html"
    return dl.fetch(url, cache).text_force


def db_prefix() -> Path:
    return api.prefix("plantregmap")


def _prefix_mirror() -> Path:
    return _rsrc.db_root(_FTP_HOST) / "plantregmap"


def rglob(pattern: str, species: str = ".") -> Iterator[Path]:
    for species_dir in db_prefix().iterdir():
        if re.search(species, species_dir.name, re.IGNORECASE):
            yield from species_dir.rglob(pattern)


def net_to_maf(net_gz: Path) -> Path:
    stem = net_gz.with_suffix("").stem
    target, query = extract_species(stem)
    sing_maf = net_gz.parent / (stem + ".maf")
    chain_gz = sing_maf.with_suffix(".chain.gz")
    return kent.net_to_maf(net_gz, chain_gz, sing_maf, target, query)


def to_cram(maf: Path) -> Path:
    target, _query = extract_species(maf.stem)
    reference = api.genome_fa(target)
    outfile = mafs2cram.maf2cram(maf, reference)
    htslib.index(outfile)
    return outfile


def download_via_ftp() -> list[cli.Future[Path]]:
    targets = [
        "Oryza_sativa_Japonica_Group",
        "Solanum_lycopersicum",
    ]
    fts: list[cli.Future[Path]] = []
    with FTPplantregmap() as ftp:
        for species in targets:
            fts.extend(ftp.download(species))
    return fts


def extract_species(stem: str) -> tuple[str, str]:
    short_target, short_query = stem.split("_")
    target = lengthen(short_target)
    query = lengthen(short_query)
    return (target, query)


def lengthen(species: str) -> str:
    try:
        return lengthen.abbr[species]  # pyright: ignore[reportFunctionMemberAccess]
    except AttributeError:
        rev_abbr = {v: k for k, v in species_abbr().items()}
        rev_abbr["Osj"] = "Oryza_sativa_Japonica_Group"
        lengthen.abbr = rev_abbr  # pyright: ignore[reportFunctionMemberAccess]
        return lengthen(species)


def shorten(species: str) -> str:
    try:
        return shorten.abbr[species]  # pyright: ignore[reportFunctionMemberAccess]
    except AttributeError:
        shorten.abbr = species_abbr()  # pyright: ignore[reportFunctionMemberAccess]
        return shorten(species)


def species_abbr() -> dict[str, str]:
    toml = db_prefix() / "Species_abbr.toml"
    _log.info(toml)
    if toml.exists():
        with toml.open("rb") as fin:
            obj = tomllib.load(fin)
    else:
        with FTPplantregmap() as ftp, ftp.species_abbr_list().open() as fin:
            lines = fin.readlines()
        obj: dict[str, str] = {"Oryza_sativa_Japonica_Group": "Osj"}
        for line in lines[1:]:
            (short, long) = line.split("\t")
            obj[long.strip().replace(" ", "_")] = short
        with toml.open("wb") as fout:
            tomli_w.dump(obj, fout)
    return obj


class FTPplantregmap(dl.LazyFTP):
    def __init__(self) -> None:
        super().__init__(
            _FTP_HOST,
            "/pub/database/PlantRegMap",
            _prefix_mirror(),
            timeout=3600,
        )

    def ls_cache(self, species: str = "") -> None:
        self.nlst_cache("")
        self.nlst_cache("08-download")
        self.nlst_cache("08-download/FTP")
        self.nlst_cache("08-download/FTP/pairwise_alignments")
        if species:
            self.nlst_cache(f"08-download/{species}")

    def species_abbr_list(self) -> Path:
        return self.retrieve("Species_abbr.list")

    def download(self, species: str) -> list[cli.Future[Path]]:
        fts: list[cli.Future[Path]] = []
        self.ls_cache(species)
        for bedgraph in self.download_conservation(species):
            fasize = api.fasize(species)
            fts.append(cli.thread_submit(kent.bedGraphToBigWig, bedgraph, fasize))
        for _ in self.download_multiple_alignments(species):
            pass
        for _ in self.download_pairwise_alignments(species):
            pass
        return fts

    def download_pairwise_alignments(self, species: str) -> Iterator[Path]:
        sp = shorten(species)
        relpath = f"08-download/FTP/pairwise_alignments/{sp}"
        nlst = self.nlst_cache(relpath)
        yield from (self.retrieve_symlink(x, species) for x in nlst)

    def download_multiple_alignments(self, species: str) -> Iterator[Path]:
        relpath = f"08-download/{species}/multiple_alignments"
        nlst = self.nlst_cache(relpath)
        yield from (self.retrieve_symlink(x, species) for x in nlst)

    def download_conservation(self, species: str) -> Iterator[Path]:
        relpath = f"08-download/{species}/sequence_conservation"
        nlst = self.nlst_cache(relpath)
        yield from (self.retrieve_symlink(x, species) for x in nlst)

    def retrieve_symlink(self, relpath: str, species: str) -> Path:
        orig = self.retrieve(relpath, checksize=True)
        outdir = db_prefix() / species / "compara"
        if re.match(r"Osj_Hvu\.(chain|net)\.gz", orig.name):
            return filter_chr_chainnet(orig, outdir / orig.name, "_unordered")
        if re.match(r"Osj_(Obr|Obl)\.(chain|net)\.gz", orig.name):
            return filter_chr_chainnet(orig, outdir / orig.name)
        return fs.symlink(orig, outdir / orig.name, relative=True)


def filter_chr_chainnet(infile: Path, outfile: Path, remove: str = "") -> Path:
    if not fs.is_outdated(outfile, infile):
        return outfile
    filename = infile.name.removesuffix(".gz")
    stem = Path(filename).stem
    _target, query = stem.split("_")
    qspecies = lengthen(query)
    qseqs = api.chrom_sizes(qspecies).keys()
    popen_filter = kent.netFilter if filename.endswith(".net") else kent.chainFilter
    with (
        subp.popen_zcat(infile) as zcat,
        subp.popen_sd(remove, stdin=zcat.stdout) as sd,
        popen_filter(stdin=sd.stdout, q=",".join(qseqs)) as pfilter,
    ):
        return subp.gzip(pfilter.stdout, outfile)


if __name__ == "__main__":
    main()
