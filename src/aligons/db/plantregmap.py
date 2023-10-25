"""http://plantregmap.gao-lab.org/."""
import logging
import re
from collections.abc import Iterator
from pathlib import Path

from aligons import db
from aligons.db import api, jgi
from aligons.extern import htslib, kent, mafs2cram
from aligons.util import cli, dl, fs, tomli_w, tomllib

from . import tools

_log = logging.getLogger(__name__)
_HOST = "plantregmap.gao-lab.org"


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-G", "--genome", action="store_true")
    parser.add_argument("-C", "--to-cram", type=Path)
    parser.add_argument("pattern", nargs="?", default="*")
    args = parser.parse_args(argv or None)
    if args.download:
        fts: list[cli.FuturePath] = []
        fts.extend(download_via_ftp())
        fts.extend(retrieve_deploy(q) for q in iter_download_queries())
        cli.wait_raise(fts)
    elif args.genome:
        fts = tools.index_as_completed(fetch_and_bgzip())
        cli.wait_raise(fts)
    elif args.to_cram:
        maf = net_to_maf(args.to_cram)
        to_cram(maf)
    else:
        for x in fs.sorted_naturally(db_prefix().rglob(args.pattern)):
            print(x)


def fetch_and_bgzip() -> list[cli.FuturePath]:
    fts: list[cli.FuturePath] = []
    for entry in db.iter_builtin_dataset("plantregmap.toml"):
        fts.extend(tools.fetch_and_bgzip(entry, db_prefix()))
    for entry in iter_jgi_dataset():
        fts.extend(tools.fetch_and_bgzip(entry, db_prefix()))
    return fts


def iter_jgi_dataset() -> Iterator[db.DataSet]:
    long_species = [
        "Brachypodium_distachyon",
        "Solanum_lycopersicum",
    ]
    from_jgi = {jgi.shorten(x): x for x in long_species}
    for entry in jgi.iter_dataset(from_jgi.keys()):
        prm_species = from_jgi[entry["species"]]
        entry["species"] = prm_species
        entry["label"] = shorten(prm_species)
        yield entry


def retrieve_deploy(query: str) -> cli.FuturePath:
    url = f"http://{_HOST}/download_ftp.php?{query}"
    relpath = query.split("/", 1)[1]
    rawfile = db.path_mirror(_HOST) / relpath
    outfile = db_prefix() / relpath
    if outfile.suffix in (".bed", ".gff"):
        outfile = outfile.with_suffix(outfile.suffix + ".gz")
    elif outfile.name.endswith(".gtf.gz"):
        outfile = outfile.with_suffix("").with_suffix(".gff.gz")
    response = dl.fetch(url, rawfile)
    if outfile.suffix == ".gz":
        future = cli.thread_submit(tools.compress_lazy, response, outfile)
    else:
        future = cli.thread_submit(fs.symlink, response.path, outfile, relative=True)
    return cli.thread_submit(htslib.try_index, future)


def iter_download_queries() -> Iterator[str]:
    for query in iter_download_queries_all():
        if re.search(r"Oryza_sativa_Japonica|Solanum_lycopersicum", query):
            yield query


def iter_download_queries_all() -> Iterator[str]:
    content = download_php()
    for mobj in re.finditer(r"download_ftp\.php\?([^\"']+)", content):
        yield mobj[1]


def download_php() -> str:
    url = f"http://{_HOST}/download.php"
    cache = db.path_mirror(_HOST) / "download.php.html"
    return dl.fetch(url, cache).text


def db_prefix() -> Path:
    return db.path("plantregmap")


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
    outfile = maf.with_suffix(".cram")
    reference = api.genome_fa(target)
    mafs2cram.maf2cram(maf, outfile, reference)
    htslib.index(outfile)
    return outfile


def download_via_ftp() -> list[cli.FuturePath]:
    targets = [
        "Oryza_sativa_Japonica_Group",
        "Solanum_lycopersicum",
    ]
    fts: list[cli.FuturePath] = []
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
    _log.info(f"{toml}")
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
        host = "ftp.cbi.pku.edu.cn"
        super().__init__(
            host,
            "/pub/database/PlantRegMap",
            db.path_mirror(host) / "plantregmap",
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

    def download(self, species: str) -> list[cli.FuturePath]:
        fts: list[cli.FuturePath] = []
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
            return sanitize_hvu_chainnet(orig, outdir / orig.name)
        return fs.symlink(orig, outdir / orig.name, relative=True)


def sanitize_hvu_chainnet(infile: Path, outfile: Path) -> Path:
    if not fs.is_outdated(outfile, infile):
        return outfile
    notq = [
        "morex_contig_137999",
        "morex_contig_1587688",
        "morex_contig_159532",
        "morex_contig_164803",
        "morex_contig_244380",
        "morex_contig_2549682",
        "morex_contig_41360",
        "morex_contig_41603",
        "morex_contig_42365",
        "morex_contig_43017",
        "morex_contig_43925",
        "morex_contig_44379",
        "morex_contig_61860",
        "morex_contig_68624",
        "morex_contig_68638",
        "morex_contig_68864",
        "morex_contig_70023",
        "morex_contig_70567",
        "morex_contig_70976",
        "morex_contig_71199",
        "morex_contig_72677",
        "morex_contig_73216",
    ]
    content = kent.chain_net_filter(infile, notQ=",".join(notq))
    content = content.replace(b"_unordered", b"")
    content = fs.gzip_compress(content)
    with outfile.open("wb") as fout:
        fout.write(content)
    return outfile


if __name__ == "__main__":
    main()
