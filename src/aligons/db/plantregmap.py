"""http://plantregmap.gao-lab.org/."""
import logging
import re
from collections.abc import Iterator
from pathlib import Path

from aligons import db
from aligons.db import api, jgi
from aligons.extern import htslib, kent, mafs2cram
from aligons.util import cli, dl, fs, subp, tomli_w, tomllib

from . import tools

_log = logging.getLogger(__name__)
_HOST = "plantregmap.gao-lab.org"


from_jgi = {
    "Bdistachyon": "Brachypodium_distachyon",
}


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-G", "--genome", action="store_true")
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
    else:
        for x in fs.sorted_naturally(db_prefix().rglob(args.pattern)):
            print(x)


def fetch_and_bgzip():
    fts: list[cli.FuturePath] = []
    for entry in db.iter_builtin_dataset("plantregmap.toml"):
        fts.extend(tools.fetch_and_bgzip(entry, db_prefix()))
    for entry in iter_jgi_dataset():
        fts.extend(tools.fetch_and_bgzip(entry, db_prefix()))
    return fts


def iter_jgi_dataset():
    keys = from_jgi.keys()
    for entry in db.iter_dataset(jgi.dataset_toml()):
        if (jgi_sp := entry["species"]) in keys:
            prm_sp = from_jgi[jgi_sp]
            entry["species"] = prm_sp
            entry["label"] = shorten(prm_sp)
            yield entry


def retrieve_deploy(query: str):
    url = f"http://{_HOST}/download_ftp.php?{query}"
    relpath = query.split("/", 1)[1]
    rawfile = db.path_mirror(_HOST) / relpath
    outfile = db_prefix() / relpath
    if outfile.suffix in (".bed", ".gff"):
        outfile = outfile.with_suffix(outfile.suffix + ".gz")
    elif outfile.name.endswith(".gtf.gz"):
        outfile = outfile.with_suffix("").with_suffix(".gff.gz")
    content = dl.fetch(url, rawfile).content
    if outfile.suffix == ".gz":
        future = cli.thread_submit(tools.compress, content, outfile)
    else:
        future = cli.thread_submit(fs.symlink, rawfile, outfile, relative=True)
    return cli.thread_submit(htslib.try_index, future)


def iter_download_queries():
    for query in iter_download_queries_all():
        if re.search(r"Oryza_sativa_Japonica|Solanum_lycopersicum", query):
            yield query


def iter_download_queries_all():
    content = download_php()
    for mobj in re.finditer(r"download_ftp\.php\?([^\"']+)", content):
        yield mobj[1]


def download_php() -> str:
    url = f"http://{_HOST}/download.php"
    cache = db.path_mirror(_HOST) / "download.php.html"
    return dl.fetch(url, cache).text


def db_prefix():
    return db.path("plantregmap")


def rglob(pattern: str, species: str = ".") -> Iterator[Path]:
    for species_dir in db_prefix().iterdir():
        if re.search(species, species_dir.name, re.IGNORECASE):
            yield from species_dir.rglob(pattern)


def to_cram(link: Path, species: str) -> Path:
    outfile = link.with_suffix(".cram")
    maf = gunzip(link)
    reference = api.genome_fa(species)
    return mafs2cram.maf2cram(maf, outfile, reference)


def to_bigwig(link: Path, species: str) -> Path:
    bedgraph = gunzip(link)
    bigwig = bedgraph.with_suffix(".bw")
    if fs.is_outdated(bigwig, bedgraph):
        kent.bedGraphToBigWig(bedgraph, api.fasize(species))
    return bigwig


def gunzip(infile: Path):
    outfile = infile.with_name(infile.name.removesuffix(".gz"))
    if fs.is_outdated(outfile, infile):
        subp.run(["gunzip", "-fk", infile])
        subp.run(["touch", outfile])
    return outfile


def download_via_ftp() -> list[cli.FuturePath]:
    fts: list[cli.FuturePath] = []
    with FTPplantregmap() as ftp:
        for entry in db.iter_builtin_dataset("plantregmap.toml"):
            species = entry["species"]
            fts.extend(ftp.download(species))
    return fts


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
    def __init__(self):
        host = "ftp.cbi.pku.edu.cn"
        super().__init__(
            host,
            "/pub/database/PlantRegMap",
            db.path_mirror(host) / "plantregmap",
            timeout=3600,
        )

    def ls_cache(self, species: str = ""):
        self.nlst_cache("")
        self.nlst_cache("08-download")
        self.nlst_cache("08-download/FTP")
        self.nlst_cache("08-download/FTP/pairwise_alignments")
        if species:
            self.nlst_cache(f"08-download/{species}")

    def species_abbr_list(self):
        return self.retrieve("Species_abbr.list")

    def download(self, species: str) -> list[cli.FuturePath]:
        fts: list[cli.FuturePath] = []
        self.ls_cache(species)
        for bedgraph in self.download_conservation(species):
            fts.append(cli.thread_submit(to_bigwig, bedgraph, species))  # noqa: PERF401
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
        return fs.symlink(orig, outdir / orig.name, relative=True)


if __name__ == "__main__":
    main()
