"""https://plants.ensembl.org/
"""
import functools
import logging
import os
import re
from collections.abc import Iterable
from ftplib import FTP
from pathlib import Path

from ..util import cli, fs, subp

_log = logging.getLogger(__name__)
LOCAL_DB_ROOT = Path("~/db/ensemblgenomes/plants").expanduser()
VERSION = os.environ["ENSEMBLGENOMES_VERSION"]
PREFIX = LOCAL_DB_ROOT / f"release-{VERSION}"


def main(argv: list[str] | None = None):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-V", "--versions", action="store_true")
    parser.add_argument("-a", "--all", action="store_true")
    parser.add_argument("-f", "--files", action="store_true")
    parser.add_argument("-g", "--glob", default="*")
    parser.add_argument("--name", action="store_true")
    parser.add_argument("species", nargs="*")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.versions:
        for x in sorted(list_versions()):
            print(x)
        return
    if args.all and not args.files:
        species = species_names_all()
    else:
        species = species_names()
    if args.species:
        species = list(filter_by_shortname(species, args.species))
    if not args.files:
        for sp in species:
            print(sp)
        return
    for x in rglob(args.glob, species):
        if args.name:
            print(x.name)
        else:
            print(x)


def make_newicks():
    ehrhartoideae = "(oryza_sativa,leersia_perrieri)ehrhartoideae"
    pooideae = "(brachypodium_distachyon,(aegilops_tauschii,hordeum_vulgare))pooideae"
    andropogoneae = "(sorghum_bicolor,zea_mays)andropogoneae"
    paniceae = "(setaria_italica,panicum_hallii_fil2)paniceae"
    bep = f"({ehrhartoideae},{pooideae})bep"
    pacmad = f"({andropogoneae},{paniceae})pacmad"
    poaceae = f"({bep},{pacmad})poaceae"
    monocot = f"(({poaceae},musa_acuminata),dioscorea_rotundata)monocot"  # noqa: F841
    if int(VERSION) > 50:
        monocot = f"({poaceae},dioscorea_rotundata)monocot"

    _solanum = "(solanum_lycopersicum,solanum_tuberosum)"
    solanaceae = f"(({_solanum},capsicum_annuum),nicotiana_attenuata)solanaceae"
    _convolvulaceae = "ipomoea_triloba"
    solanales = f"({solanaceae},{_convolvulaceae})solanales"
    _lamiales = "olea_europaea_sylvestris"
    if int(VERSION) > 51:
        _lamiales = f"({_lamiales},sesamum_indicum)"
    lamiids = f"(({solanales},coffea_canephora),{_lamiales})lamiids"
    _asteraceae = "helianthus_annuus"
    if int(VERSION) > 52:
        _asteraceae = f"({_asteraceae},lactuca_sativa)"
    _companulids = f"({_asteraceae},daucus_carota)"
    _core_asterids = f"({lamiids},{_companulids})"
    asterids = f"({_core_asterids},actinidia_chinensis)asterids"
    eudicots = f"({asterids},arabidopsis_thaliana)eudicots"

    angiospermae = f"({eudicots},{monocot})angiospermae"  # noqa: F841
    # pyright: reportUnusedVariable=false
    return {k: v + ";" for k, v in locals().items() if not k.startswith("_")}


def list_versions():
    _log.debug(f"{LOCAL_DB_ROOT=}")
    return LOCAL_DB_ROOT.glob("release-*")


@functools.cache
def species_names_all():
    cache = PREFIX / "species.tsv"
    if not cache.exists():
        assert cache.parent.exists(), f"{cache.parent} exists"
        with FTPensemblgenomes() as ftp:
            species = [Path(x).name for x in ftp.nlst("fasta")]
        with open(cache, "wt") as fout:
            fout.write("\n".join(species) + "\n")
    _log.debug(f"{cache=}")
    with open(cache, "r") as fin:
        return [x.rstrip() for x in fin.readlines()]


@functools.cache
def species_names(format: str = "fasta"):
    return [x.name for x in species_dirs(format)]


def species_dirs(format: str = "fasta", species: list[str] = []):
    assert (root := PREFIX / format).exists(), root
    requests = set(species)
    for path in root.iterdir():
        if not path.is_dir():
            continue
        if not species or (path.name in requests):
            requests.discard(path.name)  # TODO: search twice
            yield path
    assert not requests, f"directory not found: {requests}"


def get_file(pattern: str, species: str):
    found = list(rglob(pattern, [species]))
    _log.debug(f"{found=}")
    assert len(found) == 1
    return found[0]


def rglob(pattern: str, species: list[str] = []):
    for path in species_dirs("fasta", species):
        for x in fs.sorted_naturally(path.rglob(pattern)):
            yield x
    for path in species_dirs("gff3", species):
        for x in fs.sorted_naturally(path.rglob(pattern)):
            yield x


def expand_shortnames(shortnames: list[str]):
    return filter_by_shortname(species_names(), shortnames)


def filter_by_shortname(species: Iterable[str], queries: Iterable[str]):
    return (x for x in species if shorten(x) in queries)


def shorten(name: str):
    """Oryza_sativa -> osat"""
    split = name.lower().split("_")
    return split[0][0] + split[1][:3]


def sanitize_queries(target: str, queries: list[str]):
    queries = list(dict.fromkeys(queries))
    try:
        queries.remove(target)
    except ValueError:
        pass
    assert queries
    _log.debug(f"{queries=}")
    assert set(queries) <= set(species_names())
    return queries


class FTPensemblgenomes(FTP):
    def __init__(self):
        host = "ftp.ensemblgenomes.org"
        _log.info(f"FTP({host})")
        super().__init__(host)
        _log.debug("ftp.login()")
        _log.info(self.login())
        _log.debug("ftp.getwelcome()")
        _log.info(self.getwelcome())
        path = f"/pub/plants/release-{VERSION}"
        _log.info(f"ftp.cwd({path})")
        _log.info(self.cwd(path))
        self.orig_wd = os.getcwd()
        _log.info(f"os.chdir({PREFIX})")
        os.chdir(PREFIX)

    def download_fasta(self, species: str):
        pattern = r"/CHECKSUMS|/README|_sm\.chromosome\..+fa\.gz$"
        pattern += r"|_sm\.primary_assembly\..+fa\.gz$"
        for x in self.nlst_search(f"fasta/{species}/dna", pattern):
            print(self.retrieve(x))

    def download_gff3(self, species: str):
        pattern = r"/CHECKSUMS|/README|\.chromosome\..+gff3\.gz$"
        pattern += r"|\.primary_assembly\..+gff3\.gz$"
        for x in self.nlst_search(f"gff3/{species}", pattern):
            print(self.retrieve(x))

    def nlst_search(self, dir: str, pattern: str):
        _log.info(f"ftp.nlst({dir})")  # ensembl does not support mlsd
        rex = re.compile(pattern)
        for x in self.nlst(dir):
            if rex.search(x):
                yield x

    def retrieve(self, path: str):
        outfile = Path(path)
        if not outfile.exists() and not cli.dry_run:
            outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
            with open(outfile, "wb") as fout:
                cmd = f"RETR {path}"
                _log.info(cmd)
                _log.info(self.retrbinary(cmd, fout.write))
        return outfile

    def quit(self):
        _log.info(f"os.chdir({self.orig_wd})")
        os.chdir(self.orig_wd)
        _log.debug("ftp.quit()")
        _log.info(ret := super().quit())
        return ret


def rsync(relpath: str, options: str = ""):
    server = "ftp.ensemblgenomes.org"
    remote_prefix = f"rsync://{server}/all/pub/plants/release-{VERSION}"
    src = f"{remote_prefix}/{relpath}/"
    dst = PREFIX / relpath
    return subp.run(f"rsync -auv {options} {src} {dst}")


if __name__ == "__main__":
    main()
