"""https://plants.ensembl.org/
"""
import functools
import logging
import os
import re
from collections.abc import Callable, Iterable
from ftplib import FTP
from pathlib import Path

from ..util import cli, config, fs, subp

_log = logging.getLogger(__name__)
VERSION = int(os.getenv("ENSEMBLGENOMES_VERSION", config["ensemblgenomes"]["version"]))
assert VERSION > 0
LOCAL_DB_ROOT = config["db"]["root"] / "ensemblgenomes/plants"
PREFIX = LOCAL_DB_ROOT / f"release-{VERSION}"


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-V", "--versions", action="store_true")
    parser.add_argument("-a", "--all", action="store_true")
    parser.add_argument("-f", "--files", action="store_true")
    parser.add_argument("-g", "--glob", default="*")
    parser.add_argument("--name", action="store_true")
    parser.add_argument("species", nargs="*")
    args = parser.parse_args(argv or None)
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
    for x in glob(args.glob, species):
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
    if VERSION > 50:
        monocot = f"({poaceae},dioscorea_rotundata)monocot"

    _solanum = "(solanum_lycopersicum,solanum_tuberosum)"
    solanaceae = f"(({_solanum},capsicum_annuum),nicotiana_attenuata)solanaceae"
    _convolvulaceae = "ipomoea_triloba"
    solanales = f"({solanaceae},{_convolvulaceae})solanales"
    _lamiales = "olea_europaea_sylvestris"
    if VERSION > 51:
        _lamiales = f"({_lamiales},sesamum_indicum)"
    lamiids = f"(({solanales},coffea_canephora),{_lamiales})lamiids"
    _asteraceae = "helianthus_annuus"
    if VERSION > 52:
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
    with FTPensemblgenomes() as ftp:
        lst = ftp.nlst_cache("fasta")
    return [Path(x).name for x in lst]


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


def get_file(pattern: str, species: str, subdir: str = ""):
    found = list(glob(pattern, [species], subdir))
    _log.debug(f"{found=}")
    assert len(found) == 1
    return found[0]


def glob(pattern: str, species: list[str] = [], subdir: str = ""):
    for path in species_dirs("fasta", species):
        for x in fs.sorted_naturally((path / "dna" / subdir).glob(pattern)):
            yield x
    for path in species_dirs("gff3", species):
        for x in fs.sorted_naturally(path.glob(pattern)):
            yield x


def expand_shortnames(shortnames: list[str]):
    return filter_by_shortname(species_names_all(), shortnames)


def filter_by_shortname(species: Iterable[str], queries: Iterable[str]):
    return (x for x in species if shorten(x) in queries)


def shorten(name: str):
    """Oryza_sativa -> osat"""
    if name.lower() == "olea_europaea_sylvestris":
        return "oesy"
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


def consolidate_compara_mafs(indir: Path):
    _log.debug(f"{indir=}")
    mobj = re.search(r"([^_]+)_.+?\.v\.([^_]+)", indir.name)
    assert mobj
    target_short = mobj.group(1)
    query_short = mobj.group(2)
    target = list(expand_shortnames([target_short]))[0]
    query = list(expand_shortnames([query_short]))[0]
    outdir = Path("compara") / target / query
    pat = re.compile(r"lastz_net\.([^_]+)_\d+\.maf$")
    infiles_by_seq: dict[str, list[Path]] = {}
    for maf in fs.sorted_naturally(indir.glob("*_*.maf")):
        mobj = pat.search(maf.name)
        assert mobj
        seq = mobj.group(1)
        infiles_by_seq.setdefault(seq, []).append(maf)
    for seq, infiles in infiles_by_seq.items():
        if seq == "supercontig":
            continue
        chrdir = outdir / f"chromosome.{seq}"
        chrdir.mkdir(0o755, parents=True, exist_ok=True)
        sing_maf = chrdir / "sing.maf"
        _log.info(str(sing_maf))
        if not fs.is_outdated(sing_maf, infiles):
            continue
        lines: list[str] = ["##maf version=1 scoring=LASTZ_NET\n"]
        for maf in infiles:
            lines.extend(readlines_compara_maf(maf))
        with open(sing_maf, "wb") as fout:
            cmd = f"sed -e 's/{target}/{target_short}/' -e 's/{query}/{query_short}/'"
            sed = subp.popen(cmd, stdin=subp.PIPE, stdout=subp.PIPE)
            maff = subp.popen("mafFilter stdin", stdin=sed.stdout, stdout=fout)
            # for padding, not for filtering
            assert sed.stdout
            sed.stdout.close()
            sed.communicate("".join(lines).encode())
            maff.communicate()
    return outdir


def readlines_compara_maf(file: Path):
    """maf files of ensembl compara have broken "a" lines

    a# id: 0000000
     score=9999
    s aaa.1
    s bbb.1
    """
    lines: list[str] = []
    with open(file, "r") as fin:
        for line in fin:
            if line.startswith("#") or line.startswith("a#"):
                continue
            elif line.startswith(" score"):
                line = "a" + line
            lines.append(line)
    return lines


class FTPensemblgenomes(FTP):
    def __init__(self):
        _log.info("FTP()")
        super().__init__()

    def quit(self):
        _log.info(f"os.chdir({self.orig_wd})")
        os.chdir(self.orig_wd)
        _log.info("ftp.quit()")
        super().quit()

    def lazy_init(self):
        if self.sock is not None:
            return
        host = "ftp.ensemblgenomes.org"
        _log.debug(f"ftp.connect({host})")
        _log.info(self.connect(host))
        _log.debug("ftp.login()")
        _log.info(self.login())
        path = f"/pub/plants/release-{VERSION}"
        _log.info(f"ftp.cwd({path})")
        _log.info(self.cwd(path))
        _log.info(f"os.chdir({PREFIX})")
        self.orig_wd = os.getcwd()
        os.chdir(PREFIX)  # for RETR only

    def download_fasta(self, species: str):
        pattern = r"/CHECKSUMS|/README"
        pattern += r"|_sm\.chromosome\..+\.fa\.gz$"
        pattern += r"|_sm\.primary_assembly\..+\.fa\.gz$"
        dir = self.download(f"fasta/{species}/dna", pattern)
        fs.checksums(dir / "CHECKSUMS")
        return dir

    def download_gff3(self, species: str):
        pattern = r"/CHECKSUMS|/README"
        pattern += r"|\.chromosome\..+\.gff3\.gz$"
        pattern += r"|\.primary_assembly\..+\.gff3\.gz$"
        dir = self.download(f"gff3/{species}", pattern)
        fs.checksums(dir / "CHECKSUMS")
        return dir

    def download_maf(self, species: str):
        dir = "maf/ensembl-compara/pairwise_alignments"
        sp = shorten(species)
        path = self.download(dir, f"/{sp}_")
        _log.debug(f"{path=}")
        dirs: list[Path] = []
        for targz in path.glob("*.tar.gz"):
            expanded = PREFIX / targz.with_suffix("").with_suffix("")
            tar = ["tar", "xzf", targz, "-C", path]
            subp.run_if(fs.is_outdated(expanded / "README.maf"), tar)
            # TODO: MD5SUM
            dirs.append(expanded.resolve())
        return dirs

    def download(self, dir: str, pattern: str):
        for x in self.nlst_search(dir, pattern):
            print(self.retrieve(x))
        return PREFIX / dir

    def nlst_search(self, dir: str, pattern: str):
        rex = re.compile(pattern)
        for x in self.nlst_cache(dir):
            if rex.search(x):
                yield x

    def nlst_cache(self, dir: str):
        cache = PREFIX / dir / ".ftp_nlst_cache"
        if cache.exists():
            _log.info(f"{cache=}")
            with cache.open("r") as fin:
                names = fin.read().rstrip().splitlines()
            lst = [str(Path(dir) / x) for x in names]
        else:
            lst = self.nlst(dir)
            cache.parent.mkdir(0o755, parents=True, exist_ok=True)
            with cache.open("w") as fout:
                fout.write("\n".join([Path(x).name for x in lst]) + "\n")
        return lst

    def retrieve(self, path: str):
        outfile = PREFIX / path
        if not outfile.exists() and not cli.dry_run:
            outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
            with open(outfile, "wb") as fout:
                cmd = f"RETR {path}"
                _log.info(self.retrbinary(cmd, fout.write))
        common = Path(str(outfile).replace("primary_assembly", "chromosome"))
        if not common.exists():
            common.symlink_to(outfile)
        return common

    # pyright: reportIncompatibleMethodOverride=false
    def nlst(self, dir: str):  # ensembl does not support mlsd
        self.lazy_init()
        _log.info(f"ftp.nlst({dir})")
        return super().nlst(dir)

    # pyright: reportIncompatibleMethodOverride=false
    def retrbinary(self, cmd: str, callback: Callable[[bytes], int]):
        self.lazy_init()
        _log.info(f"ftp.retrbinary({cmd})")
        return super().retrbinary(cmd, callback)


def rsync(relpath: str, options: str = ""):
    server = "ftp.ensemblgenomes.org"
    remote_prefix = f"rsync://{server}/all/pub/plants/release-{VERSION}"
    src = f"{remote_prefix}/{relpath}/"
    dst = PREFIX / relpath
    return subp.run(f"rsync -auv {options} {src} {dst}")


if __name__ == "__main__":
    main()
