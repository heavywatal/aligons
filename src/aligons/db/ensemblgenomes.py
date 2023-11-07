"""https://plants.ensembl.org/.

- {db.root}/ftp.ensemblgenomes.org/pub/plants/release-{version}/
  - fasta/{species}/dna/{stem}.{version}.dna_sm.chromosome.{seqid}.fa.gz
  - gff3/{species}/{stem}.{version}.chromosome.{seqid}.gff3.gz
  - maf/ensembl-compara/pairwise_alignments/
"""
import logging
import os
import re
from collections.abc import Iterable
from pathlib import Path

from aligons.util import cli, config, dl, fs, subp

from . import _rsrc, api, phylo, tools

_log = logging.getLogger(__name__)
_HOST = "ftp.ensemblgenomes.org"


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-V", "--versions", action="store_true")
    parser.add_argument("-a", "--all", action="store_true")
    parser.add_argument("--fmt", default="fasta", choices=("fasta", "gff3", "maf"))
    args = parser.parse_args(argv or None)
    if args.versions:
        for x in sorted(_list_versions()):
            print(x)
    elif args.all:
        with FTPensemblgenomes() as ftp:
            for sp in ftp.available_species():
                print(sp)
        _log.info(f"{version()=}")
    elif (fmt_dir := _prefix_mirror() / args.fmt).exists():
        for x in fmt_dir.iterdir():
            if x.is_dir():
                print(x)
    else:
        _log.warning(f"No local mirror of release-{version()}")


def _list_versions() -> Iterable[Path]:
    _log.debug(f"{_prefix_mirror()=}")
    return _prefix_mirror().parent.glob("release-*")


def download_compara(species: str):
    assert species in phylo.list_species(), species
    with FTPensemblgenomes() as ftp:
        dirs = ftp.download_maf(species)
    pool = cli.ThreadPool()
    fts = [pool.submit(consolidate_compara_mafs, d) for d in dirs]
    cli.wait_raise(fts)


def consolidate_compara_mafs(indir: Path) -> Path:
    _log.debug(f"{indir=}")
    mobj = re.search(r"([^_]+)_.+?\.v\.([^_]+)", indir.name)
    assert mobj, indir.name
    target_short = mobj.group(1)
    query_short = mobj.group(2)
    target = next(phylo.expand_shortnames([target_short]))
    query = next(phylo.expand_shortnames([query_short]))
    outdir = Path("compara") / target / query
    pat = re.compile(r"lastz_net\.([^_]+)_\d+\.maf$")
    infiles_by_seq: dict[str, list[Path]] = {}
    for maf in fs.sorted_naturally(indir.glob("*_*.maf")):
        mobj = pat.search(maf.name)
        assert mobj, maf.name
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
        with (
            sing_maf.open("wb") as fout,
            subp.popen("mafFilter stdin", stdin=subp.PIPE, stdout=fout) as maff,
        ):
            sed = f"sed -e 's/{target}/{target_short}/' -e 's/{query}/{query_short}/'"
            subp.run(sed, input="".join(lines).encode(), stdout=maff.stdin)
            maff.communicate()
            # for padding, not for filtering
    _log.info(f"{outdir}")
    return outdir


def readlines_compara_maf(file: Path) -> Iterable[str]:
    """MAF files of ensembl compara have broken "a" lines.

    a# id: 0000000
     score=9999
    s aaa.1
    s bbb.1
    """
    with file.open("r") as fin:
        for line in fin:
            if line.startswith(("#", "a#")):
                continue
            if line.startswith(" score"):
                yield "a" + line
            else:
                yield line


def download_via_ftp(species: list[str]) -> None:
    with FTPensemblgenomes() as ftp:
        species = ftp.remove_unavailable(species)
        fasizes: list[cli.FuturePath] = []
        for sp in species:
            fasta_originals = ftp.download_chr_sm_fasta(sp)
            fasta_copies = [_ln_or_bgzip(x, sp) for x in fasta_originals]
            fasizes.append(tools.index_fasta(fasta_copies))
        futures: list[cli.FuturePath] = []
        for sp, future in zip(species, fasizes, strict=True):
            gff3_originals = ftp.download_chr_gff3(sp)
            gff3_copies = [_ln_or_bgzip(x, sp) for x in gff3_originals]
            fasize = future.result()
            futures.append(tools.index_gff3(gff3_copies, fasize))
        cli.wait_raise(futures)
    if config["db"]["kmer"]:
        cli.wait_raise([tools.softmask(sp) for sp in species])


def _ln_or_bgzip(src: Path, species: str) -> Path:
    fmt = "fasta" if src.name.removesuffix(".gz").endswith(".fa") else "gff3"
    dstname = src.name.replace("primary_assembly", "chromosome")
    dst = prefix() / fmt / species / dstname
    if ".chromosome." in dstname:
        fs.symlink(src, dst, relative=True)
    else:
        tools.bgzip_or_symlink(src, dst)
    return dst


class FTPensemblgenomes(dl.LazyFTP):
    def __init__(self) -> None:
        super().__init__(
            _HOST,
            f"/{_relpath_release()}",
            _prefix_mirror(),
        )

    def remove_unavailable(self, species: list[str]) -> list[str]:
        available = self.available_species()
        filtered: list[str] = []
        for sp in species:
            if sp in available:
                filtered.append(sp)
            else:
                _log.info(f"{sp} not found in ensemblegenomes {version()}")
        return filtered

    def available_species(self) -> list[str]:
        return [Path(x).name for x in self.nlst_cache("fasta")]

    def download_chr_sm_fasta(self, species: str) -> list[Path]:
        relpath = f"fasta/{species}/dna"
        outdir = self.prefix / relpath
        nlst = self.nlst_cache(relpath)
        files = [self.retrieve(x) for x in self.remove_duplicates(nlst, "_sm.")]
        fs.checksums(outdir / "CHECKSUMS")
        return [f for f in files if f.suffix == ".gz"]

    def download_chr_gff3(self, species: str) -> list[Path]:
        relpath = f"gff3/{species}"
        outdir = self.prefix / relpath
        nlst = self.nlst_cache(relpath)
        files = [self.retrieve(x) for x in self.remove_duplicates(nlst)]
        fs.checksums(outdir / "CHECKSUMS")
        return [f for f in files if f.suffix == ".gz"]

    def download_maf(self, species: str) -> list[Path]:
        relpath = "maf/ensembl-compara/pairwise_alignments"
        outdir = self.prefix / relpath
        nlst = self.nlst_cache(relpath)
        sp = phylo.shorten(species)
        for x in nlst:
            if f"/{sp}_" in x:
                self.retrieve(x)
        _log.debug(f"{outdir=}")
        dirs: list[Path] = []
        for targz in outdir.glob("*.tar.gz"):
            expanded = self.prefix / targz.with_suffix("").with_suffix("")
            tar = ["tar", "xzf", targz, "-C", outdir]
            subp.run(tar, if_=fs.is_outdated(expanded / "README.maf"))
            # TODO: MD5SUM
            dirs.append(expanded.resolve())
        return dirs

    def remove_duplicates(self, nlst: list[str], substr: str = "") -> list[str]:
        matched = [x for x in nlst if "chromosome" in x]
        if not matched:
            matched = [x for x in nlst if "primary_assembly" in x]
        if not matched:
            matched = [x for x in nlst if "toplevel" in x]
        if not matched:
            matched = [x for x in nlst if f"{version()}.gff3" in x]
        matched = [x for x in matched if "musa_acuminata_v2" not in x]  # v52
        if substr:
            matched = [x for x in matched if substr in x]
        assert matched, substr
        misc = [x for x in nlst if re.search("CHECKSUMS$|README$", x)]
        return matched + misc


def download_via_rsync(species: list[str]):
    for sp in species:
        options = "--include *_sm.chromosome.*.fa.gz --exclude *.gz"
        rsync(f"fasta/{sp}/dna", options)
        options = "--include *.chromosome.*.gff3.gz --exclude *.gz"
        rsync(f"gff3/{sp}", options)


def rsync(relpath: str, options: str = "") -> subp.CompletedProcess[bytes]:
    remote_prefix = f"rsync://{_HOST}/all/{_relpath_release()}"
    src = f"{remote_prefix}/{relpath}/"
    dst = _prefix_mirror() / relpath
    return subp.run(f"rsync -auv {options} {src} {dst}")


def prefix() -> Path:
    return api.prefix(f"ensembl-{version()}")


def _prefix_mirror() -> Path:
    return _rsrc.db_root(_HOST) / _relpath_release()


def _relpath_release() -> Path:
    return Path(f"pub/plants/release-{version()}")


def version() -> int:
    return int(os.getenv("ENSEMBLGENOMES_VERSION", config["ensemblgenomes"]["version"]))


if __name__ == "__main__":
    main()
