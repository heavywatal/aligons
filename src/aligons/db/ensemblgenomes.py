"""https://plants.ensembl.org/.

- {db.mirror}/ensemblgenomes.org/plants/release-{version}/
  - fasta/{species}/dna/{stem}.dna_sm.chromosome.{seqid}.fa.gz
  - gff3/{species}/
  - maf/ensembl-compara/pairwise_alignments/
"""
import logging
import os
import re
from collections.abc import Iterable
from pathlib import Path

from aligons import db
from aligons.util import cli, config, dl, fs, subp

from . import phylo

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
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
    _log.debug(f"{_prefix_mirror_root()=}")
    return _prefix_mirror_root().glob("release-*")


def consolidate_compara_mafs(indir: Path):
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
        with sing_maf.open("wb") as fout:
            cmd = f"sed -e 's/{target}/{target_short}/' -e 's/{query}/{query_short}/'"
            sed = subp.popen(cmd, stdin=subp.PIPE, stdout=subp.PIPE)
            maff = subp.popen("mafFilter stdin", stdin=sed.stdout, stdout=fout)
            # for padding, not for filtering
            assert sed.stdout, cmd
            sed.stdout.close()
            sed.communicate("".join(lines).encode())
            maff.communicate()
    _log.info(f"{outdir}")
    return outdir


def readlines_compara_maf(file: Path):
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


class FTPensemblgenomes(dl.LazyFTP):
    def __init__(self):
        super().__init__(
            "ftp.ensemblgenomes.org",
            f"/pub/plants/release-{version()}",
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

    def download_fasta(self, species: str):
        relpath = f"fasta/{species}/dna"
        outdir = self.prefix / relpath
        nlst = self.nlst_cache(relpath)
        files = [self.retrieve(x) for x in self.remove_duplicates(nlst, "_sm.")]
        fs.checksums(outdir / "CHECKSUMS")
        return [f for f in files if f.suffix == ".gz"]

    def download_gff3(self, species: str):
        relpath = f"gff3/{species}"
        outdir = self.prefix / relpath
        nlst = self.nlst_cache(relpath)
        files = [self.retrieve(x) for x in self.remove_duplicates(nlst)]
        fs.checksums(outdir / "CHECKSUMS")
        return [f for f in files if f.suffix == ".gz"]

    def download_maf(self, species: str):
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

    def remove_duplicates(self, nlst: list[str], substr: str = ""):
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


def rsync(relpath: str, options: str = ""):
    server = "ftp.ensemblgenomes.org"
    remote_prefix = f"rsync://{server}/all/pub/plants/release-{version()}"
    src = f"{remote_prefix}/{relpath}/"
    dst = _prefix_mirror() / relpath
    return subp.run(f"rsync -auv {options} {src} {dst}")


def prefix():
    return db.path(f"ensembl-{version()}")


def _prefix_mirror():
    return _prefix_mirror_root() / f"release-{version()}"


def _prefix_mirror_root():
    return db.path_mirror("ensemblgenomes.org/plants")


def version():
    return int(os.getenv("ENSEMBLGENOMES_VERSION", config["ensemblgenomes"]["version"]))


if __name__ == "__main__":
    main()
