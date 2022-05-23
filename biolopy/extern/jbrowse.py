"""
src: {ensemblgenomes.prefix}/fasta/{species}/dna/*.fa.gz
src: {ensemblgenomes.prefix}/gff3/{species}/*.chr.gff3.gz
src: {vNN}/pairwise/{species}/{query}/cram/genome.cram
src: {vNN}/multiple/{species}/{clade}/phastcons.bw
dst: {vNN}/{jbrowse_XYZ}/{species}/config.json
dst: {document_root}/{jbrowse_XYZ}/{vNN}/{species} ->
"""
import json
import logging
import os
from pathlib import Path
from typing import Any, TypeAlias

from ..db import ensemblgenomes
from ..util import cli, subp

StrPath: TypeAlias = str | Path[str]

_log = logging.getLogger(__name__)
_load = "symlink"
env_version = os.environ["JBROWSE_VERSION"]


def main(argv: list[str] = []):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-a", "--admin", action="store_true")
    parser.add_argument("-u", "--upgrade", action="store_true")
    parser.add_argument("-d", "--deploy", action="store_true")
    parser.add_argument("-o", "--outdir", type=Path, default=Path("."))
    parser.add_argument("indir", type=Path)  # multiple/oryza_sativa/
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    jb = JBrowse(args.outdir)
    if args.admin:
        jb.admin_server()
        return
    if args.upgrade:
        jb.upgrade()
        return
    if args.deploy:
        for target in iter_targets(args.indir):
            jb.deploy(target)
        return
    jbc = JBrowseConfig(args.indir, args.outdir)
    jbc.add()
    jbc.configure()


class JBrowse:
    def __init__(self, document_root: Path, prefix: StrPath = ""):
        self.document_root = document_root
        self.prefix = Path(prefix)
        self.root = document_root / self.prefix / f"jbrowse-{env_version}"

    def deploy(self, target: Path):
        _log.debug(f"{target=}")
        jbrowse_XYZ = target.parent.name
        vNN = target.parent.parent.name
        relpath = Path(vNN, target.name)
        dst = self.document_root / self.prefix / jbrowse_XYZ / relpath
        self.create(jbrowse_XYZ.split("-")[1])
        _log.info(f"{dst} -> {target}")
        slug = self.prefix / jbrowse_XYZ
        url = f"/{slug}/?config={relpath}/config.json"
        if not cli.dry_run:
            if not dst.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.symlink_to(target)
            with open(target / "index.html", "w") as fout:
                fout.write(redirect_html(url))
        print(f"http://localhost/{slug / relpath}/ -> {url}")

    def create(self, version: str = env_version):
        root = self.document_root / self.prefix / f"jbrowse-{version}"
        args = ["create", root]
        args.append(f"--tag=v{version}")
        jbrowse(args, not root.exists())
        with open(root / "version.txt", "r") as fin:
            assert fin.read().strip() == version

    def admin_server(self):
        jbrowse(["admin-server", "--root", self.root])

    def upgrade(self):
        jbrowse(["upgrade", self.root])


class JBrowseConfig:
    def __init__(self, multialign_species: Path, outdir: Path):
        self.jb = JBrowse(outdir)
        self.multiple_dir = multialign_species
        species_name = self.multiple_dir.name
        vNN_dir = multialign_species.parent.parent
        self.pairwise_dir = vNN_dir / "pairwise" / species_name
        self.target = self.jb.root / species_name
        self.tracks: list[str] = []
        _log.info(f"{self.target}")

    def add(self, force: bool = False):
        self.jb.create()  # for preview and test
        self.target.mkdir(0o755, parents=True, exist_ok=True)
        species = self.multiple_dir.name
        self.add_assembly(species)
        self.add_track_gff(species)
        for wig in self.multiple_dir.rglob("phastcons.bw"):
            clade = wig.parent.name
            self.add_track(wig, "conservation", id=clade, subdir=clade)
        for cram in self.pairwise_dir.rglob("genome.cram"):
            query = cram.parent.parent.name
            self.add_track(cram, "alignment", id=query, subdir=query)
        self.set_default_session()

    def add_assembly(self, species: str):
        """--alias, --name, --displayName"""
        genome = ensemblgenomes.get_file("*.genome.fa.gz", species)
        args: subp.Args = ["add-assembly"]
        args.extend(["--target", self.target])
        args.extend(["--load", _load])
        args.append(genome)
        jbrowse(args, not (self.target / genome.name).exists())

    def add_track_gff(self, species: str):
        gff = ensemblgenomes.get_file("*.genome.gff3.gz", species)
        ngff = gff.with_suffix("").with_suffix("").with_suffix(".name.gff3.gz")
        if ngff.exists():
            gff = ngff
        self.add_track(gff, id=gff.parent.name + ".gff3")
        self.text_index()

    def add_track(self, file: Path, category: str = "", id: str = "", subdir: str = ""):
        """--description, --config"""
        args: subp.Args = ["add-track"]
        args.extend(["--target", self.target])
        args.extend(["--load", _load])
        if subdir:
            args.extend(["--subDir", subdir])
        if id:
            args.extend(["--trackId", id])
            self.tracks.append(id)
        if category:
            args.extend(["--category", category])
        args.append(file)
        jbrowse(args, not (self.target / subdir / file.name).exists())

    def text_index(self):
        """--attributes, --exclude, --file,  --perTrack, --tracks, --dryrun"""
        args: subp.Args = ["text-index"]
        args.extend(["--target", self.target])
        jbrowse(args)

    def set_default_session(self):
        args: subp.Args = ["set-default-session"]
        args.extend(["--target", self.target])
        args.extend(["--name", f"New {self.target.name} session"])
        args.extend(["--view", "LinearGenomeView"])
        args.extend(["--tracks", ",".join(self.tracks)])
        jbrowse(args)

    def configure(self):
        config_json = self.target / "config.json"
        with open(config_json) as fin:
            cfg = json.load(fin)
        session = cfg["defaultSession"]
        view = session["views"][0]
        view["offsetPx"] = 5002000
        view["bpPerPx"] = 2.0
        view["displayedRegions"] = [
            {
                "refName": "1",
                "start": 0,
                "end": 43270923,
                "reversed": False,
                "assemblyName": "Oryza_sativa.IRGSP-1.0.dna_sm.genome",
            }
        ]
        for track in view["tracks"]:
            track["displays"] = [make_display(track)]
        with open(config_json, "w") as fout:
            json.dump(cfg, fout, indent=2)


def make_display(track: dict[str, Any]):
    clade_color = {
        "bep": "#C82828",
        "poaceae": "#C8641E",
        "monocot": "#C8B414",
    }
    if track["type"] == "FeatureTrack":
        return {
            "type": "LinearBasicDisplay",
            "height": 80,
        }
    elif track["type"] == "QuantitativeTrack":
        return {
            "type": "LinearWiggleDisplay",
            "height": 40,
            "color": clade_color.get(track["configuration"], "#888888"),
            "constraints": {"max": 1, "min": 0},
        }
    elif track["type"] == "AlignmentsTrack":
        return {
            "type": "LinearPileupDisplay",
            "height": 20,
        }


def jbrowse(args: subp.Args, cond: bool = True):
    cmd: subp.Args = ["jbrowse"]
    cmd.extend(args)
    subp.run_if(cond, cmd)


def iter_targets(path: Path):
    for config in path.rglob("config.json"):
        if "test_data" in str(config):
            continue
        yield config.parent.resolve()


def redirect_html(url: str):
    meta = f"""<meta http-equiv="refresh" content="1; URL={url}">"""
    return f"""<html><head>{meta}</head><body>Redirecting</body></html>"""


if __name__ == "__main__":
    main()