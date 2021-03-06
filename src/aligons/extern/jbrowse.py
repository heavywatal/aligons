"""
src: {ensemblgenomes.prefix}/fasta/{species}/dna/*.fa.gz
src: {ensemblgenomes.prefix}/gff3/{species}/*.chr.gff3.gz
src: {vNN}/pairwise/{species}/{query}/cram/genome.cram
src: {vNN}/multiple/{species}/{clade}/phastcons.bw
dst: {vNN}/{jbrowse_XYZ}/{species}/config.json
dst: {document_root}/{jbrowse_XYZ}/{vNN}/{species} ->
"""
import importlib.resources as resources
import json
import logging
import os
import re
from pathlib import Path
from typing import Any, TypeAlias

from ..db import ensemblgenomes, phylo, plantregmap, stat
from ..util import cli, config, fs, subp

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
        clades = [x.name for x in self.multiple_dir.iterdir()]
        clades = phylo.sorted_by_len_newicks(clades, reverse=True)
        for clade in clades:
            wig = self.multiple_dir / clade / "phastcons.bw"
            self.add_track(wig, "conservation", id=clade, subdir=clade)
        for bed in self.multiple_dir.rglob("cns.bed.gz"):
            clade = bed.parent.name
            self.add_track(bed, "conservation", id="CNS-" + clade, subdir=clade)
        gen = self.pairwise_dir.rglob("genome.cram")
        crams = {cram.parent.parent.name: cram for cram in gen}
        for query in phylo.extract_tip_names(phylo.newicks[clades[0]]):
            if cram := crams.pop(query, None):
                self.add_track(cram, "alignment", id=query, subdir=query)
        for query, cram in crams.items():
            self.add_track(cram, "alignment", id=query, subdir=query)
        if self.target.name == "oryza_sativa":
            f = config["db"]["root"] / "suzuemon/SV.bed.gz"
            if f.exists():
                self.add_track(f, category="misc", id="suzuemon-SV", subdir="suzuemon")
        for path in fs.sorted_naturally(plantregmap.rglob("*.gff.gz", species)):
            id = re.sub(r"_[^_]+\.gff\.gz$", "", path.name)
            self.add_track(path, "plantregmap", id=id)
        for path in fs.sorted_naturally(plantregmap.rglob("*.bed.gz", species)):
            id = re.sub(r"(_normal)?\.bed\.gz$", "", path.name)
            self.add_track(path, "plantregmap", id=id)
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
        if (csi := Path(str(file) + ".csi")).exists():
            args.extend(["--indexFile", csi])
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
        rex = re.compile(r"_inProm|_CE_genome-wide")  # redundant subsets
        tracks = [x for x in self.tracks if not rex.search(x)]
        args.extend(["--tracks", ",".join(tracks)])
        jbrowse(args)

    def configure(self):
        chrom_sizes = {x[0]: x[1] for x in stat.chrom_sizes(self.target.name)}
        config_json = self.target / "config.json"
        with open(config_json) as fin:
            cfg = json.load(fin)
        assembly = cfg["assemblies"][0]
        session = cfg["defaultSession"]
        view = session["views"][0]
        chrom = "6"
        start = 27475500
        view["bpPerPx"] = 5.0
        view["offsetPx"] = int(start / view["bpPerPx"])
        view["displayedRegions"] = [
            {
                "refName": chrom,
                "start": 0,
                "end": int(chrom_sizes[chrom]),
                "reversed": False,
                "assemblyName": assembly["name"],
            }
        ]
        for track in view["tracks"]:
            track["displays"] = [make_display(track)]
        if refnamealiases := self.make_refnamealiases():
            assembly["refNameAliases"] = refnamealiases
        cfg["configuration"] = make_configuration()
        with open(config_json, "w") as fout:
            json.dump(cfg, fout, indent=2)

    def make_refnamealiases(self):
        species = self.target.name
        filename = f"{species}.chromAlias.txt"
        if not resources.is_resource("aligons.data", filename):
            return None
        with open(self.target / filename, "w") as fout:
            fout.write(resources.read_text("aligons.data", filename))
        return {
            "adapter": {
                "type": "RefNameAliasAdapter",
                "location": {
                    "uri": filename,
                    "locationType": "UriLocation",
                },
            }
        }


def make_display(track: dict[str, Any]):
    clade_color = {
        "bep": "#C82828",
        "poaceae": "#C8641E",
        "monocot": "#C8B414",
    }
    item = {}
    if track["type"] == "FeatureTrack":
        if "gff3" in track["type"]:
            item = {
                "type": "LinearBasicDisplay",
                "height": 80,
            }
        else:  # bed
            item = {
                "type": "LinearBasicDisplay",
                "height": 30,
                "trackShowLabels": False,
                "renderer": {
                    "type": "SvgFeatureRenderer",
                    "height": 10,
                    "color": "#800000",
                },
            }
    elif track["type"] == "QuantitativeTrack":
        item = {
            "type": "LinearWiggleDisplay",
            "height": 40,
            "color": clade_color.get(track["configuration"], "#888888"),
            "constraints": {"max": 1, "min": 0},
        }
    elif track["type"] == "AlignmentsTrack":
        item = {
            "type": "LinearPileupDisplay",
            "height": 20,
        }
    item["configuration"] = "-".join([track["configuration"], item["type"]])
    return item


def make_configuration():
    return {"theme": make_theme()}


def make_theme():
    return {
        "palette": {
            "primary": {"main": "#186038"},
            "secondary": {"main": "#009259"},
            "tertiary": {"main": "#8fc21f"},
            "quaternary": {"main": "#d9e000"},
        }
    }


def jbrowse(args: subp.Args, cond: bool = True):
    cmd: subp.Args = ["jbrowse"]
    cmd.extend(args)
    subp.run_if(cond, cmd)


def iter_targets(path: Path):
    for config_json in path.rglob("config.json"):
        if "test_data" in str(config_json):
            continue
        config_json = config_json.resolve()
        if config_json.parent.parent.name.startswith("jbrowse-"):
            yield config_json.parent.resolve()


def redirect_html(url: str):
    meta = f"""<meta http-equiv="refresh" content="1; URL={url}">"""
    return f"""<html><head>{meta}</head><body>Redirecting</body></html>"""


if __name__ == "__main__":
    main()
