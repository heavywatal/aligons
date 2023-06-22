"""https://jbrowse.org.

src: {db.api.prefix}/fasta/{species}/*.fa.gz
src: {db.api.prefix}/gff3/{species}/*.chr.gff3.gz
src: {vNN}/pairwise/{species}/{query}/cram/genome.cram
src: {vNN}/multiple/{species}/{clade}/phastcons.bw
dst: {vNN}/{jbrowse_XYZ}/{species}/config.json
dst: {document_root}/{jbrowse_XYZ}/{vNN}/{species} ->
"""
import json
import logging
import re
from pathlib import Path
from typing import Any, TypeAlias

from aligons import db
from aligons.db import api, phylo, plantdhs, plantregmap, stat
from aligons.util import cli, config, fs, resources_data, subp

StrPath: TypeAlias = str | Path

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-a", "--admin", action="store_true")
    parser.add_argument("-u", "--upgrade", action="store_true")
    parser.add_argument("-d", "--deploy", action="store_true")
    parser.add_argument("-o", "--outdir", type=Path, default=Path("."))
    parser.add_argument("indir", type=Path)  # multiple/oryza_sativa/
    args = parser.parse_args(argv or None)
    if cli.dry_run:
        jbrowse(["version"])
        return
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
        self.load = config["jbrowse"]["load"]
        self.version = config["jbrowse"]["version"]
        self.document_root = document_root
        self.prefix = Path(prefix)
        self.root = document_root / self.prefix / f"jbrowse-{self.version}"

    def deploy(self, target: Path):
        _log.debug(f"{target=}")
        jbrowse_xyz = target.parent.name
        vnn = target.parent.parent.name
        relpath = Path(vnn, target.name)
        dst = self.document_root / self.prefix / jbrowse_xyz / relpath
        self.create(jbrowse_xyz.split("-")[1])
        _log.info(f"{dst} -> {target}")
        slug = self.prefix / jbrowse_xyz
        url = f"/{slug}/?config={relpath}/config.json"
        if not cli.dry_run:
            if not dst.exists():
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.symlink_to(target)
            with (target / "index.html").open("w") as fout:
                fout.write(redirect_html(url))
        print(f"http://localhost/{slug / relpath}/ -> {url}")

    def create(self, version: str = ""):
        version = version or self.version
        root = self.document_root / self.prefix / f"jbrowse-{version}"
        args = ["create", root]
        args.append(f"--tag=v{version}")
        if not root.exists():
            jbrowse(args)
        with (root / "version.txt").open() as fin:
            version_txt = fin.read().strip()
            assert version_txt == version, f"{version_txt=} != {version=}"

    def admin_server(self):
        jbrowse(["admin-server", "--root", self.root])

    def upgrade(self):
        jbrowse(["upgrade", self.root])


class JBrowseConfig:
    def __init__(self, multialign_species: Path, outdir: Path):
        self.jb = JBrowse(outdir)
        self.multiple_dir = multialign_species
        species_name = self.multiple_dir.name
        vnn_dir = multialign_species.parent.parent
        self.pairwise_dir = vnn_dir / "pairwise" / species_name
        self.target = self.jb.root / species_name
        self.tracks: list[str] = []
        _log.info(f"{self.target}")

    def add(self):
        self.jb.create()  # for preview and test
        self.target.mkdir(0o755, parents=True, exist_ok=True)
        species = self.multiple_dir.name
        self.add_assembly(species)
        self.add_track_gff(species)
        clades = [x.name for x in self.multiple_dir.iterdir() if "-" not in x.name]
        _log.info(f"{clades}")
        clades = phylo.sorted_by_len_newicks(clades, reverse=True)
        for clade in clades:
            wig = self.multiple_dir / clade / "phastcons.bw"
            self.add_track(wig, "conservation", trackid=clade, subdir=clade)
        for bed in self.multiple_dir.rglob("cns.bed.gz"):
            clade = bed.parent.name
            self.add_track(bed, "conservation", trackid="CNS-" + clade, subdir=clade)
        gen = self.pairwise_dir.rglob("genome.cram")
        crams = {cram.parent.parent.name: cram for cram in gen}
        for query in phylo.list_species(clades[0]):
            if cram := crams.pop(query, None):
                self.add_track(cram, "alignment", trackid=query, subdir=query)
        for query, cram in crams.items():
            self.add_track(cram, "alignment", trackid=query, subdir=query)
        if self.target.name == "oryza_sativa":
            self.add_papers_data()
            self.add_plantdhs()
        self.add_plantregmap(species)
        self.set_default_session()

    def add_papers_data(self):
        for path in fs.sorted_naturally(db.path("papers").glob("*.bed.gz")):
            self.add_track(path, "papers", trackid=path.with_suffix("").stem)
        suzuemon = db.path("suzuemon")
        if (f := suzuemon / "sv_with_DEG.bed.gz").exists():
            self.add_track(f, "papers", trackid="SV_DEG-qin2021", subdir="suzuemon")
        if (f := suzuemon / "SV.bed.gz").exists():
            self.add_track(f, "papers", trackid="SV_all-qin2021", subdir="suzuemon")

    def add_plantdhs(self):
        for path in fs.sorted_naturally(plantdhs.db_prefix().glob("Rice_*.bw")):
            trackid = path.stem.removeprefix("Rice_")
            self.add_track(path, "plantdhs", trackid=trackid)
        for path in fs.sorted_naturally(plantdhs.db_prefix().glob("*.gff.gz")):
            self.add_track(path, "plantdhs", trackid=path.stem)

    def add_plantregmap(self, species: str):
        for path in fs.sorted_naturally(plantregmap.rglob("*.gff.gz", species)):
            trackid = re.sub(r"_[^_]+\.gff\.gz$", "", path.name)
            self.add_track(path, "plantregmap", trackid=trackid)
        for path in fs.sorted_naturally(plantregmap.rglob("*.bed.gz", species)):
            trackid = re.sub(r"(_normal)?\.bed\.gz$", "", path.name)
            self.add_track(path, "plantregmap", trackid=trackid)

    def add_assembly(self, species: str):
        # --alias, --name, --displayName
        genome = api.genome_fa(species)
        args: subp.Args = ["add-assembly"]
        args.extend(["--target", self.target])
        args.extend(["--load", self.jb.load])
        args.append(genome)
        if not (self.target / genome.name).exists():
            jbrowse(args)

    def add_track_gff(self, species: str):
        gff = api.genome_gff3(species)
        ngff = gff.with_suffix("").with_suffix("").with_suffix(".name.gff3.gz")
        if ngff.exists():
            gff = ngff
        self.add_track(gff, trackid=gff.parent.name + ".gff3")
        self.text_index()

    def add_track(
        self, file: Path, category: str = "", trackid: str = "", subdir: str = ""
    ):
        # --description, --config
        args: subp.Args = ["add-track"]
        args.extend(["--target", self.target])
        args.extend(["--load", self.jb.load])
        if subdir:
            args.extend(["--subDir", subdir])
        if trackid:
            args.extend(["--trackId", trackid])
            self.tracks.append(trackid)
        if category:
            args.extend(["--category", category])
        if (csi := Path(str(file) + ".csi")).exists():
            args.extend(["--indexFile", csi])
        args.append(file)
        if not (self.target / subdir / file.name).exists():
            jbrowse(args)

    def text_index(self):
        # --attributes, --exclude, --file,  --perTrack, --tracks, --dryrun
        args: subp.Args = ["text-index"]
        args.extend(["--target", self.target])
        jbrowse(args)

    def set_default_session(self):
        args: subp.Args = ["set-default-session"]
        args.extend(["--target", self.target])
        args.extend(["--name", f"New {self.target.name} session"])
        args.extend(["--view", "LinearGenomeView"])
        patt = r"_inProm|_CE_genome-wide"  # redundant subsets
        patt += r"|_H\dK\d"
        patt += r"|SV_all-qin"
        rex = re.compile(patt)
        tracks = [x for x in self.tracks if not rex.search(x)]
        args.extend(["--tracks", ",".join(tracks)])
        jbrowse(args)

    def configure(self):
        chrom_sizes = {x[0]: x[1] for x in stat.chrom_sizes(self.target.name)}
        config_json = self.target / "config.json"
        with config_json.open() as fin:
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
            },
        ]
        for track in view["tracks"]:
            track["displays"] = [make_display(track)]
        if refnamealiases := self.make_refnamealiases():
            assembly["refNameAliases"] = refnamealiases
        cfg["configuration"] = make_configuration()
        with config_json.open("w") as fout:
            json.dump(cfg, fout, indent=2)

    def make_refnamealiases(self):
        species = self.target.name
        filename = f"{species}.chromAlias.txt"
        resources_alias = resources_data(filename)
        if not resources_alias.is_file():
            return None
        with (self.target / filename).open("w") as fout:
            fout.write(resources_alias.read_text())
        return {
            "adapter": {
                "type": "RefNameAliasAdapter",
                "location": {
                    "uri": filename,
                    "locationType": "UriLocation",
                },
            },
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
        if track["configuration"] not in clade_color:
            del item["constraints"]
    elif track["type"] == "AlignmentsTrack":
        item = {
            "type": "LinearPileupDisplay",
            "height": 20,
        }
    item["configuration"] = "-".join([track["configuration"], str(item["type"])])
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
        },
    }


def jbrowse(args: subp.Args):
    subp.run(["jbrowse", *args])


def npx_jbrowse(args: subp.Args, version: str = ""):
    pkg = "@jbrowse/cli"
    pkg = f"{pkg}@{version}" if version else pkg
    subp.run(["npx", pkg, *args])


def iter_targets(path: Path):
    for config_json in path.rglob("config.json"):
        if "test_data" in str(config_json):
            continue
        abs_config_json = config_json.resolve()
        if abs_config_json.parent.parent.name.startswith("jbrowse-"):
            yield abs_config_json.parent


def redirect_html(url: str):
    meta = f"""<meta http-equiv="refresh" content="1; URL={url}">"""
    return f"""<html><head>{meta}</head><body>Redirecting</body></html>"""


if __name__ == "__main__":
    main()
