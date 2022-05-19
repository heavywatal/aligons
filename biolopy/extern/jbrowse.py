"""
src: {ensemblgenomes.prefix}/fasta/{species}/dna/*.fa.gz
src: {ensemblgenomes.prefix}/gff3/{species}/*.chr.gff3.gz
src: ./pairwise/{species}/{query}/cram/genome.cram
src: ./multiple/{species}/{clade}/phastcons.bw
dst: ./jbrowse/{species}/config.json
"""
import json
import logging
import os
from pathlib import Path
from typing import Any

from ..db import ensemblgenomes
from ..util import cli, subp

_log = logging.getLogger(__name__)
_load = "symlink"
version = os.environ["JBROWSE_VERSION"]


def main(argv: list[str] = []):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-a", "--admin", action="store_true")
    parser.add_argument("-o", "--outdir", type=Path, default=Path(f"jbrowse-{version}"))
    parser.add_argument("multialign_species", type=Path)  # multiple/oryza_sativa/
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    jb = JBrowse(args.multialign_species, args.outdir)
    _log.info(f"{jb.target}")
    if args.admin:
        jb.admin_server()
        return
    jb.create()
    jb.add()
    jb.configure()


class JBrowse:
    def __init__(self, multialign_species: Path, root: Path):
        self.multiple_dir = multialign_species
        species_name = self.multiple_dir.name
        vNN_dir = multialign_species.parent.parent
        self.pairwise_dir = vNN_dir / "pairwise" / species_name
        self.root = root
        self.target = self.root / species_name
        self.genome = ensemblgenomes.get_file("*.genome.fa.gz", species_name)
        self.gff = ensemblgenomes.get_file("*.genome.gff3.gz", species_name)
        ngff = self.gff.with_suffix("").with_suffix("").with_suffix(".name.gff3.gz")
        if ngff.exists():
            self.gff = ngff

    def create(self):
        if not self.root.exists():
            jbrowse(["create", self.root])
        with open(self.root / "version.txt", "r") as fin:
            assert fin.read().strip() == version

    def add(self, force: bool = False):
        # if self.target.exists() and not force:
        #     return
        self.target.mkdir(0o755, parents=True, exist_ok=True)
        self.add_assembly()
        tracks: list[str] = []
        tracks.append(track_id := self.gff.parent.name)
        self.add_track(self.gff, id=track_id)
        self.text_index()
        for clade in self.multiple_dir.iterdir():
            wig = clade / "phastcons.bw"
            if wig.exists():
                tracks.append(track_id := clade.name)
                self.add_track(wig, "conservation", id=track_id, subdir=track_id)
        for query in self.pairwise_dir.iterdir():
            cram = query / "cram" / "genome.cram"
            if cram.exists():
                tracks.append(track_id := query.name)
                self.add_track(cram, "alignment", id=track_id, subdir=track_id)
        self.set_default_session(tracks)

    def add_assembly(self):
        """--alias, --name, --displayName"""
        args: subp.Args = ["add-assembly"]
        args.extend(["--target", self.target])
        args.extend(["--load", _load])
        args.append(self.genome)
        jbrowse(args, not (self.target / self.genome.name).exists())

    def add_track(self, file: Path, category: str = "", id: str = "", subdir: str = ""):
        """--description, --config"""
        args: subp.Args = ["add-track"]
        args.extend(["--target", self.target])
        args.extend(["--load", _load])
        if subdir:
            args.extend(["--subDir", subdir])
        if id:
            args.extend(["--trackId", id])
        if category:
            args.extend(["--category", category])
        args.append(file)
        jbrowse(args, not (self.target / subdir).exists())

    def text_index(self):
        """--attributes, --exclude, --file,  --perTrack, --tracks, --dryrun"""
        args: subp.Args = ["text-index"]
        args.extend(["--target", self.target])
        jbrowse(args)

    def set_default_session(self, tracks: list[str]):
        args: subp.Args = ["set-default-session"]
        args.extend(["--target", self.target])
        args.extend(["--name", f"New {self.target.name} session"])
        args.extend(["--view", "LinearGenomeView"])
        args.extend(["--tracks", ",".join(tracks)])
        jbrowse(args)

    def admin_server(self):
        jbrowse(["admin-server", "--root", self.root])

    def upgrade(self):
        jbrowse(["upgrade", self.root])

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


if __name__ == "__main__":
    main()
