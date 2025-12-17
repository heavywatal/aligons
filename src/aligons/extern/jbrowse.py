"""The next-generation genome browser.

<https://jbrowse.org>

src: {vNN}/pairwise/{species}/{query}/cram/genome.cram
src: {vNN}/conservation/{species}/{clade}/phastcons.bw
dst: {document_root}/{jbrowse_XYZ}/{vNN}/{species}
"""

import json
import logging
import re
from pathlib import Path
from typing import Any

from aligons.db import api, cart, phylo, plantdhs, plantregmap, riceencode
from aligons.util import cli, config, fs, resources_data, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    """CLI to manage JBrowse instances."""
    parser = cli.ArgumentParser()
    parser.add_argument("-a", "--admin", action="store_true")
    parser.add_argument("-u", "--upgrade", action="store_true")
    parser.add_argument("-l", "--list", action="store_true")
    parser.add_argument("indir", type=Path, help="e.g., conservation/oryza_sativa/")
    args = parser.parse_args(argv or None)
    jb = JBrowse()
    if args.admin:
        jb.admin_server(args.indir)
        return
    if args.upgrade:
        jb.upgrade()
        return
    if args.list:
        for d in args.indir.iterdir():
            jb.print_url(d)
        return
    jb.config(args.indir)


class JBrowse:
    """Wrapper for JBrowse CLI."""

    def __init__(self) -> None:
        """Initialize version and root directory.

        - `slug`: Relative URL path from the root: `jbrowse-{version}`
        - `root`: Root directory of a JBrowse instance: `{document_root}/{slug}/`
        """
        self.version = self._version()
        document_root = Path(config["jbrowse"]["document_root"]).expanduser()
        self.slug = f"jbrowse-{self.version}"
        self.root = document_root / self.slug

    def admin_server(self, indir: Path) -> None:
        """Run `jbrowse admin-server` and print the URL."""
        jbc = JBrowseConfig(self.root, self.slug, indir)
        cmd = ["jbrowse", "admin-server", "--root", self.root]
        p = subp.Popen(cmd, stdout=subp.PIPE)
        assert p.stdout is not None
        for line in p.stdout:
            if mobj := re.search(rb"https?://\S+", line):
                url = mobj.group(0).decode()
                _log.info(f"{url}&config={jbc.relpath}/config.json")

    def upgrade(self) -> None:
        """Run `jbrowse upgrade`."""
        _jbrowse(["upgrade", self.root])

    def config(self, indir: Path) -> None:
        """Configure a JBrowse instance for a species."""
        self._create()
        jbc = JBrowseConfig(self.root, self.slug, indir)
        jbc.add()
        jbc.add_external()
        jbc.configure()
        jbc.set_default_session()
        jbc.write_redirect_html()

    def print_url(self, indir: Path) -> None:
        """Print the URL of the JBrowse instance for the given directory."""
        jbc = JBrowseConfig(self.root, self.slug, indir)
        if jbc.target.exists():
            _log.info(jbc.url)

    def _create(self) -> None:
        if not self.root.exists():
            args = ["create", self.root]
            args.append(f"--tag=v{self.version}")
            _jbrowse(args)
        with (self.root / "version.txt").open() as fin:
            version_txt = fin.read().strip()
            if version_txt != self.version:
                msg = f"{version_txt=} != {self.version=}"
                raise ValueError(msg)

    def _version(self) -> str:
        mobj = re.search(r"@jbrowse/cli/(\S+)", self._help())
        assert mobj, self._help()
        return mobj.group(1)

    def _help(self) -> str:
        p = _jbrowse(["help"], stdout=subp.PIPE, quiet=True)
        return p.stdout.decode()


class JBrowseConfig:
    """Configuration for a JBrowse instance."""

    def __init__(self, root: Path, slug: str, cons_species: Path) -> None:
        """Initialize with directory paths.

        :param root: Root directory of a JBrowse instance: `{document_root}/{slug}/`.
        :param slug: Relative URL path from the root: `jbrowse-{version}`.
        :param cons_species: Directory with phast results: `conservation/{species}`.
        """
        self.load = config["jbrowse"]["load"]
        self.cons_dir = cons_species
        self.species = self.cons_dir.name
        vnn_dir = cons_species.parent.parent.resolve()
        self.pairwise_dir = vnn_dir / "pairwise" / self.species
        self.relpath = Path(vnn_dir.name) / self.species
        self.target = root / self.relpath
        self.tracks: list[str] = []
        self.slug = slug
        base_url = config["jbrowse"]["base_url"].rstrip("/")
        self.url = f"{base_url}/{slug}/{self.relpath}/"
        _log.debug(self.target)

    def write_redirect_html(self) -> None:
        """Write `{target}/index.html` to redirect to the config URL."""
        url = f"/{self.slug}/?config={self.relpath}/config.json"
        if not cli.dry_run:
            with (self.target / "index.html").open("w") as fout:
                fout.write(_redirect_html(url))
        _log.info(f"{self.url} -> {url}")

    def add(self) -> None:
        """Add assemblies and basic tracks."""
        self.target.mkdir(0o755, parents=True, exist_ok=True)
        self._add_assembly()
        self._add_track_gff()
        clades = [x.name for x in self.cons_dir.iterdir() if "-" not in x.name]
        _log.info(clades)
        clades = phylo.sorted_by_len_newicks(clades, reverse=True)
        for clade in clades:
            self._add_phast(clade)
        iter_crai = self.pairwise_dir.rglob("*.cram.crai")
        crams = {cram.parent.parent.name: cram.with_suffix("") for cram in iter_crai}
        for query in phylo.list_species(clades[0]):
            if cram := crams.pop(query, None):
                self.add_track(cram, "alignment", trackid=query, subdir=query)
        for query, cram in crams.items():
            self.add_track(cram, "alignment", trackid=query, subdir=query)

    def _add_phast(self, clade: str) -> None:
        clade_dir = self.cons_dir / clade
        for wig in clade_dir.glob("*.bw"):
            stem = wig.stem
            self.add_track(wig, "conservation", trackid=f"{stem}-{clade}", subdir=clade)
        for csi in clade_dir.glob("*.csi"):
            bed = csi.with_suffix("")
            stem = bed.with_suffix("").stem
            self.add_track(bed, "conservation", trackid=f"{stem}-{clade}", subdir=clade)

    def add_external(self) -> None:
        """Add supplementary tracks from external databases."""
        if self.species == "oryza_sativa":
            _add_plantregmap(self, self.species)
            _add_papers_data(self)
            _add_plantdhs(self)
            _add_cart(self)
            _add_riceencode(self)
        elif self.species.endswith("_mh63") or self.species.endswith("_zs97"):
            _add_riceencode(self, self.species)
        elif self.species == "solanum_lycopersicum":
            _add_plantregmap(self, self.species)

    def _add_assembly(self) -> None:
        # --alias, --name, --displayName
        genome = api.genome_fa(self.species)
        name = genome.name.split(".dna_sm.", 1)[0]
        args: subp.Args = ["add-assembly", "--overwrite"]
        args.extend(["--name", name])
        args.extend(["--displayName", name])
        args.extend(["--target", self.target])
        args.extend(["--load", self.load])
        args.append(genome)
        if not (self.target / genome.name).exists():
            _jbrowse(args)

    def _add_track_gff(self) -> None:
        gff = api.genome_gff3(self.species)
        named = gff.with_suffix("").with_suffix("").with_suffix(".name.gff3.gz")
        if named.exists():
            gff = named
        self.add_track(gff, trackid=gff.parent.name + ".gff3")
        self._text_index()

    def add_track(
        self,
        file: Path,
        category: str = "",
        trackid: str = "",
        subdir: str = "",
        *,
        session: bool = True,
    ) -> None:
        """Call `jbrowse add-track`.

        :param file: The data file to add.
        :param category: Category name for the track.
        :param trackid: Track ID. The file name is used if empty.
        :param subdir: Subdirectory under the target directory.
        :param session: Add the track to the default session if True.
        """
        # --description, --config
        args: subp.Args = ["add-track"]
        args.extend(["--target", self.target])
        args.extend(["--load", self.load])
        if subdir:
            args.extend(["--subDir", subdir])
        if trackid:
            args.extend(["--trackId", trackid])
            if session:
                self.tracks.append(trackid)
        if category:
            args.extend(["--category", category])
        if (csi := Path(str(file) + ".csi")).exists():
            args.extend(["--indexFile", csi])
        args.append(file)
        if not (self.target / subdir / file.name).exists():
            _jbrowse(args)

    def _text_index(self) -> None:
        # --attributes, --exclude, --file,  --perTrack, --tracks, --dryrun
        args: subp.Args = ["text-index"]
        args.extend(["--target", self.target])
        _jbrowse(args)

    def set_default_session(self, session: Path | None = None) -> None:
        """Write `{target}/session.json` and call `jbrowse set-default-session`."""
        if session is None:
            obj = self._create_session_dict()
            session = self.target / "session.json"
            with session.open("w") as fout:
                json.dump(obj, fout, indent=2)
        args: subp.Args = ["set-default-session"]
        args.extend(["--target", self.target])
        args.extend(["--session", session])
        _jbrowse(args)

    def _select_tracks(self) -> list[str]:
        patt = r"_inProm|_CE_genome-wide"  # redundant subsets
        patt += r"|_H\dK\d"
        patt += r"|SV_all-qin"
        patt += r"|PhyloP|\.net$"
        patt += r"|_DNase$|_NPS$"
        patt += r"|^cns0|^most-cons|cns-poaceae|cns-monocot"
        patt += r"|^NIP_"
        rex = re.compile(patt)
        return [x for x in self.tracks if not rex.search(x)]

    def configure(self) -> None:
        """Write `{target}/config.json`."""
        config_json = self.target / "config.json"
        with config_json.open() as fin:
            cfg = json.load(fin)
        for track in cfg["tracks"]:
            display = _make_display(track["trackId"], track["adapter"]["type"])
            track["displays"] = [display]
        assembly = cfg["assemblies"][0]
        assembly["sequence"]["name"] = assembly["name"] + ".fa"
        assembly["sequence"]["displays"] = [_LinearGCContentDisplay(assembly["name"])]
        if refnamealiases := self._make_refnamealiases():
            assembly["refNameAliases"] = refnamealiases
        cfg["configuration"] = config["jbrowse"]["configuration"]
        with config_json.open("w") as fout:
            json.dump(cfg, fout, indent=2)

    def _create_session_dict(self) -> dict[str, Any]:
        config_json = self.target / "config.json"
        with config_json.open() as fin:
            cfg = json.load(fin)
        assembly = cfg["assemblies"][0]
        view = _create_view(self.target.name, assembly["name"])
        track_types: dict[str, str] = {}
        for track in cfg["tracks"]:
            track_types[track["trackId"]] = track["type"]
        view_tracks: list[dict[str, Any]] = []
        view_tracks.append(_session_track_refseq(assembly["name"]))
        for track_id in self._select_tracks():
            track_type = track_types[track_id]
            view_tracks.append(_session_track(track_id, track_type))
        view["tracks"] = view_tracks
        return {
            "name": f"New {self.target.name} session",
            "views": [view],
        }

    def _make_refnamealiases(self) -> dict[str, Any] | None:
        path = f"chromAlias/{self.species}.chromAlias.txt"
        resources_alias = resources_data(path)
        if not resources_alias.is_file():
            return None
        filename = Path(path).name
        with (self.target / filename).open("w") as fout:
            fout.write(resources_alias.read_text())
        fs.print_if_exists(self.target / filename)
        return {
            "adapter": {
                "type": "RefNameAliasAdapter",
                "location": {
                    "uri": filename,
                    "locationType": "UriLocation",
                },
            },
        }


def _make_display(track_id: str, adapter_type: str) -> dict[str, Any]:
    if adapter_type.startswith("Bed"):
        return _LinearBasicDisplay(track_id)
    if adapter_type.startswith("Gff"):
        if "gff3" in track_id:
            return _LinearBasicDisplay(track_id, 60, labels=True, desc=True)
        return _LinearBasicDisplay(track_id, 40)
    if adapter_type.startswith("BigWig"):
        return _LinearWiggleDisplay(track_id)
    if adapter_type.startswith("Cram"):
        return _LinearPileupDisplay(track_id)
    return {}


def _LinearGCContentDisplay(track_id: str, height: int = 20) -> dict[str, Any]:  # noqa: N802
    return _display("LinearGCContent", track_id, height) | {
        "minScore": 0,
        "maxScore": 1,
        "renderers": {
            "DensityRenderer": {
                "type": "DensityRenderer",
                "posColor": "#000",
            }
        },
    }


def _LinearBasicDisplay(  # noqa: N802
    track_id: str, height: int = 30, *, labels: bool = False, desc: bool = False
) -> dict[str, Any]:
    clade = track_id.rsplit("-", 1)[-1]
    return _display("LinearBasic", track_id, height) | {
        "renderer": _SvgFeatureRenderer(_palette_get(clade), labels=labels, desc=desc),
    }


def _LinearWiggleDisplay(track_id: str, height: int = 40) -> dict[str, Any]:  # noqa: N802
    item = _display("LinearWiggle", track_id, height)
    clade = track_id.rsplit("-", 1)[-1]
    item["renderers"] = {
        "XYPlotRenderer": {"type": "XYPlotRenderer", "color": _palette_get(clade)}
    }
    if "phast" in track_id.lower():
        item["constraints"] = {"max": 1, "min": 0}
    return item


def _LinearPileupDisplay(track_id: str, height: int = 20) -> dict[str, Any]:  # noqa: N802
    return _display("LinearPileup", track_id, height)


def _display(prefix: str, track_id: str, height: int) -> dict[str, Any]:
    return {
        "type": f"{prefix}Display",
        "displayId": f"{track_id}-{prefix}Display",
        "height": height,
    }


def _SvgFeatureRenderer(  # noqa: N802
    color1: str = "", *, labels: bool = False, desc: bool = False
) -> dict[str, Any]:
    item = {
        "type": "SvgFeatureRenderer",
        "showLabels": labels,
        "showDescriptions": desc,
    }
    if color1:
        item["color1"] = color1
    return item


def _session_track_refseq(asm_name: str) -> dict[str, Any]:
    track_type = "ReferenceSequenceTrack"
    display = _session_display(asm_name, track_type)
    display["rendererTypeNameState"] = "density"
    return {
        "type": track_type,
        "configuration": f"{asm_name}-{track_type}",
        "displays": [display],
    }


def _session_track(track_id: str, track_type: str) -> dict[str, Any]:
    return {
        "type": track_type,
        "configuration": track_id,
        "displays": [_session_display(track_id, track_type)],
    }


def _session_display(track_id: str, track_type: str) -> dict[str, Any]:
    display_map = {
        "ReferenceSequenceTrack": "LinearGCContentDisplay",
        "FeatureTrack": "LinearBasicDisplay",
        "QuantitativeTrack": "LinearWiggleDisplay",
        "AlignmentsTrack": "LinearPileupDisplay",
    }
    display_type = display_map[track_type]
    return {
        "type": display_type,
        "configuration": f"{track_id}-{display_type}",
    }


def _add_riceencode(jbc: JBrowseConfig, species: str = "") -> None:
    for fmt in riceencode.db_prefix().iterdir():
        for strain in fs.sorted_naturally(fmt.iterdir()):
            session = species.endswith(strain.name.lower())
            for histone in fs.sorted_naturally(strain.iterdir()):
                for path in fs.sorted_naturally(histone.iterdir()):
                    if path.suffix not in (".gz", ".bw"):
                        continue
                    category = f"RiceENCODE,{strain.name},{histone.name},{fmt.name}"
                    jbc.add_track(path, category, trackid=path.stem, session=session)


def _add_cart(jbc: JBrowseConfig) -> None:
    for path in fs.sorted_naturally(cart.db_prefix("bw").glob("NIP_*.bw")):
        jbc.add_track(path, "CART,BigWig", trackid=path.stem)
    for path in fs.sorted_naturally(cart.db_prefix("RNA").glob("NIP_*.bw")):
        jbc.add_track(path, "CART,RNA", trackid=path.stem)
    for path in fs.sorted_naturally(cart.db_prefix("Peaks").glob("NIP_*.bed.gz")):
        jbc.add_track(path, "CART,Peaks", trackid=path.stem)


def _add_plantdhs(jbc: JBrowseConfig) -> None:
    for path in fs.sorted_naturally(plantdhs.db_prefix().glob("Rice_*.bw")):
        trackid = path.stem.removeprefix("Rice_")
        jbc.add_track(path, "plantdhs", trackid=trackid)
    for path in fs.sorted_naturally(plantdhs.db_prefix().glob("*.gff.gz")):
        jbc.add_track(path, "plantdhs", trackid=path.stem)


def _add_plantregmap(jbc: JBrowseConfig, species: str) -> None:
    for path in fs.sorted_naturally(plantregmap.rglob("*.bw", species)):
        trackid = path.with_suffix(".bedGraph").name
        jbc.add_track(path, "plantregmap", trackid=trackid)
    for path in fs.sorted_naturally(plantregmap.rglob("*.cram", species)):
        trackid = path.with_suffix(".net").name
        jbc.add_track(path, "plantregmap", trackid=trackid)
    for path in fs.sorted_naturally(plantregmap.rglob("*.gff.gz", species)):
        trackid = re.sub(r"_[^_]+\.gff\.gz$", "", path.name)
        jbc.add_track(path, "plantregmap", trackid=trackid)
    for path in fs.sorted_naturally(plantregmap.rglob("*.bed.gz", species)):
        if "_work" in str(path):
            continue
        trackid = re.sub(r"(_normal)?\.bed\.gz$", "", path.name)
        jbc.add_track(path, "plantregmap", trackid=trackid)


def _add_papers_data(jbc: JBrowseConfig) -> None:
    for path in fs.sorted_naturally(api.prefix("papers").glob("*.bed.gz")):
        jbc.add_track(path, "papers", trackid=path.with_suffix("").stem)
    suzuemon = api.prefix("suzuemon")
    if (f := suzuemon / "sv_with_DEG.bed.gz").exists():
        jbc.add_track(f, "papers", trackid="SV_DEG-qin2021", subdir="suzuemon")
    if (f := suzuemon / "SV.bed.gz").exists():
        jbc.add_track(f, "papers", trackid="SV_all-qin2021", subdir="suzuemon")


def _jbrowse(args: subp.Args, **kwargs: Any) -> subp.CompletedProcess[Any]:
    return subp.run(["jbrowse", *args], **kwargs)


def _redirect_html(url: str) -> str:
    meta = f"""<meta http-equiv="refresh" content="1; URL={url}">"""
    return f"""<html><head>{meta}</head><body>Redirecting</body></html>"""


def _create_view(species: str, assembly: str) -> dict[str, Any]:
    asm_conf = _find_config_assembly(species)
    _log.info(asm_conf)
    location = asm_conf["location"]
    mobj = re.search(r"(\w+)\s*:\s*([\d,]+)\s*(?::|\.{2,})\s*([\d,]+)", location)
    assert mobj is not None, location
    chrom = mobj.group(1)
    start = int(mobj.group(2).replace(",", ""))
    end = int(mobj.group(3).replace(",", ""))
    chrom_sizes = api.chrom_sizes(species)
    region = {
        "refName": chrom,
        "start": 0,
        "end": chrom_sizes[chrom],
        "reversed": False,
        "assemblyName": assembly,
    }
    window_width = config["jbrowse"]["window_width"]
    bp_per_px = round((end - start) / window_width, 2)
    view_type = "LinearGenomeView"
    return {
        "id": f"{view_type}-0",
        "type": view_type,
        "bpPerPx": bp_per_px,
        "offsetPx": int(start / bp_per_px),
        "displayedRegions": [region],
        "tracks": [],
    }


def _find_config_assembly(species: str) -> dict[str, Any]:
    for entry in config["jbrowse"]["assemblies"]:
        if entry["species"] == species:
            return entry
    msg = f"{species} not in config.jbrowse.assemblies"
    _log.warning(msg)
    chrom_sizes = api.chrom_sizes(species)
    key, value = next(iter(chrom_sizes.items()))
    end = min(10000, value)
    return {"species": species, "location": f"{key}:1..{end}"}


def _palette_get(name: str, default: str = "#888") -> str:
    return config["jbrowse"]["palette"].get(name, default)


if __name__ == "__main__":
    main()
