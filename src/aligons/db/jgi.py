import logging
import re
import xml.etree.ElementTree as ET
from collections.abc import Iterable, Iterator
from pathlib import Path

from aligons.util import cli, config, dl, fs, tomli_w

from . import _rsrc, api, phylo, tools

_log = logging.getLogger(__name__)
_HOST = "genome-downloads.jgi.doe.gov"

session = dl.LazySession(
    "https://signon.jgi.doe.gov/signon/create",
    data=config["jgi"].get("auth", {}),
)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    args = parser.parse_args(argv or None)
    if args.download:
        download()
    else:
        list(_iter_available())


def prefix() -> Path:
    return api.prefix(config["jgi"]["organism"].replace("V", "-").lower())


def download() -> None:
    queries: set[str] = {
        # "arabidopsis_lyrata",  # v2 but scaffold
        "sorghum_bicolor",
        "solanum_lycopersicum",
        "solanum_tuberosum",
    }
    fts: list[cli.Future[Path]] = []
    for entry in _iter_available():
        if entry["species"] not in queries:
            continue
        ft_fa, ft_gff = tools.fetch_and_bgzip(entry, prefix())
        fts.append(ft_gff)
        fts.extend(tools.genome_to_twobits(ft_fa))
    cli.wait_raise(fts)


def _iter_available() -> Iterator[_rsrc.DataSet]:
    nicknames = {shorten(long): long for long in phylo.list_species()}
    found: set[str] = set()
    for entry in _rsrc.iter_dataset(_dataset_toml()):
        if species := nicknames.get(entry["species"], None):
            entry["species"] = species
            if species in found:
                _log.debug(f"# {entry['label']}")
                continue
            _log.debug(f"{entry['label']}")
            found.add(species)
            yield entry
    not_found = set(nicknames.values()) - found
    _log.debug(f"{found = }")
    _log.debug(f"{not_found = }")


def _dataset_toml(organism: str = config["jgi"]["organism"]) -> Path:
    outfile = prefix() / (organism + ".toml")
    if not outfile.exists():
        xml = _fetch_xml(organism)
        dataset = {"dataset": list(_iter_dataset_xml(xml, organism))}
        if not cli.dry_run:
            outfile.parent.mkdir(0o755, exist_ok=True)
            with outfile.open("wb") as fout:
                tomli_w.dump(dataset, fout)
    return fs.print_if_exists(outfile)


def _iter_dataset_xml(xml: Path, organism: str) -> Iterable[dict[str, str | list[str]]]:
    tree = ET.parse(xml)  # noqa: S314
    xpath = f"folder[@name='{organism}']/folder"
    for sp_folder in tree.iterfind(xpath):
        _log.debug(sp_folder.attrib)
        for ver_folder in sp_folder.iterfind("folder"):
            d = _as_dict(ver_folder)
            if d:
                yield d


def _as_dict(folder: ET.Element) -> dict[str, str | list[str]]:
    try:
        elem_assembly = next(_finditer(r"softmasked\.fa\.gz$", folder, "filename"))
        elem_annot = next(_finditer(r"gene\.gff3?\.gz$", folder, "filename"))
    except StopIteration:
        return {}
    fa_url = _simplify_url(elem_assembly.attrib["url"])
    gff_url = _simplify_url(elem_annot.attrib["url"])
    species, rev = _parse_filename_gff(Path(gff_url).name)
    version = folder.get("name") or rev
    return {
        "url_prefix": f"https://{_HOST}",
        "species": species,
        "version": version,
        "label": f"{species}_{version}",
        "clade": "",
        "sequences": [fa_url],
        "annotation": gff_url,
    }


def _simplify_url(url: str) -> str:
    return url.split("&url=", 1)[1]


def _parse_filename_gff(name: str) -> list[str]:
    stem = name.split(".gene.gff")[0]
    return stem.split("_", 1)


def _finditer(pattern: str, folder: ET.Element, attrib: str) -> Iterator[ET.Element]:
    for elem in folder.iter("file"):
        if re.search(pattern, elem.attrib[attrib]):
            yield elem


def _fetch_xml(organism: str) -> Path:
    outfile = _rsrc.db_root(_HOST) / (organism + ".xml")
    all_xml = _fetch_xml_impl("Phytozome")
    if fs.is_outdated(outfile, all_xml):
        tree = ET.parse(all_xml)  # noqa: S314
        xpath = f"folder[@name!='{organism}']"
        for e in tree.getroot().findall(xpath):
            tree.getroot().remove(e)
        tree.write(outfile)
    return outfile


def _fetch_xml_impl(organism: str) -> Path:
    outfile = _rsrc.db_root(_HOST) / (organism + ".xml")
    query = f"?organism={organism}"
    url = f"https://{_HOST}/portal/ext-api/downloads/get-directory"
    return session.fetch(url + query, outfile).path


def shorten(species: str) -> str:
    (genus, specific_epithet) = species.split("_", maxsplit=1)
    return f"{genus[0].upper()}{specific_epithet}"


if __name__ == "__main__":
    main()
