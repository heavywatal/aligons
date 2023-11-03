import logging
import re
from collections.abc import Iterable, Iterator
from pathlib import Path
from xml.etree import ElementTree

from aligons.util import cli, config, dl, tomli_w

from . import _rsrc, api

_log = logging.getLogger(__name__)
_HOST = "genome.jgi.doe.gov"

session = dl.LazySession(
    "https://signon.jgi.doe.gov/signon/create",
    data=config["jgi"].get("auth", {}),
)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    _args = parser.parse_args(argv or None)
    dataset_toml()


def iter_dataset(queries: Iterable[str]) -> Iterator[_rsrc.DataSet]:
    yet_to_find = set(queries)
    for entry in _rsrc.iter_dataset(dataset_toml()):
        if queries:
            try:
                yet_to_find.remove(entry["species"])
                yield entry
            except KeyError:
                pass
        else:
            yield entry
    if yet_to_find:
        msg = f"{yet_to_find} not found in {dataset_toml()}"
        raise ValueError(msg)


def dataset_toml(organism: str = config["jgi"]["organism"]) -> Path:
    outfile = api.prefix(_HOST) / (organism + ".toml")
    if not outfile.exists():
        dataset = {"dataset": list(iter_dataset_xml(organism))}
        if not cli.dry_run:
            outfile.parent.mkdir(0o755, exist_ok=True)
            with outfile.open("wb") as fout:
                tomli_w.dump(dataset, fout)
    _log.info(f"{outfile}")
    return outfile


def iter_dataset_xml(organism: str) -> Iterable[dict[str, str | list[str]]]:
    xml = fetch_xml(organism).path
    tree = ElementTree.parse(xml)  # noqa: S314
    for efolder in iter_species_folder(tree.getroot()):
        d = as_dict(efolder)
        if d:
            yield d


def as_dict(folder: ElementTree.Element) -> dict[str, str | list[str]]:
    _log.debug(f"{folder.attrib}")
    try:
        elem_assem = next(finditer(r"softmasked\.fa\.gz$", folder, "filename"))
        elem_annot = next(finditer(r"gene\.gff3?\.gz$", folder, "filename"))
    except StopIteration:
        return {}
    fa_url = _simplify_url(elem_assem.attrib["url"])
    gff_url = _simplify_url(elem_annot.attrib["url"])
    species, version = parse_filename_fa(Path(fa_url).name)
    return {
        "url_prefix": f"https://{_HOST}",
        "species": species,
        "version": version,
        "label": species,
        "clade": "",
        "sequences": [fa_url],
        "annotation": gff_url,
    }


def _simplify_url(url: str) -> str:
    return url.split("&url=", 1)[1]


def iter_species_folder(root: ElementTree.Element) -> Iterator[ElementTree.Element]:
    _log.info(f"{root.attrib}")
    for elem in root.iter("folder"):
        if elem.attrib["name"].split("_", 1)[0].istitle():
            yield elem


def parse_filename_fa(name: str) -> list[str]:
    stem = name.split(".softmasked")[0]
    return stem.split("_", 1)


def finditer(
    pattern: str, folder: ElementTree.Element, attrib: str
) -> Iterator[ElementTree.Element]:
    for elem in folder.iter("file"):
        if re.search(pattern, elem.attrib[attrib]):
            yield elem


def fetch_xml(organism: str) -> dl.Response:
    outfile = _rsrc.db_root(_HOST) / (organism + ".xml")
    query = f"?organism={organism}"
    url = f"https://{_HOST}/portal/ext-api/downloads/get-directory"
    return session.fetch(url + query, outfile)


def shorten(species: str) -> str:
    (genus, specific_epithet) = species.replace("_", " ", 1).split(maxsplit=1)
    return f"{genus[0]}{specific_epithet}"


if __name__ == "__main__":
    main()
