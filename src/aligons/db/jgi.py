import logging
import re
from pathlib import Path
from xml.etree import ElementTree

from aligons import db
from aligons.util import cli, config, dl, tomli_w

_log = logging.getLogger(__name__)
_HOST = "genome.jgi.doe.gov"

session = dl.LazySession(
    "https://signon.jgi.doe.gov/signon/create",
    data=config["jgi"].get("auth", {}),
)


def main(argv: list[str] | None = None):
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    _args = parser.parse_args(argv or None)
    dataset_toml()


def dataset_toml(organism: str = config["jgi"]["organism"]):
    outfile = db.path(_HOST) / (organism + ".toml")
    if not outfile.exists():
        xml = fetch_xml(organism).path
        tree = ElementTree.parse(xml)  # noqa: S314
        dataset = gather_dataset(tree.getroot())
        if not cli.dry_run:
            outfile.parent.mkdir(0o755, exist_ok=True)
            with outfile.open("wb") as fout:
                tomli_w.dump(dataset, fout)
    _log.info(f"{outfile}")
    return outfile


def gather_dataset(root: ElementTree.Element):
    datasets: list[dict[str, str | list[str]]] = []
    for efolder in iter_species_folder(root):
        d = as_dict(efolder)
        if d:
            datasets.append(d)
    return {"dataset": datasets}


def as_dict(efolder: ElementTree.Element) -> dict[str, str | list[str]]:
    _log.debug(f"{efolder.attrib}")
    elem_assem = search_filename(efolder, r"softmasked\.fa\.gz$")
    elem_annot = search_filename(efolder, r"repeatmasked.+\.gff3\.gz$")
    if elem_assem is None or elem_annot is None:
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


def _simplify_url(url: str):
    return url.split("&url=", 1)[1]


def iter_species_folder(root: ElementTree.Element):
    _log.info(f"{root.attrib}")
    for elem in root.iter("folder"):
        if re.match("^[A-Z]", elem.attrib["name"]):
            yield elem


def parse_filename_fa(name: str):
    stem = name.split(".softmasked")[0]
    return stem.split("_", 1)


def search_filename(
    folder: ElementTree.Element, pattern: str
) -> ElementTree.Element | None:
    for elem in folder.iter("file"):
        if re.search(pattern, elem.attrib["filename"]):
            return elem
    return None


def fetch_xml(organism: str):
    outfile = db.path_mirror(_HOST) / (organism + ".xml")
    query = f"?organism={organism}"
    url = f"https://{_HOST}/portal/ext-api/downloads/get-directory"
    return session.fetch(url + query, outfile)


if __name__ == "__main__":
    main()
