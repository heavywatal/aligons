import os
import tomllib
from pathlib import Path
from typing import TYPE_CHECKING, TypedDict

if TYPE_CHECKING:
    from collections.abc import Iterator
    from importlib.resources.abc import Traversable

from aligons.util import config, resources_data


class DataSet(TypedDict):
    url_prefix: str
    species: str
    version: str
    draft: bool
    label: str
    clade: str
    sequences: list[str]
    annotation: str


def iter_builtin_dataset(filename: str) -> Iterator[DataSet]:
    yield from iter_dataset(resources_data(filename))


def iter_dataset(toml: Path | Traversable) -> Iterator[DataSet]:
    with toml.open("rb") as fin:
        meta = tomllib.load(fin)
    for dic in meta["dataset"]:
        if dic.get("draft", False):
            continue
        yield DataSet(dic)


def db_root(relpath: str | Path = "") -> Path:
    return _expand_path(config["db"]["root"]) / relpath


def _expand_path(s: str) -> Path:
    return Path(os.path.expandvars(s)).expanduser()


# https://www.htslib.org/doc/samtools.html#REFERENCE_SEQUENCES
_ref_cache = str(db_root("hts-ref/%2s/%2s/%s"))
os.environ["REF_CACHE"] = _ref_cache
os.environ["REF_PATH"] = _ref_cache
