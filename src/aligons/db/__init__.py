import os.path
from collections.abc import Iterator
from importlib.abc import Traversable
from pathlib import Path
from typing import TypedDict

from aligons.util import config, resources_data, tomllib


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


def _expand_path(s: str):
    return Path(os.path.expandvars(s)).expanduser()


def path(relpath: str | Path = ""):
    return _expand_path(config["db"]["root"]) / relpath


def path_mirror(relpath: str | Path = ""):
    return _expand_path(config["db"]["mirror"]) / relpath
