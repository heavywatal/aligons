import os.path
from pathlib import Path
from typing import TypedDict

from aligons.util import config


class DataSet(TypedDict):
    url_prefix: str
    species: str
    version: str
    draft: bool
    label: str
    clade: str
    sequences: list[str]
    annotation: str


def _expand_path(s: str):
    return Path(os.path.expandvars(s)).expanduser()


def path(relpath: str | Path = ""):
    return _expand_path(config["db"]["root"]) / relpath


def path_mirror(relpath: str | Path = ""):
    return _expand_path(config["db"]["mirror"]) / relpath
