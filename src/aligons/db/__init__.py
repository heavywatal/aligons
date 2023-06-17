import os.path
from pathlib import Path

from aligons.util import config


def _expand_path(s: str):
    return Path(os.path.expandvars(s)).expanduser()


def path(relpath: str | Path = ""):
    return _expand_path(config["db"]["root"]) / relpath


def path_mirror(relpath: str | Path = ""):
    return _expand_path(config["db"]["mirror"]) / relpath
