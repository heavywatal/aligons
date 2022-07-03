import importlib.resources as resources
import os.path
from pathlib import Path
from typing import Any, TypeAlias

import tomli

ConfDict: TypeAlias = dict[str, Any]


def read_config(path: Path):
    with open(path, "rb") as fin:
        update_nested(config, tomli.load(fin))


def update_nested(x: ConfDict, other: ConfDict):
    for key, value in other.items():
        if isinstance(x_val := x.get(key), dict):
            update_nested(x_val, value)  # type: ignore
        else:
            x[key] = value
    return x


def _expand_path(s: str):
    return Path(os.path.expandvars(s)).expanduser()


with resources.open_binary("biolopy.data", "config.toml") as fin:
    config: ConfDict = tomli.load(fin)

config["db"]["root"] = _expand_path(config["db"]["root"])
