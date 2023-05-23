import importlib.resources as resources
from pathlib import Path
from typing import Any, TypeAlias

import tomllib

ConfDict: TypeAlias = dict[str, Any]


def read_config(path: Path):
    with open(path, "rb") as fin:
        update_nested(config, tomllib.load(fin))


def update_nested(x: ConfDict, other: ConfDict):
    for key, value in other.items():
        if isinstance(x_val := x.get(key), dict):
            update_nested(x_val, value)  # type: ignore
        else:
            x[key] = value
    return x


with resources.files("aligons.data").joinpath("config.toml").open("rb") as fin:
    config: ConfDict = tomllib.load(fin)
