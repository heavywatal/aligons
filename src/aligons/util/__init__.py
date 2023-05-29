try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib
from importlib import resources
from pathlib import Path
from types import MappingProxyType
from typing import Any, TypeAlias

ConfDict: TypeAlias = MappingProxyType[str, Any]


def read_config(path: Path):
    with path.open("rb") as fin:
        update_nested(_config_src, tomllib.load(fin))


def update_nested(x: dict[str, Any], other: dict[str, Any]):
    for key, value in other.items():
        if isinstance(x_val := x.get(key), dict):
            update_nested(x_val, value)  # type: ignore[reportUnknownArgumentType]
        else:
            x[key] = value
    return x


def resources_data(child: str = ""):
    return resources.files("aligons.data").joinpath(child)


with resources_data("config.toml").open("rb") as fin:
    _config_src: dict[str, Any] = tomllib.load(fin)

config = ConfDict(_config_src)
empty_options = ConfDict({})
