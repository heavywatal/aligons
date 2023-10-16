import logging
import tomllib
from collections.abc import Mapping
from importlib import resources
from pathlib import Path
from types import MappingProxyType
from typing import Any, TypeAlias

import tomli_w

_log = logging.getLogger(__name__)

ConfDict: TypeAlias = MappingProxyType[str, Any]


def read_config(path: Path):
    with path.open("rb") as fin:
        update_nested(_config_src, tomllib.load(fin))


def log_config(path: Path = Path(".log.aligons.toml")):
    _log.info(f"{config}")
    if path.exists():
        with path.open("rb") as fin:
            reference = ConfDict(tomllib.load(fin))
        if _diff(reference, config):
            msg = f"config differs from the previous run: {path.absolute()}"
            raise ValueError(msg)
    else:
        with path.open("wb") as fout:
            tomli_w.dump(_config_src, fout)


def update_nested(x: dict[str, Any], other: Mapping[str, Any]):
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

_config_user = Path.home() / ".aligons.toml"
_config_pwd = Path(".aligons.toml")
for file in [_config_user, _config_pwd]:
    if file.exists():
        read_config(file)


config = ConfDict(_config_src)
empty_options = ConfDict({})


def _diff(lhs: Mapping[str, Any], rhs: Mapping[str, Any]):
    n = 0
    for key, lvalue in lhs.items():
        rvalue = rhs.get(key, {})
        if lvalue != rvalue:
            n += 1
            _log.warning(f"-[{key}] {lvalue}")
            _log.warning(f"+[{key}] {rvalue}")
    for key in set(rhs.keys()) - set(lhs.keys()):
        n += 1
        _log.warning(f"-[{key}] {{}}")
        _log.warning(f"+[{key}] {rhs[key]}")
    return n
