import logging
import tomllib
from importlib import resources
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Mapping
    from importlib.resources.abc import Traversable

import tomli_w

_log = logging.getLogger(__name__)


def resources_data(child: str = "") -> Traversable:
    return resources.files("aligons.data").joinpath(child)


def update_config_if_exists(path: Path) -> None:
    if path.exists():
        update_nested(_config_src, _read_config(path))


def update_nested(x: dict[str, Any], other: Mapping[str, Any]) -> dict[str, Any]:
    for key, value in other.items():
        if isinstance(x_val := x.get(key), dict):
            update_nested(x_val, value)  # type: ignore[reportUnknownArgumentType]
        elif key == "assemblies" and isinstance(value, list):
            x[key].extend(value)
        else:
            x[key] = value
    return x


def _read_config(path: Path | Traversable) -> dict[str, Any]:
    with path.open("rb") as fin:
        return tomllib.load(fin)


_config_src: dict[str, Any] = _read_config(resources_data("config.toml"))
update_config_if_exists(Path.home() / ".aligons.toml")
update_config_if_exists(Path(".aligons.toml"))
config = MappingProxyType[str, Any](_config_src)


def log_config(path: Path = Path(".log.aligons.toml")) -> None:
    _log.info(config)
    if path.exists():
        reference = _read_config(path)
        if _diff(reference, config):
            msg = f"config differs from the previous run: {path.absolute()}"
            raise ValueError(msg)
    else:
        with path.open("wb") as fout:
            tomli_w.dump(config, fout)


def _diff(lhs: Mapping[str, Any], rhs: Mapping[str, Any]) -> int:
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
