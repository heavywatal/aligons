import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

import pytest

from aligons.util import config, log_config, update_config_if_exists, update_nested


def test_read_config(tmp_path: Path):
    assert config["db"]["root"]
    file = tmp_path / "aligons.toml"
    with file.open("wt") as fout:
        fout.write(f"[db]\nroot = '{tmp_path}'\n")
    update_config_if_exists(file)
    assert config["db"]["root"] == str(tmp_path)


def test_update_nested():
    x = {"both": 1, "x": 1}
    y = {"both": 2, "y": 2}
    update_nested(x, y)
    assert x == {"both": 2, "x": 1, "y": 2}
    assert y == {"both": 2, "y": 2}


def test_write_config(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    file = tmp_path / "aligons-log.toml"
    log_config(file)
    assert "+" not in caplog.text
    log_config(file)
    assert "+" not in caplog.text
    with file.open("r") as fin:
        content = fin.read()
    content = re.sub(r"(inner) = \d+", r"\1 = 999", content)
    content = re.sub(r"\[multiz\][^[]+", "\n[MISSING]\nx = 999\n\n", content)
    modified = tmp_path / "aligons-mod.toml"
    with modified.open("w") as fout:
        fout.write(content)
    with pytest.raises(ValueError, match="config"):
        log_config(modified)
    assert "'inner': 999" in caplog.text
    assert "-[multiz] {}" in caplog.text
    assert "+[multiz] {" in caplog.text
    assert "-[MISSING] {'x': 999}" in caplog.text
    assert "+[MISSING] {}" in caplog.text
