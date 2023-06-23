from pathlib import Path

from aligons.util import config, read_config, update_nested


def test_read_config(tmp_path: Path):
    assert config["db"]["root"]
    file = tmp_path / "aligons.toml"
    with file.open("wt") as fout:
        fout.write(f"[db]\nroot = '{tmp_path}'\n")
    read_config(file)
    assert config["db"]["root"] == str(tmp_path)


def test_update_nested():
    x = {"both": 1, "x": 1}
    y = {"both": 2, "y": 2}
    update_nested(x, y)
    assert x == {"both": 2, "x": 1, "y": 2}
    assert y == {"both": 2, "y": 2}
