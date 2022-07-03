import pytest
from pathlib import Path

from aligons.util import fs


def test_is_outdated(tmp_path: Path):
    tmp_file = tmp_path / "newfile.txt"
    assert fs.is_outdated(tmp_file) and not tmp_file.exists()
    open(tmp_file, "w").close()
    assert fs.is_outdated(tmp_file) and tmp_file.exists()
    this_file = Path(__file__)
    caches = list((this_file.parent / "__pycache__").iterdir())
    assert fs.is_outdated(this_file, caches) and caches
    assert not fs.is_outdated(this_file)


def test_sorted_naturally(capsys: pytest.CaptureFixture[str]):
    names = ["chr.1.fa", "chr.2.fa", "chr.10.fa", "chr.Mt.fa"]
    files = [Path(x) for x in names]
    assert sorted(files) != files
    assert fs.sorted_naturally(names) == names
    assert fs.sorted_naturally(files) == files
    assert fs.sorted_naturally(reversed(files)) == files
    fs.main(names)
    captured = capsys.readouterr()
    assert captured.out == "\n".join(names) + "\n"
