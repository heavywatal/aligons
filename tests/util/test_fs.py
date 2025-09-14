from pathlib import Path

import pytest

from aligons.util import fs

data_dir = Path(__file__).parent.parent / "_data"


def test_is_outdated(tmp_path: Path):
    tmp_file = tmp_path / "new_file.txt"
    assert not tmp_file.exists()
    assert fs.is_outdated(tmp_file)
    tmp_file.open("w").close()
    assert tmp_file.exists()
    assert fs.is_outdated(tmp_file)
    this_file = Path(__file__)
    caches = list((this_file.with_name("__pycache__")).iterdir())
    assert caches
    assert fs.is_outdated(this_file, caches)
    assert not fs.is_outdated(this_file)


def test_sorted_naturally():
    names = ["chr.1.fa", "chr.2.fa", "chr.10.fa", "chr.Mt.fa"]
    files = [Path(x) for x in names]
    assert sorted(files) != files
    assert fs.sorted_naturally(names) == names
    assert fs.sorted_naturally(files) == files
    assert fs.sorted_naturally(reversed(files)) == files


def test_relpath(tmp_path: Path):
    sib = fs.relpath(tmp_path / "sib", tmp_path / "noexist")
    assert sib == Path("sib")
    parent = fs.relpath(tmp_path, tmp_path / "noexist")
    assert parent == Path("..")
    child = fs.relpath(tmp_path / "noexist", tmp_path)
    assert child == Path("noexist")


def test_symlink(tmp_path: Path):
    new_link = tmp_path / "link"
    assert fs.symlink(tmp_path, new_link) == new_link
    assert new_link.exists()
    assert new_link.is_symlink()
    assert fs.symlink(tmp_path, new_link) == new_link
    assert len(list(tmp_path.iterdir())) == 1
    broken = fs.symlink(tmp_path / "noexist", tmp_path / "broken")
    assert broken.is_symlink()
    assert not broken.exists()
    assert fs.symlink(new_link, broken).exists()
    rel_sib = fs.symlink(tmp_path / "noexist", tmp_path / "rel_sib", relative=True)
    assert rel_sib.readlink() == Path("noexist")
    rel_parent = fs.symlink(tmp_path, tmp_path / "rel_parent", relative=True)
    assert rel_parent.readlink() == Path("..")


def test_expect_suffix():
    fs.expect_suffix(Path("hello.txt"), ".txt")
    with pytest.raises(ValueError, match=r"expected suffix is \.gz"):
        fs.expect_suffix(Path("hello.txt"), ".gz")
    with pytest.raises(ValueError, match=r"unexpected suffix \.gz"):
        fs.expect_suffix(Path("hello.txt.gz"), ".gz", negate=True)


def test_chdir(tmp_path: Path):
    assert Path.cwd() != tmp_path
    with fs.chdir(tmp_path):
        assert Path.cwd() == tmp_path
    assert Path.cwd() != tmp_path


def test_checksums(caplog: pytest.LogCaptureFixture):
    fs.checksums(data_dir / "CHECKSUMS")
    assert not caplog.text
    fs.checkline("42 1 sorted.gff3", data_dir)
    assert "expected: 42" in caplog.text
