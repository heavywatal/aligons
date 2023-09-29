import gzip
from pathlib import Path
from zipfile import ZipFile

import pytest
from aligons.util import fs


def test_is_outdated(tmp_path: Path):
    tmp_file = tmp_path / "newfile.txt"
    assert not tmp_file.exists()
    assert fs.is_outdated(tmp_file)
    tmp_file.open("w").close()
    assert tmp_file.exists()
    assert fs.is_outdated(tmp_file)
    this_file = Path(__file__)
    caches = list((this_file.parent / "__pycache__").iterdir())
    assert caches
    assert fs.is_outdated(this_file, caches)
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


def test_symlink(tmp_path: Path):
    newlink = tmp_path / "link"
    assert fs.symlink(tmp_path, newlink) == newlink
    assert newlink.exists()
    assert newlink.is_symlink()
    assert fs.symlink(tmp_path, newlink) == newlink
    assert len(list(tmp_path.iterdir())) == 1
    broken = fs.symlink(tmp_path / "noexist", tmp_path / "broken")
    assert broken.is_symlink()
    assert not broken.exists()
    assert fs.symlink(newlink, broken).exists()


def test_compress(tmp_path: Path):
    content = b"hello\n"
    gz_file = tmp_path / "plain.txt.gz"
    with gzip.open(gz_file, "wb") as fout:
        fout.write(content)
    with gz_file.open("rb") as fin:
        gz_content = fin.read()
    # contains filename
    assert gz_content != gzip.compress(content)
    assert fs.is_gz(gz_content)
    assert not fs.is_zip(gz_content)
    decompressed = fs.gzip_decompress(gz_content)
    assert not fs.is_gz(decompressed)
    assert decompressed == content
    assert decompressed == fs.gzip_decompress(decompressed)
    compressed = fs.gzip_compress(decompressed)
    assert fs.is_gz(compressed)
    assert gzip.decompress(compressed) == content
    assert compressed == fs.gzip_compress(compressed)


def test_zipfile(tmp_path: Path):
    content = b"hello\n"
    txt_file = tmp_path / "file.txt"
    with txt_file.open("wb") as fout:
        fout.write(content)
    zip_file = txt_file.with_suffix(txt_file.suffix + ".zip")
    with ZipFile(zip_file, "w") as zf:
        zf.write(txt_file)
    with zip_file.open("rb") as fin:
        zip_content = fin.read()
    assert fs.is_zip(zip_content)
    assert not fs.is_gz(zip_content)
    decompressed = fs.zip_decompress(zip_content)
    assert not fs.is_zip(decompressed)
    assert decompressed == content


def test_chdir(tmp_path: Path):
    assert Path.cwd() != tmp_path
    with fs.chdir(tmp_path):
        assert Path.cwd() == tmp_path
    assert Path.cwd() != tmp_path


def test_checksums(caplog: pytest.LogCaptureFixture):
    data_dir = Path(__file__).parent / "data"
    fs.checksums(data_dir / "CHECKSUMS")
    assert not caplog.text
    fs.checkline("42 1 sorted.gff3", data_dir)
    assert "expected: 42" in caplog.text
