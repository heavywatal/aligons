"""Chromatin Accessibility of Rice Tissues.

https://biobigdata.nju.edu.cn/cart/
"""

import logging
import re
import tarfile
from pathlib import Path

from aligons.extern import htslib
from aligons.util import cli, dl, fs, subp

from . import _rsrc, api

_log = logging.getLogger(__name__)
_HOST = "biobigdata.nju.edu.cn"
_files = [
    "Gene_Expression_Matrix_Log2.csv.gz",
    "NIP_ATAC_quantification_CPM.tsv.gz",
    "NIP_MH63_ZS97_bw.tar.gz",
    "NIP_MH63_ZS97_NIP_ref_ATAC_quantification_CPM.tsv.gz",
    "NIP_MH63_ZS97_Peaks.tar.gz",
    "NIP_MH63_ZS97_RNA.tar.gz",
]


def db_prefix(key: str = "") -> Path:
    res = api.prefix("cart")
    subdir = res / f"NIP_MH63_ZS97_{key}"
    if subdir.exists():
        return subdir
    return res


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-D", "--download", action="store_true")
    args = parser.parse_args(argv or None)
    if args.download:
        for f in _files:
            _retrieve_deploy(f)


def _retrieve_deploy(query: str) -> None:
    query_path = f"cart/static/download/{query}"
    url = f"https://{_HOST}/{query_path}"
    rawfile = _rsrc.db_root(_HOST) / query_path
    response = dl.fetch(url, rawfile)
    if query.endswith(".tar.gz"):
        _untar(response)


def _untar(response: dl.Response) -> None:
    tar_gz = response.path
    _log.info(tar_gz)
    dirname = tar_gz.with_suffix("").with_suffix("").name
    outdir = db_prefix()
    rex = re.compile(r"/NIP_[^/]+$")
    with tarfile.open(tar_gz) as tar:
        members = [x for x in tar.getmembers() if rex.search(x.name) is not None]
        if dirname not in members[0].name:
            outdir /= dirname
        to_extract: list[tarfile.TarInfo] = []
        for info in members:
            outfile = outdir / info.name
            if outfile.name.endswith(".bw.gz"):
                outfile = outfile.with_suffix("")
            if fs.is_outdated(outfile):
                _log.info(f"{info.name}: {info.size} -> {outfile}")
                to_extract.append(info)
        if to_extract and not cli.dry_run:
            tar.extractall(outdir, to_extract, filter="data")
        for info in to_extract:
            outfile = outdir / info.name
            outfile.chmod(0o644)
            if info.name.endswith(".bed.gz"):
                p = subp.run_zcat(outfile)
                htslib.bgzip(p.stdout, outfile)
                htslib.tabix(outfile)
            if info.name.endswith(".bw.gz"):
                fs.print_if_exists(_gunzip(outfile))


def _gunzip(infile: Path, *, if_: bool = True) -> Path:
    assert infile.suffix == ".gz", infile
    outfile = infile.with_suffix("")
    subp.run(["unzstd", "--rm", "-f", infile, "-o", outfile], if_=if_)
    return outfile


if __name__ == "__main__":
    main()
