"""Download genomes from NCBI.

<https://www.ncbi.nlm.nih.gov/datasets/docs/>
"""

import hashlib
import json
import re
import zipfile
from pathlib import Path

from aligons.extern import htslib
from aligons.util import cli, fs, logging, subp

from . import _rsrc, api, tools

_log = logging.getLogger(__name__)
_HOST = "ftp.ncbi.nlm.nih.gov"


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("-C", "--check", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-M", "--mask", type=str, help="taxon for RepeatMasker")
    parser.add_argument("accession", nargs="*")
    args = parser.parse_args(argv or None)
    if args.download:
        for acc in args.accession:
            _download_genome(acc)
    if args.check:
        for file in _prefix_mirror().glob("*.zip"):
            _check_zip(file)
    if args.mask:
        fts_fa = [_index_genome_fa(x) for x in _prefix_mirror().glob("*.zip")]
        fts: list[cli.Future[Path]] = []
        for ft in cli.as_completed(fts_fa):
            if args.mask:
                masked = tools.softmask(ft.result(), args.mask)
                fts.extend(tools.genome_to_twobits(masked))
            else:
                fs.print_if_exists(ft.result())
        cli.wait_raise(fts)


def db_prefix() -> Path:
    return api.prefix("ncbi")


def _prefix_mirror() -> Path:
    return _rsrc.db_root(_HOST) / "genomes" / "all"


def _download_genome(accession: str, *, if_: bool = True) -> None:
    """Download a genome data package via `datasets` CLI.

    `--include`: genome, rna, protein, cds, gff3, gtf, gbff, seq-report, none

    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/
    """
    outdir = _prefix_mirror()
    if not cli.dry_run:
        outdir.mkdir(0o755, parents=True, exist_ok=True)
    outfile = outdir / f"{accession}.zip"
    args: subp.Args = ["datasets", "download", "genome"]
    args.extend(["accession", accession])
    args.extend(["--filename", outfile])
    args.extend(["--include", "genome,gff3"])
    subp.run(args, if_=if_ and fs.is_outdated(outfile))


def _index_genome_fa(infile: Path) -> cli.Future[Path]:
    """Make indexed genome FASTA from a genome data package zip."""
    info = _get_info(infile)
    _log.debug(info)
    species: str = info["species"]  # pyright: ignore[reportAssignmentType]
    species = species.replace(" ", "_").lower()
    version = info["version"]
    stem = f"{species}_{version}"
    name = info["sequences"][0]
    out_fa = db_prefix() / species / f"{stem}.dna.genome.fa.gz"
    return cli.thread_submit(_index_bgzip, infile, name, out_fa)


def _index_bgzip(infile: Path, name: str, outfile: Path) -> Path:
    """Make bgzip and index of a specific file inside a zip archive."""
    if fs.is_outdated(outfile, infile):
        if not cli.dry_run:
            outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        with subp.popen_zcat([infile, Path(name)]) as zcat:
            htslib.bgzip(zcat.stdout, outfile)
    htslib.try_index(outfile)
    return outfile


def _get_info(file: Path) -> dict[str, str | list[str]]:
    """Extract metadata from a genome package zip."""
    res: dict[str, str | list[str]] = {}
    try:
        with zipfile.ZipFile(file) as zf:
            for name in zf.namelist():
                _log.debug(name)
            sequences = [x for x in zf.namelist() if "_genomic.f" in x]
            res["accession"], res["version"] = _get_asm_version(sequences[0])
            res["sequences"] = sequences
            with zf.open("ncbi_dataset/data/assembly_data_report.jsonl") as fin:
                for line in fin:
                    data = json.loads(line)
                    organism = data["organism"]
                    _log.debug(f"{organism}")
                    res["species"] = organism["organismName"]
            annotation = [x for x in zf.namelist() if ".gff" in x]
            if annotation:
                res["annotation"] = annotation[0]
    except zipfile.BadZipFile:
        _log.error(f"Bad zip file: {file}")
    return res


def _get_asm_version(name: str) -> tuple[str, str]:
    """Extract accession and assembly version from a file name inside a zip."""
    accession, name = name.removeprefix("ncbi_dataset/data/").split("/", 1)
    x = name.removeprefix(f"{accession}_")
    return accession, re.sub(r"_genomic\..+", "", x)


def _check_zip(file: Path) -> None:
    """Check integrity of a zip file and its MD5 checksums."""
    try:
        with zipfile.ZipFile(file) as zf:
            if zf.testzip() is not None:
                _log.error(f"Corrupted zip file: {file}")
            _check_md5(zf)
    except zipfile.BadZipFile:
        _log.error(f"Bad zip file: {file}")


def _check_md5(zf: zipfile.ZipFile, sum_file: str = "md5sum.txt") -> None:
    """Check MD5 checksums inside a zip file."""
    with zf.open(sum_file) as sum_fin:
        for line in sum_fin:
            hex_exp, path = line.decode().strip().split(None, 1)
            with zf.open(path) as fin:
                md5_obs = hashlib.md5(fin.read())  # noqa: S324
                if md5_obs.hexdigest() != hex_exp:
                    _log.error(f"MD5 mismatch: {path}")


if __name__ == "__main__":
    main()
