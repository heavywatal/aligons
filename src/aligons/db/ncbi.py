"""Download genomes from NCBI.

<https://www.ncbi.nlm.nih.gov/datasets/docs/>
"""

import hashlib
import json
import re
import zipfile
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

from aligons.util import cli, fs, logging, subp

from . import _rsrc, api, tools

_log = logging.getLogger(__name__)
_HOST = "ftp.ncbi.nlm.nih.gov"


def main(argv: list[str] | None = None) -> None:
    """CLI for downloading and processing NCBI genomes."""
    parser = cli.ArgumentParser()
    parser.add_argument("-C", "--check", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-M", "--deploy", action="store_true")
    parser.add_argument("accession", nargs="*")
    args = parser.parse_args(argv or None)
    if args.download:
        for acc in args.accession:
            _download_genome(acc)
    if args.check:
        for file in _prefix_mirror().glob("*.zip"):
            _check_zip(file)
    if args.deploy:
        fts: list[cli.Future[Path]] = []
        for acc in args.accession:
            genome_zip = _download_genome(acc)
            fts.extend(_index_genome(genome_zip))
        fts2b: list[cli.Future[Path]] = []
        for ft in cli.as_completed(fts):
            if (res := ft.result()).name.endswith(".fa.gz"):
                fts2b.extend(tools.genome_to_twobits(res))
        cli.wait_raise(fts2b)


def db_prefix() -> Path:
    """Directory of preprocessed NCBI genomes."""
    return api.prefix("ncbi")


def _prefix_mirror() -> Path:
    return _rsrc.db_root(_HOST) / "genomes" / "all"


def _download_genome(accession: str) -> Path:
    """Download a genome data package via `datasets` CLI.

    `--include`: genome, rna, protein, cds, gff3, gtf, gbff, seq-report, none

    <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/>
    """
    outdir = _prefix_mirror()
    if not cli.dry_run:
        outdir.mkdir(0o755, parents=True, exist_ok=True)
    outfile = outdir / f"{accession}.zip"
    args: subp.Args = ["datasets", "download", "genome", "--no-progressbar"]
    args.extend(["accession", accession])
    args.extend(["--filename", outfile])
    args.extend(["--include", "genome,gff3"])
    subp.run(args, if_=fs.is_outdated(outfile))
    return outfile


def _index_genome(genome_zip: Path) -> list[cli.Future[Path]]:
    """Make indexed genome FASTA and GFF from a genome data package zip.

    > All genome sequences are softmasked using WindowMasker or RepeatMasker.

    <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/>

    :param genome_zip: NCBI genome package file: `{accession}.zip`.
    :returns: Future of the bgzipped files.
    """
    info = _get_info(genome_zip)
    _log.debug(info)
    species: str = info["species"]  # pyright: ignore[reportAssignmentType]
    species = species.replace(" ", "_").lower()
    version = info["version"]
    stem = f"{species}_{version}"
    name = info["sequences"][0]
    out_fa = db_prefix() / species / f"{stem}.dna_sm.genome.fa.gz"
    fts = [cli.thread_submit(tools.unzip_index_bgzip, genome_zip, name, out_fa)]
    if isinstance(gff := info.get("annotation"), str) and gff:
        _log.debug(f"{gff = }")
        out_gff = db_prefix() / species / f"{stem}.genome.gff3.gz"
        fts.append(cli.thread_submit(tools.unzip_index_bgzip, genome_zip, gff, out_gff))
    return fts


def _get_info(genome_zip: Path) -> dict[str, str | list[str]]:
    """Extract metadata from a genome package zip."""
    res: dict[str, str | list[str]] = {}
    try:
        with zipfile.ZipFile(genome_zip) as zf:
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
        _log.error(f"Bad zip file: {genome_zip}")
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
