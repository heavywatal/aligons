"""Download genomes from NCBI.

<https://www.ncbi.nlm.nih.gov/datasets/docs/>
"""

import hashlib
import json
import re
import tomllib
import zipfile
from typing import TYPE_CHECKING

import polars as pl

if TYPE_CHECKING:
    from pathlib import Path

from aligons.extern import htslib
from aligons.util import cli, config, dl, fs, logging, resources_data, subp

from . import _rsrc, api, tools

_log = logging.getLogger(__name__)
_HOST = "ftp.ncbi.nlm.nih.gov"


def main(argv: list[str] | None = None) -> None:
    """CLI for downloading and processing NCBI genomes."""
    parser = cli.ArgumentParser()
    parser.add_argument("-C", "--check", action="store_true")
    parser.add_argument("-D", "--download", action="store_true")
    parser.add_argument("-M", "--deploy", action="store_true")
    parser.add_argument("-l", "--list", action="store_true")
    parser.add_argument("--geo", action="store_true")
    parser.add_argument("accession", nargs="*")
    args = parser.parse_args(argv or None)
    if args.list:
        _ls(verbose=args.verbose > 1)
    elif args.geo:
        for file in args.accession:
            cache = _download_geo_supplementary(file)
            _deploy_cluster_bed(cache)
    elif args.download or args.check or args.deploy:
        fts = [cli.thread_submit(_download_genome, acc) for acc in args.accession]
        if args.check:
            # not so parallel due to GIL
            fts_check = [
                cli.thread_submit(_check_zip, ft) for ft in cli.as_completed(fts)
            ]
            fts = fts_check
        if args.deploy:
            fts_idx: list[cli.Future[Path]] = []
            for ft in cli.as_completed(fts):
                fts_idx.extend(_index_genome(ft))
            fts: list[cli.Future[Path]] = []
            for ft in cli.as_completed(fts_idx):
                if (res := ft.result()).name.endswith(".fa.gz"):
                    fts.extend(tools.genome_to_twobits(res))
        cli.wait_raise(fts)
    else:
        _view_downloaded_genomes()


def _view_downloaded_genomes() -> None:
    for genome_zip in _prefix_mirror().glob("*.zip"):
        info = _get_info(genome_zip)
        seq_report = info.pop("seq-report", [])
        _log.info(info)
        _log.info(_filter_sequences(seq_report))  # pyright: ignore[reportArgumentType]


def db_prefix() -> Path:
    """Directory of preprocessed NCBI genomes."""
    return api.prefix("ncbi")


def _prefix_mirror() -> Path:
    return _rsrc.db_root(_HOST) / "genomes" / "all"


def _download_geo_supplementary(filename: str) -> Path:
    """Download a GEO supplementary file.

    :param file: File name like `GSE131202_ZS97_H3K4me3.cluster.bed.gz`.
    :returns: Path to the downloaded file.
    """
    accession = filename.split("_", 1)[0]
    cache_dir = _rsrc.db_root(_HOST) / "geo" / accession
    if not cli.dry_run:
        cache_dir.mkdir(0o755, parents=True, exist_ok=True)
    url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={accession}&format=file&file={filename}"
    response = dl.fetch(url, cache_dir / filename)
    return response.path


def _deploy_cluster_bed(file: Path) -> Path:
    """Deploy a GEO cluster BED file to the database.

    :param file: Path to a downloaded file.
    :returns: Path to the deployed file.
    """
    stem = file.name.split(".", 1)[0]
    accession = stem.split("_", 1)[0]
    outdir = db_prefix() / "geo" / accession
    if not cli.dry_run:
        outdir.mkdir(0o755, parents=True, exist_ok=True)
    outfile = htslib.bgzip(_flatten_cluster(file), outdir / f"{stem}.bed.gz")
    htslib.tabix(outfile)
    return outfile


def _flatten_cluster(file: Path) -> bytes:
    """Flatten GSE131202_*.cluster.bed.

    :param file: Path to `GSE131202_*.cluster.bed.gz`.
    :returns: TSV bytes in three-column BED format.
    """
    cols = ["chr1", "start1", "end1", "chr2", "start2", "end2", "count", "p", "p_adj"]
    lf = pl.scan_csv(file, separator="\t", has_header=False, new_columns=cols)
    lf1 = lf.select(
        pl.col("chr1").alias("chr"),
        pl.col("start1").alias("start"),
        pl.col("end1").alias("end"),
    )
    lf2 = lf.select(
        pl.col("chr2").alias("chr"),
        pl.col("start2").alias("start"),
        pl.col("end2").alias("end"),
    )
    lf = pl.concat([lf1, lf2], how="vertical").unique().sort(["chr", "start", "end"])
    return lf.collect().write_csv(separator="\t", include_header=False).encode()


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
    if config["db"].get("ncbi", None).get("api-key", None):
        args.extend(["--api-key", config["db"]["ncbi"]["api-key"]])
    args.extend(["accession", accession])
    args.extend(["--filename", outfile])
    args.extend(["--include", "genome,gff3,seq-report"])
    subp.run(args, if_=fs.is_outdated(outfile))
    return outfile


def _index_genome(genome_zip: Path | cli.Future[Path]) -> list[cli.Future[Path]]:
    """Make indexed genome FASTA and GFF from a genome data package zip.

    > All genome sequences are softmasked using WindowMasker or RepeatMasker.

    <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/>

    :param genome_zip: NCBI genome package file: `{accession}.zip`.
    :returns: Future of the bgzipped files.
    """
    if isinstance(genome_zip, cli.Future):
        genome_zip = genome_zip.result()
    info = _get_info(genome_zip)
    _log.debug(info)
    species: str = info["species"]  # pyright: ignore[reportAssignmentType]
    species = species.replace(" ", "_").lower()
    version = info["version"]
    stem = f"{species}_{version}"
    name = info["sequences"][0]
    seq_report = info.get("seq-report", [])
    seq_names = _filter_sequences(seq_report)  # pyright: ignore[reportArgumentType]
    _log.debug(f"{seq_names = }")  # TODO: filter and rename sequences
    out_fa = db_prefix() / species / f"{stem}.dna_sm.genome.fa.gz"
    fts = [cli.thread_submit(tools.unzip_index_bgzip, genome_zip, name, out_fa)]
    if isinstance(gff := info.get("annotation"), str) and gff:
        _log.debug(f"{gff = }")
        out_gff = db_prefix() / species / f"{stem}.genome.gff3.gz"
        fts.append(cli.thread_submit(tools.unzip_index_bgzip, genome_zip, gff, out_gff))
    return fts


def _get_info(
    genome_zip: Path,
) -> dict[str, str | list[str] | list[dict[str, str | int]]]:
    """Extract metadata from a genome package zip."""
    res: dict[str, str | list[str] | list[dict[str, str | int]]] = {}
    try:
        with zipfile.ZipFile(genome_zip) as zf:
            for name in zf.namelist():
                _log.debug(name)
            sequences = [x for x in zf.namelist() if "_genomic.f" in x]
            acc, res["version"] = _get_asm_version(sequences[0])
            res["accession"] = acc
            res["sequences"] = sequences
            with zf.open("ncbi_dataset/data/assembly_data_report.jsonl") as fin:
                for line in fin:
                    data = json.loads(line)
                    organism = data["organism"]
                    _log.debug(f"{organism}")
                    res["species"] = organism["organismName"]
            seq_report = f"ncbi_dataset/data/{acc}/sequence_report.jsonl"
            if seq_report in zf.namelist():
                with zf.open(seq_report) as fin:
                    res["seq-report"] = [_read_seq_report(line) for line in fin]
            annotation = [x for x in zf.namelist() if ".gff" in x]
            if annotation:
                res["annotation"] = annotation[0]
    except zipfile.BadZipFile:
        _log.error(f"Bad zip file: {genome_zip}")
    if res.get("species", "") == "Eustoma russellianum":
        res["species"] = "Eustoma grandiflorum"  # NCBI and Wikipedia/en are wrong
    if res.get("species", "") == "Gentiana dahurica var. dahurica":
        res["species"] = "Gentiana dahurica"  # for simplicity
    return res


def _read_seq_report(line: bytes) -> dict[str, str | int]:
    """Parse a line of `sequence_report.jsonl`."""
    obj = json.loads(line)
    keys = (
        "chrName",
        "length",
        "refseqAccession",
        "genbankAccession",
        "role",
        "sequenceName",
    )
    return {k: v for k, v in obj.items() if k in keys}


def _filter_sequences(
    seq_report: list[dict[str, str | int]], min_length: int = 1000000
) -> dict[str, str]:
    """Filter sequences by role and length.

    :returns: A dictionary of {accession: name}.
    """
    seq_names: dict[str, str] = {}
    for rec in seq_report:
        if (
            rec.get("role", "") == "assembled-molecule"
            or int(rec["length"]) > min_length
        ):
            key, name = _get_sequence_name(rec)
            seq_names[key] = name
    return seq_names


def _get_sequence_name(rec: dict[str, str | int]) -> tuple[str, str]:
    key = rec.get("refseqAccession") or rec.get("genbankAccession")
    if not key or not isinstance(key, str):
        msg = f"Invalid sequence report: {rec}"
        raise ValueError(msg)
    name = str(rec.get("chrName")) or "Un"
    if name == "Un":
        name = ""
    return (key, name)


def _get_asm_version(name: str) -> tuple[str, str]:
    """Extract accession and assembly version from a file name inside a zip."""
    accession, name = name.removeprefix("ncbi_dataset/data/").split("/", 1)
    x = name.removeprefix(f"{accession}_")
    return accession, re.sub(r"_genomic\..+", "", x)


def _check_zip(file: Path | cli.Future[Path]) -> Path:
    """Check integrity of a zip file and its MD5 checksums."""
    if isinstance(file, cli.Future):
        file = file.result()
    try:
        with zipfile.ZipFile(file) as zf:
            if zf.testzip() is not None:
                _log.error(f"Corrupted zip file: {file}")
            _check_md5(zf)
    except zipfile.BadZipFile:
        _log.error(f"Bad zip file: {file}")
    return file


def _check_md5(zf: zipfile.ZipFile, sum_file: str = "md5sum.txt") -> None:
    """Check MD5 checksums inside a zip file."""
    with zf.open(sum_file) as sum_fin:
        for line in sum_fin:
            hex_exp, path = line.decode().strip().split(None, 1)
            with zf.open(path) as fin:
                md5_obs = hashlib.md5(fin.read())  # noqa: S324
                if md5_obs.hexdigest() != hex_exp:
                    _log.error(f"MD5 mismatch: {path}")


def _ls(*, verbose: bool = False) -> None:
    """List built-in datasets."""
    for ds in _builtin_datasets():
        if verbose:
            _log.info(f"{ds.get('accession', '')} {ds['species']}")
        else:
            _log.info(ds.get("accession", ""))


def ucsc_resources(species: str) -> list[str]:
    """List UCSC Genome Browser resources for a given species.

    :param species: Species name like `danio_rerio`.
    :returns: List of URLs.
    """
    for dataset in _builtin_datasets():
        if dataset["species"].replace(" ", "_").lower() == species:
            break
    else:
        msg = f"Species not found in built-in datasets: {species}"
        raise ValueError(msg)
    accession = dataset.get("accession", "")
    name = dataset["version"]
    slashed = accession.replace("_", "/")
    for pos in (7, 10):
        slashed = _insert_str(slashed, pos, "/")
    host = "hgdownload.soe.ucsc.edu"
    prefix = f"https://{host}/hubs/{slashed}/{accession}/bbi/{accession}_{name}."
    files = [
        "tandemDups.bb",
        "cpgIslandExt.bb",
        "rmsk.bb",
        "simpleRepeat.bb",
        "windowMasker.bb",
    ]
    return [prefix + filename for filename in files]


def _insert_str(string: str, pos: int, value: str) -> str:
    return string[:pos] + value + string[pos:]


def _builtin_datasets() -> list[_rsrc.DataSet]:
    """Return datasets defined in the bundled `ncbi.toml`."""
    with resources_data("ncbi.toml").open("rb") as fin:
        return tomllib.load(fin)["dataset"]


if __name__ == "__main__":
    main()
