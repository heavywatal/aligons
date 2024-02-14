import logging
import os
import re
from collections.abc import Iterable
from pathlib import Path

from aligons.extern import htslib, jellyfish, kent
from aligons.util import cli, dl, fs, gff, subp

from . import _rsrc, jgi, mask

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    parser = cli.ArgumentParser()
    parser.add_argument("url")
    args = parser.parse_args(argv or None)
    dl_mirror_db(args.url)


def fetch_and_bgzip(
    entry: _rsrc.DataSet, prefix: Path
) -> tuple[cli.FuturePath, cli.FuturePath]:
    url_prefix = entry["url_prefix"]
    species = entry["species"]
    annotation = entry["annotation"]
    sequences = entry["sequences"]
    softmasked = re.search(r"dna_sm|softmasked", sequences[0])
    dna_sm = "dna_sm" if softmasked else "dna"
    stem = f"{species}_{ver}" if (ver := entry.get("version", None)) else species
    out_fa = prefix / f"fasta/{species}/{stem}.{dna_sm}.toplevel.fa.gz"
    out_gff = prefix / f"gff3/{species}/{stem}.gff3.gz"
    responses_fa = [dl_mirror_db(url_prefix + s) for s in sequences]
    ft_fa = cli.thread_submit(index_compress_concat, responses_fa, out_fa)
    response_gff = dl_mirror_db(url_prefix + annotation)
    ft_gff = cli.thread_submit(index_compress, response_gff, out_gff)
    return ft_fa, ft_gff


def dl_mirror_db(url: str) -> dl.Response:
    if "jgi.doe.gov" in url:
        return jgi.session.mirror(url, _rsrc.db_root())
    return dl.mirror(url, _rsrc.db_root())


def process_genome(futures: tuple[cli.FuturePath, cli.FuturePath]) -> cli.FuturePath:
    ft_fa, ft_gff = futures
    fa = ft_fa.result()
    ft_fasize = index_fasta([fa]) if ".dna_sm." in fa.name else split_mask_index(fa)
    return index_gff3([ft_gff.result()], ft_fasize.result())


def split_mask_index(toplevel_fa_gz: Path) -> cli.FuturePath:
    future_chromosomes = _split_toplevel_fa_work(toplevel_fa_gz)
    future_masked = [mask.submit(f) for f in future_chromosomes]
    links = [_symlink_masked(f) for f in future_masked]
    return index_fasta(links)


def _symlink_masked(ft: cli.FuturePath) -> Path:
    masked = ft.result()
    link = masked.parent.parent / masked.name
    return fs.symlink(masked, link, relative=True)


def index_compress(response: dl.Response, outfile: Path) -> Path:
    htslib.try_index(bgzip_or_symlink(response, outfile))
    return outfile


def index_compress_concat(responses: Iterable[dl.Response], outfile: Path) -> Path:
    htslib.try_index(concat_bgzip(responses, outfile))
    return outfile


def index_fasta(paths: list[Path]) -> cli.FuturePath:
    """Create bgzipped and indexed genome.fa."""
    if len(paths) == 1:
        if "chromosome" in paths[0].name:
            _log.warning(f"splitting chromosome? {paths[0]}")
        paths = [f.result() for f in _split_toplevel_fa(paths[0])]
    return cli.thread_submit(_create_genome_bgzip, paths)


def index_gff3(paths: list[Path], fasize: Path | None = None) -> cli.FuturePath:
    """Create bgzipped and indexed genome.gff3."""
    if len(paths) == 1:
        if "chromosome" in paths[0].name:
            _log.warning(f"splitting chromosome? {paths[0]}")
        assert fasize is not None, paths
        paths = htslib.split_gff3(paths[0], fasize)
    return cli.thread_submit(_create_genome_bgzip, paths)


def _create_genome_bgzip(files: list[Path]) -> Path:
    """Combine chromosome files and bgzip it."""
    files = fs.sorted_naturally(files)
    _log.debug(str(files))
    if cli.dry_run and not files:
        return Path(os.devnull)
    name = files[0].name
    ext = files[0].with_suffix("").suffix
    (outname, count) = re.subn(rf"\.chromosome\..+{ext}", rf".genome{ext}", name)
    assert count == 1, name
    outfile = files[0].with_name(outname)
    htslib.concat_bgzip(files, outfile)
    htslib.try_index(outfile)
    if outfile.name.endswith(".fa.gz"):
        return kent.faSize(outfile)
    return outfile


def _split_toplevel_fa_work(fa_gz: Path) -> list[cli.FuturePath]:
    fmt = "{stem}.{seqid}.fa"
    outdir = fa_gz.with_name("_work")
    return htslib.split_fa_gz(fa_gz, fmt, (r"toplevel", "chromosome"), outdir)


def _split_toplevel_fa(fa_gz: Path) -> list[cli.FuturePath]:
    if "toplevel" not in fa_gz.name:
        _log.warning(f"'toplevel' not in {fa_gz.name}")
    fmt = "{stem}.{seqid}.fa.gz"
    return htslib.split_fa_gz(fa_gz, fmt, (r"toplevel", "chromosome"))


def softmask(species: str) -> cli.FuturePath:
    masked = jellyfish.run(species)
    return index_fasta(masked)


def bgzip_or_symlink(infile: Path | dl.Response, outfile: Path) -> Path:
    """Decompress/compress/symlink depending on file names."""
    if isinstance(infile, dl.Response):
        infile = infile.path
    if fs.is_outdated(outfile, infile) and not cli.dry_run:
        outfile.parent.mkdir(0o755, parents=True, exist_ok=True)
        if htslib.to_be_bgzipped(outfile.name):
            fs.expect_suffix(outfile, ".gz")
            with subp.popen_zcat(infile) as zcat:
                if any(s in (".gff", ".gff3") for s in outfile.suffixes):
                    assert zcat.stdout is not None
                    if infile.name.startswith("PGSC_DM_V403"):
                        content = clean_seqid_pgsc(zcat.stdout)
                    else:
                        content, _ = zcat.communicate()
                    content = gff.sort(content)
                    htslib.bgzip(content, outfile)
                else:
                    htslib.bgzip(zcat.stdout, outfile)
        elif outfile.suffix == infile.suffix:
            fs.symlink(infile, outfile, relative=True)
        elif outfile.suffix == ".gz":
            _log.debug(f"rare case: gzip {infile = } to {outfile = }")
            with subp.popen_zcat(infile) as zcat:
                subp.gzip(zcat.stdout, outfile)
        else:
            msg = f"unexpected formats: {infile = }, {outfile = }"
            raise ValueError(msg)
    _log.info(f"{outfile}")
    return outfile


def concat_bgzip(responses: Iterable[dl.Response], outfile: Path) -> Path:
    return htslib.concat_bgzip([res.path for res in responses], outfile)


def clean_seqid_pgsc(stdin: subp.FILE) -> bytes:
    return subp.popen_sd("^ST4.03ch", "chr", stdin=stdin).communicate()[0]


if __name__ == "__main__":
    main()
