"""Multiple genome alignment

src: ./pairwise/{target}/{query}/{chromosome}/sing.maf
lnk: ./multiple/{target}/{clade}/{chromosome}/{target}.{query}.sing.maf
dst: ./multiple/{target}/{clade}/{chromosome}/multiz.maf

https://github.com/multiz/multiz
"""
import concurrent.futures as confu
import itertools
import logging
import os
import shutil
from pathlib import Path

from ..db import phylo
from ..util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = cli.logging_argparser()
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("-n", "--dry-run", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("indir", type=Path)  # pairwise/oryza_sativa
    parser.add_argument("clade", choices=phylo.newicks.keys())
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    cli.dry_run = args.dry_run
    if args.config:
        cli.read_config(args.config)
    if args.clean:
        clean(Path("multiple") / args.indir.name)
        return
    run(args.indir, args.clade, args.jobs)


def run(indir: Path, clade: str, jobs: int):
    target = indir.name
    outdir = Path("multiple") / target / clade
    prepare(indir, outdir)
    chromodirs = outdir.glob("chromosome.*")
    multiz_opts = cli.config["multiz"]
    with confu.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [executor.submit(multiz, p, multiz_opts) for p in chromodirs]
    for future in confu.as_completed(futures):
        if (multiz_maf := future.result()).exists():
            print(multiz_maf)
    return outdir


def multiz(path: Path, options: subp.Optdict = {}):
    sing_mafs = list(path.glob("*.sing.maf"))
    clade = path.parent.name
    tmpdir = path / "_tmp"
    outfile = path / "multiz.maf"
    roasted = roast(sing_mafs, clade, tmpdir.name, outfile.name, options)
    script = "set -eu\n" + roasted.stdout
    is_to_run = fs.is_outdated(outfile, sing_mafs)
    if is_to_run and not cli.dry_run:
        tmpdir.mkdir(0o755, exist_ok=True)
        with open(path / "roasted.sh", "wt") as fout:
            fout.write(script)
    try:
        comp = subp.run_if(
            is_to_run,
            script,
            shell=True,
            cwd=path,
            stdout=subp.PIPE,
            stderr=subp.STDOUT,
        )
    except subp.CalledProcessError as perr:
        for line in perr.stdout.strip().splitlines():
            _log.error(f"{path.name}:{line}")
        _log.error(outfile)
        raise perr
    if is_to_run and not cli.dry_run:
        (tmpdir / outfile.name).replace(outfile)
        try:
            tmpdir.rmdir()
        except OSError as err:
            _log.warning(str(err))
    if out := comp.stdout.strip():
        _log.info(out)
    return outfile


def roast(
    sing_mafs: list[Path],
    clade: str,
    tmpdir: str,
    outfile: str,
    options: subp.Optdict = {},
):
    """Generate shell script to execute multiz"""
    radius = options.get("R", 30)
    min_width = options.get("M", 1)
    ref_label = sing_mafs[0].name.split(".", 1)[0]
    tree = phylo.newicks[clade]
    tree = phylo.shorten_names(tree).replace(",", " ").rstrip(";")
    args = (
        f"roast - R={radius} M={min_width} T={tmpdir} E={ref_label} '{tree}' "
        + " ".join([x.name for x in sing_mafs])
        + f" {tmpdir}/{outfile}"
    )
    return subp.run(args, stdout=subp.PIPE, text=True)


def prepare(indir: Path, outdir: Path):
    target = indir.name
    clade = outdir.name
    tree = phylo.newicks[clade]
    species = phylo.extract_names(tree)
    assert target in species
    if not cli.dry_run:
        outdir.mkdir(0o755, parents=True, exist_ok=True)
    for query in species:
        if query == target:
            continue
        querypath = indir / query
        if not querypath.exists():
            _log.warning(f"not found {querypath}")
        dstname = f"{phylo.shorten(target)}.{phylo.shorten(query)}.sing.maf"
        for chrdir in querypath.glob("chromosome.*"):
            src = chrdir / "sing.maf"
            dstdir = outdir / chrdir.name
            if not cli.dry_run:
                dstdir.mkdir(0o755, exist_ok=True)
            dst = dstdir / dstname
            relsrc = os.path.relpath(src, dstdir)
            _log.info(f"{dst}@\n -> {relsrc}")
            if not dst.exists() and not cli.dry_run:
                dst.symlink_to(relsrc)
    return outdir


def clean(path: Path):
    it = itertools.chain(
        path.rglob("multiz.maf"),
        path.rglob("roasted.sh"),
        path.rglob("_tmp"),
    )
    for entry in it:
        print(entry)
        if not cli.dry_run:
            rm_rf(entry)


def rm_rf(path: Path):
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


if __name__ == "__main__":
    main()
