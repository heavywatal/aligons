"""Multiple genome alignment.

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
from collections.abc import Sequence
from pathlib import Path

from aligons.db import phylo
from aligons.util import ConfDict, cli, config, empty_options, fs, read_config, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
    nodes_all = phylo.extract_names(phylo.newicks_with_inner["angiospermae"])
    parser = cli.ArgumentParser()
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("-j", "--jobs", type=int, default=os.cpu_count())
    parser.add_argument("-c", "--config", type=Path)
    parser.add_argument("indir", type=Path)  # pairwise/oryza_sativa
    parser.add_argument("query", choices=nodes_all)
    args = parser.parse_args(argv or None)
    if args.config:
        read_config(args.config)
    if args.clean:
        clean(Path("multiple") / args.indir.name)
        return
    run(args.indir, args.query, args.jobs)


def run(indir: Path, query: Sequence[str], jobs: int):
    target = indir.name
    assert target in query
    tree = phylo.get_newick(query)
    if len(query) > 1:
        dirname = "-".join(phylo.shorten(x) for x in query)
    else:
        dirname = query[0]
        query = phylo.extract_names(tree)
    outdir = Path("multiple") / target / dirname
    prepare(indir, outdir, query)
    chromodirs = outdir.glob("chromosome.*")
    multiz_opts = config["multiz"]
    multiz_opts["tree"] = tree
    with confu.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [executor.submit(multiz, p, multiz_opts) for p in chromodirs]
    for future in confu.as_completed(futures):
        if (multiz_maf := future.result()).exists():
            print(multiz_maf)
    return outdir


def multiz(path: Path, options: ConfDict = empty_options):
    sing_mafs = list(path.glob("*.sing.maf"))
    tmpdir = path / "_tmp"
    outfile = path / "multiz.maf"
    roasted = roast(sing_mafs, tmpdir.name, outfile.name, options)
    script = "set -eu\n" + roasted.stdout
    is_to_run = fs.is_outdated(outfile, sing_mafs)
    if is_to_run and not cli.dry_run:
        tmpdir.mkdir(0o755, exist_ok=True)
        with (path / "roasted.sh").open("w") as fout:
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
        raise
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
    tmpdir: str,
    outfile: str,
    options: ConfDict = empty_options,
):
    """Generate shell script to execute multiz."""
    tree = options["tree"]
    tree = phylo.shorten_names(tree).replace(",", " ").rstrip(";")
    radius = options.get("R", 30)
    min_width = options.get("M", 1)
    ref_label = sing_mafs[0].name.split(".", 1)[0]
    args = (
        f"roast - R={radius} M={min_width} T={tmpdir} E={ref_label} '{tree}' "
        + " ".join([x.name for x in sing_mafs])
        + f" {tmpdir}/{outfile}"
    )
    return subp.run(args, stdout=subp.PIPE, text=True)


def prepare(indir: Path, outdir: Path, queries: Sequence[str]):
    target = indir.name
    assert target in queries
    if not cli.dry_run:
        outdir.mkdir(0o755, parents=True, exist_ok=True)
    for query in queries:
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
