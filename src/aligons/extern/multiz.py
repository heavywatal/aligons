"""Multiple genome alignment.

src: ./pairwise/{target}/{query}/{chromosome}/sing.maf
lnk: ./multiple/{target}/{clade}/{chromosome}/{target}.{query}.sing.maf
dst: ./multiple/{target}/{clade}/{chromosome}/multiz.maf

https://github.com/multiz/multiz
"""
import itertools
import logging
import os
import shutil
from collections.abc import Sequence
from pathlib import Path

from aligons.db import phylo
from aligons.util import cli, config, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    nodes_all = phylo.extract_names(phylo.get_tree())
    parser = cli.ArgumentParser()
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("indir", type=Path)  # pairwise/oryza_sativa
    parser.add_argument("query", nargs="+", choices=nodes_all)
    args = parser.parse_args(argv or None)
    if args.clean:
        clean(Path("multiple") / args.indir.name)
        return
    run(args.indir, args.query)


def run(indir: Path, query: Sequence[str]) -> Path:
    target = indir.name
    tree = phylo.get_subtree(query)
    if len(query) > 1:
        dirname = "-".join(phylo.shorten(x) for x in query)
        if target not in query:
            msg = f"{target=} not in {query=}"
            raise ValueError(msg)
    else:
        dirname = query[0]
        query = phylo.extract_names(tree)
    outdir = Path("multiple") / target / dirname
    prepare(indir, outdir, query)
    chromodirs = outdir.glob("chromosome.*")
    pool = cli.ThreadPool()
    cli.wait_raise(pool.submit(multiz, p, tree) for p in chromodirs)
    return outdir


def multiz(path: Path, tree: str) -> Path:
    sing_mafs = list(path.glob("*.sing.maf"))
    tmpdir = path / "_tmp"
    outfile = path / "multiz.maf"
    roasted = roast(sing_mafs, tmpdir.name, outfile.name, tree)
    script = "set -eu\n" + roasted.stdout
    is_to_run = fs.is_outdated(outfile, sing_mafs)
    if is_to_run and not cli.dry_run:
        tmpdir.mkdir(0o755, exist_ok=True)
        with (path / "roasted.sh").open("w") as fout:
            fout.write(script)
    try:
        comp = subp.run(
            script,
            if_=is_to_run,
            shell=True,  # noqa: S604
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
    if outfile.exists():
        print(outfile)
    return outfile


def roast(
    sing_mafs: list[Path], tmpdir: str, outfile: str, tree: str
) -> subp.subprocess.CompletedProcess[str]:
    """Generate shell script to execute multiz."""
    options = config["multiz"]
    radius = options.get("R", 30)
    min_width = options.get("M", 1)
    ref_label = sing_mafs[0].name.split(".", 1)[0]
    tree = phylo.shorten_names(tree).replace(",", " ").rstrip(";")
    args = (
        f"roast - R={radius} M={min_width} T={tmpdir} E={ref_label} '{tree}' "
        + " ".join([x.name for x in sing_mafs])
        + f" {tmpdir}/{outfile}"
    )
    return subp.run(args, stdout=subp.PIPE, text=True)


def prepare(indir: Path, outdir: Path, queries: Sequence[str]) -> Path:
    target = indir.name
    if target not in queries:
        msg = f"{target=} not in {queries=}"
        raise ValueError(msg)
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


def clean(path: Path) -> None:
    it = itertools.chain(
        path.rglob("multiz.maf"),
        path.rglob("roasted.sh"),
        path.rglob("_tmp"),
    )
    for entry in it:
        print(entry)
        if not cli.dry_run:
            rm_rf(entry)


def rm_rf(path: Path) -> None:
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


if __name__ == "__main__":
    main()
