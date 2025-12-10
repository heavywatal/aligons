"""Multiple genome alignment.

<https://github.com/multiz/multiz>
"""

import logging
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence

from aligons.db import phylo
from aligons.util import cli, config, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    nodes_all = phylo.extract_names(phylo.get_tree())
    parser = cli.ArgumentParser()
    parser.add_argument("indir", type=Path)  # pairwise/oryza_sativa
    parser.add_argument("query", nargs="+", choices=nodes_all)
    args = parser.parse_args(argv or None)
    run(args.indir, args.query)


def run(indir: Path, query: Sequence[str]) -> Path:
    """Prepare inputs and perform multiple alignment with multiz.

    - src: `./pairwise/{target}/{query}/chromosome.{seqid}/sing.maf`
    - lnk: `./multiple/{target}/{clade}/chromosome.{seqid}/{target}.{query}.sing.maf`
    - dst: `./multiple/{target}/{clade}/chromosome.{seqid}/multiz.maf`

    :param indir: Pairwise alignment directory with target species name:
        `./pairwise/{target}/`
    :param query: Query species names or a clade name.
    :returns: Output directory for the multiple alignment:
        `./multiple/{target}/{clade}/`
    """
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
    chromo_dirs = outdir.glob("chromosome.*")
    pool = cli.ThreadPool()
    cli.wait_raise(pool.submit(multiz, p, tree) for p in chromo_dirs)
    return outdir


def multiz(path: Path, tree: str) -> Path:
    """Perform multiple alignment with multiz.

    :param path: Directory for multiz input (pairwise `*.sing.maf`) and output:
        `./multiple/{target}/{clade}/chromosome.{seqid}/`
    :param tree: Species tree in Newick format.
    :returns: Output MAF file:
        `./multiple/{target}/{clade}/chromosome.{seqid}/multiz.maf`
    """
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
        comp = subp.run(  # noqa: S604
            script,
            if_=is_to_run,
            shell=True,
            cwd=path,
            stdout=subp.PIPE,
            stderr=subp.STDOUT,
        )
    except subp.CalledProcessError as p_err:
        for line in p_err.stdout.strip().splitlines():
            _log.error(f"{path.name}:{line}")
        _log.error(outfile)
        raise
    if is_to_run and not cli.dry_run:
        (tmpdir / outfile.name).replace(outfile)
        try:
            tmpdir.rmdir()
        except OSError as err:
            _log.warning(err)
    if out := comp.stdout.strip():
        _log.info(out)
    return fs.print_if_exists(outfile)


def roast(
    sing_mafs: list[Path], tmpdir: str, outfile: str, tree: str
) -> subp.CompletedProcess[str]:
    """Generate shell script to execute multiz.

    Never use bash because it cannot pass {tree} to roast.
    > roast.v3: tree specification contains too many '('<br>
    > roast.v3: improper character in tree specification: "

    :param sing_mafs: Pairwise `*.sing.maf` files.
    :param tmpdir: Temporary directory.
    :param outfile: Output MAF file: `multiz.maf`.
    :param tree: Species tree in Newick format.
    :returns: `CompletedProcess` with stdout containing shell script.
    """
    options = config["multiz"]
    radius = options.get("R", 30)
    min_width = options.get("M", 1)
    ref_label = sing_mafs[0].name.split(".", 1)[0]
    tree = phylo.shorten_names(tree).replace(",", " ").rstrip(";")
    opts = [f"R={radius}", f"M={min_width}", f"T={tmpdir}"]
    args = ["roast", "-", *opts, f"E={ref_label}", f"{tree}"]
    args.extend([x.name for x in sing_mafs])
    args.append(f"{tmpdir}/{outfile}")
    return subp.run(args, stdout=subp.PIPE, text=True)


def prepare(indir: Path, outdir: Path, queries: Sequence[str]) -> Path:
    """Create output directory with symlinks to pairwise alignment files.

    :param indir: Pairwise alignment directory with target species name.
    :param outdir: Directory for multiz inputs and outputs.
    :param queries: Query species names.
    :returns: The same path as `outdir`.
    """
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
        for chr_dir in querypath.glob("chromosome.*"):
            src = chr_dir / "sing.maf"
            link = outdir / chr_dir.name / dstname
            fs.symlink(src, link, relative=True)
    return outdir


if __name__ == "__main__":
    main()
