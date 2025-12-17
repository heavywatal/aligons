"""RepeatMasker finds interspersed repeats and low complexity sequences.

- <https://www.repeatmasker.org/>
- <https://github.com/Dfam-consortium/TETools>
"""

import io
import logging
import re
import tempfile
from pathlib import Path

import polars as pl

from aligons.util import cli, fs, subp

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    """CLI for manual execution and testing."""
    parser = cli.ArgumentParser()
    parser.add_argument("--test", action="store_true")
    parser.add_argument("-S", "--species")
    parser.add_argument("infile", type=Path, nargs="*")
    args = parser.parse_args(argv or None)
    if args.test:
        if args.species:
            test_species(args.species)
        else:
            test_famdb_angiosperms()
        return
    for infile in args.infile:
        cli.thread_submit(repeatmasker, infile, args.species)


def repeatmasker(infile: Path, species: str = "", *, soft: bool = True) -> Path:
    """Be careful of messy and dirty output from RepeatMasker.

    - Returns 0 even if aborted, e.g., "cannot read file".
    - Prints verbose log to stdout, not stderr.
    - `*.out.gff` lacks some information stored in `*.out`.

    :param infile: Input FASTA file: `{basename}.fa`
    :param species: Species or higher taxonomic term for repeat library.
    :param soft: If `True`, mask low-complexity regions in lowercase.
    :returns: Output GFF file: `{basename}.fa.out.gff`
    """
    fs.expect_suffix(infile, ".gz", negate=True)
    outfile = infile.with_suffix(infile.suffix + ".out.gff")
    parallel: int = 2
    args: subp.Args = ["RepeatMasker", "-e", "rmblast", "-gff"]
    if soft:
        args.append("-xsmall")
    if parallel > 1:
        args.extend(["-pa", f"{parallel}"])  # RMBlast uses 4 cores per "-pa"
    if species:
        args.extend(["-species", species])
    args.append(infile)
    subp.run(args, if_=fs.is_outdated(outfile, infile))
    if not cli.dry_run:
        with outfile.open("r") as fin:  # raise FileNotFoundError if failed
            line1 = fin.readline()
            if not line1.startswith("##gff-version 3"):
                msg = f"RepeatMasker produced incomplete {outfile}:\n{line1}"
                raise ValueError(msg)
    return fs.print_if_exists(outfile)


def read_out(infile: Path) -> pl.LazyFrame:
    """Read RepeatMasker `.out` file into a LazyFrame."""
    fs.expect_suffix(infile, ".out")
    with infile.open("rb") as fin:
        content = re.sub(rb" *\n *", rb"\n", fin.read())
        content = re.sub(rb" +", rb"\t", content)
    return pl.read_csv(
        io.BytesIO(content),
        has_header=False,
        separator="\t",
        skip_rows=3,
        new_columns=[
            "sw_score",
            "perc_div",
            "perc_del",
            "perc_ins",
            "seqid",
            "begin",
            "end",
            "left",
            "strand",
            "repeat",
            "class",
            "rep_begin",
            "rep_end",
            "rep_left",
            "id",
        ],
    ).lazy()


def test_species(species: str) -> bool:
    """Test if species is recognized by RepeatMasker.

    Existence in NCBI taxonomy does not suffice.
    <https://www.ncbi.nlm.nih.gov/taxonomy>
    """
    with tempfile.TemporaryDirectory(prefix="RepeatMasker") as tmp, fs.chdir(tmp):
        _log.info(f"{Path().absolute()}")
        fasta = Path("test.fa")
        if not fasta.exists():
            with fasta.open("wt") as fout:
                fout.write(">name\nAAACCCGGGTTT\n")
        args: subp.Args = ["RepeatMasker", "-species", species, fasta]
        args.extend(["-qq", "-nolow", "-norna", "-no_is", "-noisy", "-nopost"])
        p = subp.run(args, stdout=subp.PIPE, stderr=subp.PIPE, check=False)
        _log.info(p.stdout.decode())
        stderr = p.stderr.decode()
        _log.info(stderr)
        if not p.returncode:
            return True
        mobj = re.search(r"Species .+ is not known to \w+", stderr)
        assert mobj, "Unexpected error from RepeatMasker"
        _log.warning(mobj.group(0))
        return False


def test_famdb_angiosperms() -> bool:
    """Test. TODO."""
    assert not famdb_families("angiosperms", descendants=True)
    fa = famdb_families("oryza_sativa", ancestors=True)
    for mobj in re.finditer(r">.+?\n", fa):
        seqname = mobj.group(0)
        assert "@root" in seqname, seqname
    _log.info("No angiosperms has specific repeat families registered.")
    return True


def famdb_families(
    term: str,
    *,
    ancestors: bool = False,
    descendants: bool = False,
    fmt: str = "fasta_name",
) -> str:
    """Run `famdb.py families` to get repeat families for a taxonomic term.

    :param term: Taxonomic term like species name.
    :param ancestors: If `True`, include families from ancestor taxa.
    :param descendants: If `True`, include families from descendant taxa.
    :param fmt: Output format:
        summary, hmm, hmm_species, fasta_name, fasta_acc, embl, embl_meta, embl_seq
    :returns: Families in specified format.
    """
    args: subp.Args = ["famdb.py", "families", "-f", fmt]
    if ancestors:
        args.append("-a")
    if descendants:
        args.append("-d")
    args.append(term)
    p = subp.run(args, stdout=subp.PIPE)
    return p.stdout.decode()


if __name__ == "__main__":
    main()
