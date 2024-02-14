from aligons.db import ensemblgenomes


def test_replace_label_gff3():
    v = ensemblgenomes.version()
    exp = f"stem.{v}.genome.gff3.gz"
    assert ensemblgenomes.replace_label_gff3(f"stem.{v}.chr.gff3.gz") == exp
    assert ensemblgenomes.replace_label_gff3(f"stem.{v}.gff3.gz") == exp
    assert ensemblgenomes.replace_label_gff3(f"stem.{v}.chromosome.{v}.gff3.gz") == exp
    exp = f"stem.{v}.0.{v}.genome.gff3.gz"
    assert ensemblgenomes.replace_label_gff3(f"stem.{v}.0.{v}.chr.gff3.gz") == exp
    assert ensemblgenomes.replace_label_gff3(f"stem.{v}.0.{v}.gff3.gz") == exp
    assert (
        ensemblgenomes.replace_label_gff3(f"stem.{v}.0.{v}.chromosome.{v}.gff3.gz")
        == exp
    )


def test_replace_label_fa():
    for dna in ("dna", "dna_rm", "dna_sm"):
        exp = f"stem.{dna}.genome.fa.gz"
        assert ensemblgenomes.replace_label_fa(f"stem.{dna}.chromosome.1.fa.gz") == exp
        assert ensemblgenomes.replace_label_fa(f"stem.{dna}.toplevel.fa.gz") == exp
