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


def test_match_fa_name():
    fa = "Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz"
    md = ensemblgenomes.match_fa_name(fa)
    assert md["species"] == "Oryza_sativa"
    assert md["asm"] == "IRGSP-1.0"
    assert md["dna"] == "dna"
    assert md["type"] == "toplevel"
    assert md["seqid"] is None
    fmt = "{species}.{asm}.{dna}.genome.fa.gz"
    assert fmt.format(**md) == "Oryza_sativa.IRGSP-1.0.dna.genome.fa.gz"
    fa = "Panicum_hallii_fil2.PHallii_v3.1.dna_sm.chromosome.1.fa.gz"
    md = ensemblgenomes.match_fa_name(fa)
    assert md["species"] == "Panicum_hallii_fil2"
    assert md["asm"] == "PHallii_v3.1"
    assert md["dna"] == "dna_sm"
    assert md["type"] == "chromosome"
    assert md["seqid"] == "1"
    fa = "H_vulgare.MorexV3_pseudomolecules_assembly.dna.primary_assembly.1H.fa.gz"
    md = ensemblgenomes.match_fa_name(fa)
    assert md["type"] == "primary_assembly"
    assert md["seqid"] == "1H"
