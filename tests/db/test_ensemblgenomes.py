from aligons.db import ensemblgenomes


def test_match_gff3_name():
    v = ensemblgenomes.version()
    name = f"Oryza_sativa.IRGSP-1.0.{v}.gff3.gz"
    md = ensemblgenomes.match_gff3_name(name)
    assert md["species"] == "Oryza_sativa"
    assert md["asm"] == "IRGSP-1.0"
    assert md["v"] == f"{v}"
    assert md["type"] is None
    assert md["seqid"] is None
    name = f"Oryza_sativa.IRGSP-1.0.{v}.chr.gff3.gz"
    md = ensemblgenomes.match_gff3_name(name)
    assert md["species"] == "Oryza_sativa"
    assert md["asm"] == "IRGSP-1.0"
    assert md["type"] == "chr"
    assert md["seqid"] is None
    fmt = "{species}.{asm}.{v}.genome.gff3.gz"
    assert fmt.format(**md) == f"Oryza_sativa.IRGSP-1.0.{v}.genome.gff3.gz"
    name = f"Oryza_sativa.IRGSP-1.0.{v}.chromosome.1.gff3.gz"
    md = ensemblgenomes.match_gff3_name(name)
    assert md["species"] == "Oryza_sativa"
    assert md["asm"] == "IRGSP-1.0"
    assert md["type"] == "chromosome"
    assert md["seqid"] == "1"
    name = f"H_vulgare.MorexV3_pseudomolecules_assembly.{v}.primary_assembly.1H.gff3.gz"
    md = ensemblgenomes.match_gff3_name(name)
    assert md["type"] == "primary_assembly"
    assert md["seqid"] == "1H"
    name = f"species.{v}.0.{v}.chr.gff3.gz"
    md = ensemblgenomes.match_gff3_name(name)
    assert md["species"] == "species"
    assert md["asm"] == f"{v}.0"
    assert md["type"] == "chr"
    assert md["seqid"] is None


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
