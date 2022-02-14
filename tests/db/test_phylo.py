from biolopy.db import phylo


def test_shorten():
    assert phylo.shorten("Oryza_sativa") == "osat"
