from biolopy.db import name


def test_shorten():
    assert name.shorten("Oryza_sativa") == "osat"
