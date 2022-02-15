from biolopy.db import phylo

tree = "((oryza_sativa:0.1,leersia_perrieri:0.2):0.3,panicum_hallii_fil2:0.4):0.5;"


def test_shorten():
    assert phylo.shorten("Oryza_sativa") == "osat"


def test_remove_lengths():
    exp = "((oryza_sativa,leersia_perrieri),panicum_hallii_fil2);"
    assert phylo.remove_lengths(tree) == exp
    assert phylo.remove_lengths(phylo.shorten_labels(tree)) == "((osat,lper),phal);"


def test_shorten_labels():
    assert phylo.shorten_labels(tree) == "((osat:0.1,lper:0.2):0.3,phal:0.4):0.5;"
    assert phylo.shorten_labels(phylo.remove_lengths(tree)) == "((osat,lper),phal);"


def test_extract_labels():
    labels = ["oryza_sativa", "leersia_perrieri", "panicum_hallii_fil2"]
    short_labels = ["osat", "lper", "phal"]
    short_tree = phylo.shorten_labels(tree)
    assert phylo.extract_labels(tree) == labels
    assert phylo.extract_labels(phylo.remove_lengths(tree)) == labels
    assert phylo.extract_labels(short_tree) == short_labels
    assert phylo.extract_labels(phylo.remove_lengths(short_tree)) == short_labels


def test_extract_lengths():
    assert phylo.extract_lengths(tree) == [0.1, 0.2, 0.3, 0.4, 0.5]
