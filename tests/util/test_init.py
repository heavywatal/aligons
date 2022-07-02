from biolopy.util import update_nested


def test_update_nested():
    x = {"both": 1, "x": 1}
    y = {"both": 2, "y": 2}
    update_nested(x, y)
    assert x == {"both": 2, "x": 1, "y": 2}
    assert y == {"both": 2, "y": 2}
