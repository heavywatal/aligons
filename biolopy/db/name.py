from collections.abc import Iterable


def shorten(name: str):
    """Oryza_sativa -> osat"""
    split = name.lower().split("_")
    return split[0][0] + split[1][:3]


def filter_by_shortname(species: Iterable[str], queries: Iterable[str]):
    map = {k: v for v in species if (k := shorten(v)) in queries}
    for key in queries:
        yield map[key]


def main(argv: list[str] | None = None):
    import sys

    if argv is None:
        argv = sys.argv[1:]
    for x in argv:
        print(shorten(x))


if __name__ == "__main__":
    main()
