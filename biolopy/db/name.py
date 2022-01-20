def shorten(name: str):
    """Oryza_sativa -> osat"""
    split = name.lower().split('_')
    return split[0][0] + split[1][:3]


def filter(species: list[str], short_names: list[str]):
    return [x for x in species if shorten(x) in short_names]


def main(argv: list[str] | None = None):
    if argv is None:
        import sys
        argv = sys.argv[1:]
    for x in argv:
        print(shorten(x))


if __name__ == "__main__":
    main()
