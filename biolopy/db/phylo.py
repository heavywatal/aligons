import argparse
import logging
import re

from .. import cli

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--name", action="store_true")
    parser.add_argument("-s", "--short", action="store_true")
    parser.add_argument("clade", nargs="?")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    global trees
    if args.short:
        trees = {k: shorten_labels(v) for k, v in trees.items()}
    if args.clade:
        tree = trees[args.clade]
        if args.name:
            print(" ".join(extract_labels(tree)))
        else:
            print(tree)
        return
    for key, tree in trees.items():
        print(f"{key}: {tree}")
    return


def extract_labels(tree: str):
    labels = (x.split(":")[0] for x in re.findall(r"[^ (),;]+", tree))
    return list(filter(None, labels))


def extract_lengths(tree: str):
    return [float(x) for x in re.findall(r"(?<=:)[\d.]+", tree)]


def remove_lengths(tree: str):
    return re.sub(r":[\d.]+", "", tree)


def shorten_labels(tree: str):
    return re.sub(r"[^ (),:;]+_[^ (),:;]+", lambda m: shorten(m.group(0)), tree)


def shorten(name: str):
    """Oryza_sativa -> osat"""
    split = name.lower().split("_")
    return split[0][0] + split[1][:3]


def make_trees():
    ehrhartoideae = "(oryza_sativa,leersia_perrieri)"
    pooideae = "(brachypodium_distachyon,(aegilops_tauschii,hordeum_vulgare))"
    andropogoneae = "(sorghum_bicolor,zea_mays)"
    paniceae = "(setaria_italica,panicum_hallii_fil2)"
    bep = f"({ehrhartoideae},{pooideae})"
    pacmad = f"({andropogoneae},{paniceae})"
    poaceae = f"({bep},{pacmad})"
    monocot = f"(({poaceae},musa_acuminata),dioscorea_rotundata)"  # noqa: F841
    # pyright: reportUnusedVariable=false
    return {k: v + ";" for k, v in locals().items()}


trees = make_trees()


if __name__ == "__main__":
    main()
