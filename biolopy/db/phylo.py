import argparse
import logging
import re

from .. import cli
from . import name

_log = logging.getLogger(__name__)

clades = {
    "ehrhartoideae": "(oryza_sativa,leersia_perrieri)",
    "pooideae": "(brachypodium_distachyon,(aegilops_tauschii,hordeum_vulgare))",
    "andropogoneae": "(sorghum_bicolor,zea_mays)",
    "paniceae": "(setaria_italica,panicum_hallii_fil2)",
}
clades["bep"] = "({ehrhartoideae},{pooideae})".format(**clades)
clades["pacmad"] = "({andropogoneae},{paniceae})".format(**clades)
clades["poaceae"] = "({bep},{pacmad})".format(**clades)
clades["monocot"] = "(({poaceae},musa_acuminata),dioscorea_rotundata)".format(**clades)


def main(argv: list[str] = []):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
    parser.add_argument("-n", "--name", action="store_true")
    parser.add_argument("-s", "--short", action="store_true")
    parser.add_argument("clade", nargs="?")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    global clades
    if args.short:
        clades = {k: shorten(v) for k, v in clades.items()}
    if args.clade:
        tree = clades[args.clade]
        if args.name:
            print(" ".join(extract_labels(tree)))
        else:
            print(tree)
        return
    for key, tree in clades.items():
        print(f"{key}: {tree}")
    return


def shorten(tree: str):
    return re.sub(
        r"([^ (),:;]+)(:[\d.]+)?",
        lambda m: name.shorten(m.group(1)) + (m.group(2) or ""),
        tree,
    )


def extract_labels(tree: str):
    return [x.split(":")[0] for x in re.findall(r"[^ (),;]+", tree)]


if __name__ == "__main__":
    main()
