from __future__ import annotations

import logging
import re
from collections.abc import Callable, Iterable, Iterator
from typing import NamedTuple, TypeAlias

from .. import cli

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--name", action="store_true")
    parser.add_argument("-s", "--short", action="store_true")
    parser.add_argument("-g", "--graph", action="count")
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
        elif args.graph:
            newick = full_trees[args.clade]
            root = parse_newick(newick)
            if args.graph > 1:
                gen = render_nodes(root)
            else:
                gen = render_tips(root)
            for branch, label in gen:
                print(f"{branch} {label}")
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


def remove_inner(newick: str):
    return re.sub(r"\)[^(),;]+", ")", newick)


def make_trees():
    ehrhartoideae = "(oryza_sativa,leersia_perrieri)ehrhartoideae"
    pooideae = "(brachypodium_distachyon,(aegilops_tauschii,hordeum_vulgare))pooideae"
    andropogoneae = "(sorghum_bicolor,zea_mays)andropogoneae"
    paniceae = "(setaria_italica,panicum_hallii_fil2)paniceae"
    bep = f"({ehrhartoideae},{pooideae})bep"
    pacmad = f"({andropogoneae},{paniceae})pacmad"
    poaceae = f"({bep},{pacmad})poaceae"
    monocot = f"(({poaceae},musa_acuminata),dioscorea_rotundata)monocot"  # noqa: F841
    # pyright: reportUnusedVariable=false
    return {k: v + ";" for k, v in locals().items()}


full_trees = make_trees()
trees = {k: remove_inner(v) for k, v in full_trees.items()}


class Node(NamedTuple):
    name: str
    children: list[Node] = []
    distance: float | None = None


class Bricks(NamedTuple):
    first_branch: str = "├─"
    middle_branch: str = "├─"
    last_branch: str = "└─"
    straight: str = "│ "
    blank: str = "  "


StrGen: TypeAlias = Iterator[str]
GraphGen: TypeAlias = Iterable[tuple[str, str]]
NORMAL_BRICKS = Bricks()
COMPACT_BRICKS = Bricks("┬─")


def render_nodes(node: Node, columns: list[StrGen] = []) -> GraphGen:
    for line in node.name.splitlines():
        prefix = "".join([next(gen) for gen in columns])
        yield (prefix, line)
    for x in _iter_children(node, columns, render_nodes, NORMAL_BRICKS):
        yield x


def render_tips(node: Node, columns: list[StrGen] = []) -> GraphGen:
    for x in _iter_children(node, columns, render_tips, COMPACT_BRICKS):
        yield x
    if not node.children:
        prefix = "".join([next(gen) for gen in columns])
        yield (prefix, node.name)


def _iter_children(
    node: Node,
    columns: list[StrGen],
    render_fun: Callable[[Node, list[StrGen]], GraphGen],
    bricks: Bricks,
):
    for child in (children := node.children):
        if child.name == children[-1].name:
            gen = _column_generator(bricks.last_branch, bricks.blank)
        elif child.name == children[0].name:
            gen = _column_generator(bricks.first_branch, bricks.straight)
        else:
            gen = _column_generator(bricks.middle_branch, bricks.straight)
        columns.append(gen)
        for x in render_fun(child, columns):
            yield x
        columns.pop()


def _column_generator(top: str, bottom: str) -> StrGen:
    yield top
    while True:
        yield bottom


def parse_newick(newick: str):
    nodes = {}
    newick_old = ""
    while newick != newick_old:
        newick_old = newick
        newick, nodes = _extract_tip_clade(newick, nodes)
    root = nodes.popitem()[1]
    assert not nodes
    return root


def _extract_tip_clade(tree: str, nodes: dict[str, Node]):
    for mobj in re.finditer(r"\(([^(]+?)\)([^(),;]+)?", tree):
        children: list[Node] = []
        for label in mobj.group(1).split(","):
            name, distance = _parse_node_label(label)
            children.append(nodes.pop(name, Node(name, [], distance)))
        name, distance = _parse_node_label(mobj.group(2) or "")
        name = name or "+".join([x.name for x in children])
        nodes[name] = Node(name, children, distance)
        tree = tree.replace(mobj.group(0), name)
    return tree, nodes


def _parse_node_label(label: str):
    if ":" in label:
        name, distance = label.split(":")
        return name.strip(), float(distance.strip())
    return label.strip(), None


def rectangulate(renderer: GraphGen):
    lines = list(renderer)
    widest = max(lines, key=lambda p: len(p[0]) + len(p[1]))
    max_width = len(widest[0]) + len(widest[1])
    for prefix, label in lines:
        yield (prefix.ljust(max_width - len(label), "─"), label)


def elongate(renderer: GraphGen):
    lines = list(renderer)
    deepest = max(lines, key=lambda p: len(p[0]))
    max_depth = len(deepest[0])
    for prefix, label in lines:
        yield (prefix.ljust(max_depth, "─"), label)


if __name__ == "__main__":
    main()
