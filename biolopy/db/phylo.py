from __future__ import annotations

import argparse
import logging
import re
from collections.abc import Callable, Generator, Iterable
from typing import NamedTuple, TypeAlias

from .. import cli

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = argparse.ArgumentParser(parents=[cli.logging_argparser()])
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
            newick = trees[args.clade]
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


class Node(NamedTuple):
    name: str
    children: list[Node] = []
    distance: float | None = None


class Branches(NamedTuple):
    fork: str = "├─"
    straight: str = "│ "
    corner: str = "└─"
    blank: str = "  "


StrGen: TypeAlias = Generator[str, None, None]
GraphGen: TypeAlias = Generator[tuple[str, str], None, None]
BRANCHES = Branches()
COMPACT_BRANCHES = Branches("┬─")


def render_nodes(node: Node, columns: list[StrGen] = []) -> GraphGen:
    for line in node.name.splitlines():
        prefix = "".join([next(gen) for gen in columns])
        yield (prefix, line)
    for x in _iter_children(node, columns, render_nodes, BRANCHES):
        yield x


def render_tips(node: Node, columns: list[StrGen] = []) -> GraphGen:
    for x in _iter_children(node, columns, render_tips, COMPACT_BRANCHES):
        yield x
    if not node.children:
        prefix = "".join([next(gen) for gen in columns])
        yield (prefix, node.name)


def _iter_children(
    node: Node,
    columns: list[StrGen],
    render_fun: Callable[[Node, list[StrGen]], GraphGen],
    branches: Branches,
):
    for child in (children := node.children):
        if child != children[-1]:
            gen = _column_generator(branches.fork, branches.straight)
        else:
            gen = _column_generator(branches.corner, branches.blank)
        columns.append(gen)
        for x in render_fun(child, columns):
            yield x
        columns.pop()


def _column_generator(first: str, rest: str) -> StrGen:
    yield first
    while True:
        yield rest


def parse_newick(newick: str):
    nodes = {}
    newick_old = ""
    while newick != newick_old:
        newick_old = newick
        newick, nodes = _extract_tip_trios(newick, nodes)
    root = nodes.popitem()[1]
    assert not nodes
    return root


def _extract_tip_trios(tree: str, nodes: dict[str, Node]):
    for mobj in re.finditer(r"\(([^(]+?),([^(]+?)\)([^(),;]+)?", tree):
        name1 = mobj.group(1)
        name2 = mobj.group(2)
        node1 = nodes.pop(name1, _make_node(name1))
        node2 = nodes.pop(name2, _make_node(name2))
        name3 = f"{node1.name}+{node2.name}" + (mobj.group(3) or "")
        nodes[name3] = _make_node(name3, [node1, node2])
        tree = tree.replace(mobj.group(0), name3)
    return tree, nodes


def _make_node(label: str, children: list[Node] = []) -> Node:
    if ":" in label:
        (name, dist) = label.split(":")
        return Node(name.strip(), children, float(dist.strip()))
    return Node(label.strip(), children)


def rectangulate(renderer: Iterable[tuple[str, str]]):
    lines = list(renderer)
    widest = max(lines, key=lambda p: len(p[0]) + len(p[1]))
    max_width = len(widest[0]) + len(widest[1])
    for prefix, label in lines:
        yield (prefix.ljust(max_width - len(label), "─"), label)


def elongate(renderer: Iterable[tuple[str, str]]):
    lines = list(renderer)
    deepest = max(lines, key=lambda p: len(p[0]))
    max_depth = len(deepest[0])
    for prefix, label in lines:
        yield (prefix.ljust(max_depth, "─"), label)


if __name__ == "__main__":
    main()
