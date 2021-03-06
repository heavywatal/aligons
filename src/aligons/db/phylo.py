from __future__ import annotations

import logging
import re
from collections.abc import Callable, Iterable, Iterator
from typing import NamedTuple, TypeAlias

from ..util import cli
from .ensemblgenomes import make_newicks, shorten

_log = logging.getLogger(__name__)


def main(argv: list[str] = []):
    parser = cli.logging_argparser()
    parser.add_argument("-n", "--name", action="store_true")
    parser.add_argument("-s", "--short", action="store_true")
    parser.add_argument("-g", "--graph", action="count")
    parser.add_argument("clade", nargs="?")
    args = parser.parse_args(argv or None)
    cli.logging_config(args.loglevel)
    global newicks
    if args.short:
        newicks = {k: shorten_names(v) for k, v in newicks.items()}
    if args.clade:
        tree = newicks[args.clade]
        if args.name:
            print(" ".join(extract_names(tree)))
        elif args.graph:
            newick = newicks_with_inner[args.clade]
            print_graph(newick, args.graph)
        else:
            print(tree)
        return
    for key, tree in newicks.items():
        print(f"{key}: {tree}")
    return


def sorted_by_len_newicks(clades: list[str], reverse: bool = False):
    return sorted(clades, key=lambda x: len(newicks[x]), reverse=reverse)


def extract_tip_names(newick: str):
    return extract_names(remove_inner(newick))


def extract_names(newick: str):
    names = (x.split(":")[0] for x in re.findall(r"[^\s(),;]+", newick))
    return list(filter(None, names))


def extract_lengths(newick: str):
    return [float(x.lstrip()) for x in re.findall(r"(?<=:)\s*[\d.]+", newick)]


def shorten_names(newick: str):
    return re.sub(r"[^\s(),:;]+_[^\s(),:;]+", lambda m: shorten(m.group(0)), newick)


def remove_lengths(newick: str):
    return re.sub(r":\s*[\d.]+", "", newick)


def remove_inner(newick: str):
    return re.sub(r"\)[^(),;]+", ")", newick)


def remove_inner_names(newick: str):
    return re.sub(r"\)[^(),:;]+", ")", newick)


def remove_whitespace(x: str):
    return "".join(x.split())


newicks_with_inner = make_newicks()
newicks = {k: remove_inner(v) for k, v in newicks_with_inner.items()}


def print_graph(newick: str, graph: int = 0):
    root = parse_newick(newick)
    if graph >= 4:
        gen = rectangulate(render_tips(root))
    elif graph == 3:
        gen = elongate(render_tips(root))
    elif graph == 2:
        gen = render_tips(root)
    else:
        gen = render_nodes(root)
    for branch, label in gen:
        print(f"{branch} {label}")


class Node(NamedTuple):
    name: str
    children: list[Node] = []
    distance: float | None = None


class Bricks(NamedTuple):
    first_branch: str = "??????"
    middle_branch: str = "??????"
    last_branch: str = "??????"
    straight: str = "??? "
    blank: str = "  "


StrGen: TypeAlias = Iterator[str]
GraphGen: TypeAlias = Iterable[tuple[str, str]]
NORMAL_BRICKS = Bricks()
COMPACT_BRICKS = Bricks("??????")


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
        yield (prefix.ljust(max_width - len(label), "???"), label)


def elongate(renderer: GraphGen):
    lines = list(renderer)
    deepest = max(lines, key=lambda p: len(p[0]))
    max_depth = len(deepest[0])
    for prefix, label in lines:
        yield (prefix.ljust(max_depth, "???"), label)


if __name__ == "__main__":
    main()
