"""Utilities for phylogenetic trees."""

from __future__ import annotations

import functools
import logging
import re
from collections.abc import Callable, Iterable, Iterator, Sequence
from pathlib import Path
from typing import NamedTuple

from aligons.util import cli, config, resources_data

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None) -> None:
    """CLI for manual execution and testing."""
    parser = cli.ArgumentParser()
    parser.add_argument("-N", "--name", action="store_true")
    parser.add_argument("-s", "--short", action="store_true")
    parser.add_argument("-i", "--inner", action="store_true")
    parser.add_argument("-g", "--graph", action="count")
    parser.add_argument("query", nargs="*")
    args = parser.parse_args(argv or None)
    tree = select(get_tree(), args.query)
    if args.short:
        tree = shorten_names(tree)
    if not args.inner:
        tree = remove_inner(tree)
    if args.name:
        _log.info("\n".join(extract_names(tree)))
    elif args.graph:
        print_graph(tree, args.graph)
    else:
        _log.info(tree)


def sorted_by_len_newicks(clades: list[str], *, reverse: bool = False) -> list[str]:
    """Sort clade names by the length of their subtrees.

    :param clades: Clade names.
    :param reverse: Descending order if `True`.
    :returns: Sorted clade names.
    """
    return sorted(clades, key=lambda x: len(get_subtree([x])), reverse=reverse)


def extract_tip_names(newick: str) -> list[str]:
    """Extract tip (leaf) node names from a Newick tree.

    :param newick: Species tree in Newick format.
    :returns: Tip names.
    """
    return extract_names(remove_inner(newick))


def extract_inner_names(newick: str) -> list[str]:
    """Extract inner (ancestral) node names from a Newick tree.

    :param newick: Species tree in Newick format.
    :returns: Inner node names.
    """
    return re.findall(r"(?<=\))[^\s(),;:]+", newick)


def extract_names(newick: str) -> list[str]:
    """Extract node names including inner nodes.

    :param newick: Species tree in Newick format.
    :returns: Node names including both tip and inner nodes.
    """
    names = (x.split(":")[0] for x in re.findall(r"[^\s(),;]+", newick))
    return list(filter(None, names))


def extract_lengths(newick: str) -> list[float]:
    """Extract branch lengths from a Newick tree.

    :param newick: Species tree in Newick format.
    :returns: Branch lengths.
    """
    return [float(x.lstrip()) for x in re.findall(r"(?<=:)\s*[\d.]+", newick)]


def shorten_names(newick: str) -> str:
    """Shorten node names in the Newick tree.

    :param newick: Species tree in Newick format.
    :returns: Newick tree with shortened names.
    """
    return re.sub(r"[^\s(),:;]+_[^\s(),:;]+", lambda m: shorten(m.group(0)), newick)


def remove_lengths(newick: str) -> str:
    """Remove branch lengths from a Newick tree."""
    return re.sub(r":\s*[\d.]+", "", newick)


def remove_inner(newick: str) -> str:
    """Remove inner node names and branch lengths from a Newick tree."""
    return re.sub(r"\)[^(),;]+", ")", newick)


def remove_inner_names(newick: str) -> str:
    """Remove inner node names from a Newick tree."""
    return re.sub(r"\)[^(),:;]+", ")", newick)


def remove_whitespace(x: str) -> str:
    """Remove all whitespace characters from a string."""
    return "".join(x.split())


def select(newick: str, queries: Sequence[str]) -> str:
    """Extract subtree containing specified tips or a clade.

    :param newick: Original species tree in Newick format.
    :param queries: Clade name or tip names to extract. No filtering if empty.
    :returns: Subtree in Newick format.
    """
    if len(queries) == 1:
        newick = select_clade(newick, queries[0])
    elif len(queries) > 1:
        newick = select_tips(newick, queries)
    return newick


def select_clade(newick: str, clade: str) -> str:
    """Extract subtree containing the specified clade.

    :param newick: Original species tree in Newick format.
    :param clade: Clade name to extract.
    :returns: Subtree in Newick format.
    """
    return to_newick(parse_newick(newick, clade))


def select_tips(newick: str, tips: Sequence[str]) -> str:
    """Extract subtree containing the specified tips.

    :param newick: Original species tree in Newick format.
    :param tips: Tip names to extract.
    :returns: Subtree in Newick format.
    """

    def repl(mobj: re.Match[str]) -> str:
        if (s := mobj.group(0)) in tips:
            return s
        return ""

    filtered = re.sub(r"[\w_]+", repl, newick)
    while re.search(r"\(,\)", filtered):
        filtered = re.sub(r"\(,\)", "", filtered)
    _log.debug(filtered)
    root = parse_newick(filtered)
    return to_newick(root)


def get_subtree(queries: Sequence[str], fun: Callable[[str], str] = lambda x: x) -> str:
    """Extract subtree containing specified tips or a clade.

    :param queries: Tip names or clade name to extract.
    :param fun: A function to modify the tree before extraction.
    :returns: Subtree in Newick format.
    """
    tree = fun(get_tree())
    return remove_inner(select(tree, queries))


@functools.cache
def get_tree() -> str:
    """Get species tree from `config.db.tree` or built-in file.

    :returns: Species tree in Newick format.
    """
    tree = config["db"].get("tree", "") or read_builtin_newick()
    if "," not in tree:
        with Path(tree).expanduser().open("rt") as fin:
            tree = fin.read()
    tree = re.sub(r"\s", "", tree)
    return tree.removesuffix(";")


@functools.cache
def read_builtin_newick() -> str:
    """Read built-in Newick tree from resources."""
    return resources_data("angiospermae.nhx").read_text()


def list_species(clade: str = "") -> list[str]:
    """List species names in the given clade.

    :param clade: Clade name. No filtering if empty.
    """
    tree = get_subtree([clade] if clade else [])
    names = extract_names(tree)
    exclude = config["db"]["exclude"]
    if exclude:
        names = list(filter(lambda x: x not in exclude, names))
    return names


def lengthen(species: str) -> str:
    """Expand a shortened species name to the full name."""
    try:
        return next(_expand_short_names([species]))
    except StopIteration:
        _log.warning(f"cannot expand {species = }")
        return species


def _expand_short_names(short_names: list[str]) -> Iterator[str]:
    return _filter_by_shortname(list_species(), short_names)


def _filter_by_shortname(
    species: Iterable[str], queries: Iterable[str]
) -> Iterator[str]:
    return (x for x in species if shorten(x) in queries)


def shorten(name: str) -> str:
    """Extract four letters from a species name.

    e.g., Oryza_sativa -> osat.

    :param name: Species name.
    :returns: Shortened name with four letters.
    """
    if name.lower() == "olea_europaea_sylvestris":
        return "oesy"
    if name.lower() == "oryza_sativa_mh63":
        return "mh63"
    if name.lower() == "oryza_sativa_zs97":
        return "zs97"
    split = name.lower().split("_")
    return split[0][0] + split[1][:3]


def print_graph(newick: str, graph: int = 0) -> None:
    """Print tree structure in a graphical format.

    :param newick: Species tree in Newick format.
    :param graph: Graph style level:
        1. simple tree structure
        2. with tip names
        3. elongated branches
        4. rectangular branches
    """
    root = parse_newick(newick)
    if graph >= 4:  # noqa: PLR2004
        gen = rectangular(render_tips(root, []))
    elif graph == 3:  # noqa: PLR2004
        gen = elongate(render_tips(root, []))
    elif graph == 2:  # noqa: PLR2004
        gen = render_tips(root, [])
    else:
        gen = render_nodes(root, [])
    for branch, label in gen:
        _log.info(f"{branch} {label}")


class Node(NamedTuple):
    """Recursive node to represent a tree."""

    name: str
    children: list[Node] | None = None
    distance: float | None = None


class Bricks(NamedTuple):
    """Characters for drawing tree branches."""

    first_branch: str = "├─"
    middle_branch: str = "├─"
    last_branch: str = "└─"
    straight: str = "│ "
    blank: str = "  "


type StrGen = Iterator[str]
type GraphGen = Iterable[tuple[str, str]]
NORMAL_BRICKS = Bricks()
COMPACT_BRICKS = Bricks("┬─")


def to_newick(node: Node) -> str:
    """Convert a Node tree to a Newick format string.

    :param node: Root node.
    :returns: The tree in Newick format.
    """
    return _to_newick(node) + ";"


def _to_newick(node: Node) -> str:
    ret = ""
    if node.children:
        ret += "("
        ret += ",".join([_to_newick(child) for child in node.children])
        ret += ")"
    if "+" not in node.name:
        ret += node.name
    if node.distance:
        ret += f":{node.distance}"
    return ret


def render_nodes(node: Node, columns: list[StrGen]) -> GraphGen:
    """Generate lines for rendering tree including inner nodes.

    :param node: Current node.
    :param columns: List of column generators for branch drawing.
    :returns: An iterator of (prefix, node label) tuples.
    """
    for line in node.name.splitlines():
        prefix = "".join([next(gen) for gen in columns])
        yield (prefix, line)
    yield from _iter_children(node, columns, render_nodes, NORMAL_BRICKS)


def render_tips(node: Node, columns: list[StrGen]) -> GraphGen:
    """Generate lines for rendering tree without inner nodes.

    :param node: Current node.
    :param columns: List of column generators for branch drawing.
    :returns: An iterator of (prefix, node label) tuples.
    """
    yield from _iter_children(node, columns, render_tips, COMPACT_BRICKS)
    if not node.children:
        prefix = "".join([next(gen) for gen in columns])
        yield (prefix, node.name)


def _iter_children(
    node: Node,
    columns: list[StrGen],
    render_fun: Callable[[Node, list[StrGen]], GraphGen],
    bricks: Bricks,
) -> GraphGen:
    for child in (children := node.children or []):
        if child.name == children[-1].name:
            gen = _column_generator(bricks.last_branch, bricks.blank)
        elif child.name == children[0].name:
            gen = _column_generator(bricks.first_branch, bricks.straight)
        else:
            gen = _column_generator(bricks.middle_branch, bricks.straight)
        columns.append(gen)
        yield from render_fun(child, columns)
        columns.pop()


def _column_generator(top: str, bottom: str) -> StrGen:
    yield top
    while True:
        yield bottom


def parse_newick(newick: str, inner: str = "") -> Node:
    """Parse Newick string into a tree of Nodes.

    :param newick: Species tree in Newick format.
    :param inner: Name of an ancestral node of the clade to extract.
    :returns: The root Node of the tree.
    """
    if inner not in newick:
        msg = f"{inner} not in {newick}"
        raise ValueError(msg)
    nodes = {}
    newick_old = ""
    while newick != newick_old:
        newick_old = newick
        newick, nodes = _extract_tip_clade(newick, nodes)
        if inner and (clade := nodes.get(inner)):
            return clade
    root = nodes.popitem()[1]
    if nodes:
        msg = f"{newick} may have multiple roots: {nodes}"
        raise ValueError(msg)
    return root


def _extract_tip_clade(
    tree: str, nodes: dict[str, Node]
) -> tuple[str, dict[str, Node]]:
    for mobj in re.finditer(r"\(([^(]+?)\)([^(),;]+)?", tree):
        children: list[Node] = []
        for label in mobj.group(1).split(","):
            if label:
                name, distance = _parse_node_label(label)
                children.append(nodes.pop(name, Node(name, [], distance)))
        if len(children) > 1:
            name, distance = _parse_node_label(mobj.group(2) or "")
            name = name or "+".join([x.name for x in children])
            nodes[name] = Node(name, children, distance)
        else:
            name = children[0].name
            nodes[name] = children[0]
        tree = tree.replace(mobj.group(0), name)
    return tree, nodes


def _parse_node_label(label: str) -> tuple[str, float | None]:
    if ":" in label:
        name, distance = label.split(":")
        return name.strip(), float(distance.strip())
    return label.strip(), None


def rectangular(renderer: GraphGen) -> GraphGen:
    """Generate lines with adjusted branch lengths for rectangular tree shape.

    :param renderer: An iterator of (prefix, node label) tuples.
    :returns: Modified iterator with adjusted prefixes.
    """
    lines = list(renderer)
    widest = max(lines, key=lambda p: len(p[0]) + len(p[1]))
    max_width = len(widest[0]) + len(widest[1])
    for prefix, label in lines:
        yield (prefix.ljust(max_width - len(label), "─"), label)


def elongate(renderer: GraphGen) -> Iterable[tuple[str, str]]:
    """Generate lines with elongated branches.

    :param renderer: An iterator of (prefix, node label) tuples.
    :returns: Modified iterator with elongated prefixes.
    """
    lines = list(renderer)
    deepest = max(lines, key=lambda p: len(p[0]))
    max_depth = len(deepest[0])
    for prefix, label in lines:
        yield (prefix.ljust(max_depth, "─"), label)


if __name__ == "__main__":
    main()
