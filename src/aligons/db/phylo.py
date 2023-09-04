from __future__ import annotations

import functools
import logging
import re
from collections.abc import Callable, Iterable, Iterator, Sequence
from pathlib import Path
from typing import NamedTuple, TypeAlias

from aligons.util import cli, config

_log = logging.getLogger(__name__)


def main(argv: list[str] | None = None):
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
        print("\n".join(extract_names(tree)))
    elif args.graph:
        print_graph(tree, args.graph)
    else:
        print(tree)


def sorted_by_len_newicks(clades: list[str], *, reverse: bool = False):
    return sorted(clades, key=lambda x: len(get_subtree([x])), reverse=reverse)


def extract_tip_names(newick: str):
    return extract_names(remove_inner(newick))


def extract_inner_names(newick: str) -> list[str]:
    return re.findall(r"(?<=\))[^\s(),;:]+", newick)


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


def select(newick: str, queries: Sequence[str]):
    if len(queries) == 1:
        newick = select_clade(newick, queries[0])
    elif len(queries) > 1:
        newick = select_tips(newick, queries)
    return newick


def select_clade(newick: str, clade: str):
    return newickize(parse_newick(newick, clade))


def select_tips(newick: str, tips: Sequence[str]):
    def repl(mobj: re.Match[str]):
        if (s := mobj.group(0)) in tips:
            return s
        return ""

    filtered = re.sub(r"[\w_]+", repl, newick)
    while re.search(r"\(,\)", filtered):
        filtered = re.sub(r"\(,\)", "", filtered)
    _log.debug(filtered)
    root = parse_newick(filtered)
    return newickize(root)


def get_subtree(queries: Sequence[str], fun: Callable[[str], str] = lambda x: x):
    tree = fun(get_tree())
    return remove_inner(select(tree, queries))


@functools.cache
def get_tree() -> str:
    if not (tree := config["db"].get("tree", "")):
        tree = make_newick()
    elif "," not in tree:
        with Path(tree).expanduser().open("rt") as fin:
            tree = fin.read()
    return tree


@functools.cache
def make_newick() -> str:
    ehrhartoideae = "(oryza_sativa,leersia_perrieri)ehrhartoideae"
    pooideae = "(brachypodium_distachyon,(aegilops_tauschii,hordeum_vulgare))pooideae"
    andropogoneae = "(sorghum_bicolor,zea_mays)andropogoneae"
    paniceae = "(setaria_italica,panicum_hallii_fil2)paniceae"
    bep = f"({ehrhartoideae},{pooideae})bep"
    pacmad = f"({andropogoneae},{paniceae})pacmad"
    poaceae = f"({bep},{pacmad})poaceae"
    poales = f"({poaceae},ananas_comosus)poales"
    monocot = f"(({poales},musa_acuminata),dioscorea_rotundata,)monocot"

    _lycopersicon = "(solanum_lycopersicum,solanum_pennellii)"
    solanum = f"(({_lycopersicon},solanum_tuberosum),solanum_melongena)solanum"
    solanaceae = f"(({solanum},capsicum_annuum),nicotiana_attenuata)solanaceae"
    _convolvulaceae = "ipomoea_triloba"
    solanales = f"({solanaceae},{_convolvulaceae})solanales"
    _lamiales = "(olea_europaea_sylvestris,sesamum_indicum)lamiales"
    lamiids = f"(({solanales},coffea_canephora),{_lamiales})lamiids"
    _asteraceae = "(helianthus_annuus,lactuca_sativa)asteraceae"
    _companulids = f"({_asteraceae},daucus_carota)"
    _core_asterids = f"({lamiids},{_companulids})"
    asterids = f"({_core_asterids},actinidia_chinensis)asterids"
    eudicots = f"({asterids},arabidopsis_thaliana)eudicots"

    return f"({eudicots},{monocot})angiospermae"


def list_species(clade: str = "angiospermae") -> list[str]:
    return extract_names(get_subtree([clade]))


def expand_shortnames(shortnames: list[str]):
    return _filter_by_shortname(list_species(), shortnames)


def _filter_by_shortname(species: Iterable[str], queries: Iterable[str]):
    return (x for x in species if shorten(x) in queries)


def shorten(name: str):
    """Oryza_sativa -> osat."""
    if name.lower() == "olea_europaea_sylvestris":
        return "oesy"
    split = name.lower().split("_")
    return split[0][0] + split[1][:3]


def print_graph(newick: str, graph: int = 0):
    root = parse_newick(newick)
    if graph >= 4:  # noqa: PLR2004
        gen = rectangulate(render_tips(root, []))
    elif graph == 3:  # noqa: PLR2004
        gen = elongate(render_tips(root, []))
    elif graph == 2:  # noqa: PLR2004
        gen = render_tips(root, [])
    else:
        gen = render_nodes(root, [])
    for branch, label in gen:
        print(f"{branch} {label}")


class Node(NamedTuple):
    name: str
    children: list[Node] | None = None
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


def newickize(node: Node) -> str:
    return _newickize(node) + ";"


def _newickize(node: Node) -> str:
    ret = ""
    if node.children:
        ret += "("
        ret += ",".join([_newickize(child) for child in node.children])
        ret += ")"
    if "+" not in node.name:
        ret += node.name
    if node.distance:
        ret += f":{node.distance}"
    return ret


def render_nodes(node: Node, columns: list[StrGen]) -> GraphGen:
    for line in node.name.splitlines():
        prefix = "".join([next(gen) for gen in columns])
        yield (prefix, line)
    yield from _iter_children(node, columns, render_nodes, NORMAL_BRICKS)


def render_tips(node: Node, columns: list[StrGen]) -> GraphGen:
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
    if inner:
        assert inner in newick, f"{inner} not found in {newick}"
    nodes = {}
    newick_old = ""
    while newick != newick_old:
        newick_old = newick
        newick, nodes = _extract_tip_clade(newick, nodes)
        if inner and (clade := nodes.get(inner, None)):
            return clade
    root = nodes.popitem()[1]
    assert not nodes, nodes
    return root


def _extract_tip_clade(tree: str, nodes: dict[str, Node]):
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
