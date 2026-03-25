#!/usr/bin/env python3
"""Prune bac120 tree to keep only phyla present in a filtered CSV file."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_TREE = Path("output/bac120.tree_stripped")
DEFAULT_TAXONOMY = Path("bac120_taxonomy.tsv")
DEFAULT_OUTPUT = Path("output/bac120.tree_stripped_pruned_by_filtered_phyla")


@dataclass
class Node:
    name: str = ""
    length: str | None = None
    children: list["Node"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children


def parse_newick(newick: str) -> Node:
    """Parse Newick string into a tree of Node objects."""

    def skip_ws(i: int) -> int:
        while i < len(newick) and newick[i].isspace():
            i += 1
        return i

    def read_token(i: int, stop_chars: set[str]) -> tuple[str, int]:
        i = skip_ws(i)
        start = i
        while i < len(newick) and newick[i] not in stop_chars:
            i += 1
        return newick[start:i].strip(), i

    def parse_subtree(i: int) -> tuple[Node, int]:
        i = skip_ws(i)
        if i >= len(newick):
            raise ValueError("Unexpected end of Newick while parsing subtree.")

        if newick[i] == "(":
            i += 1
            children: list[Node] = []
            while True:
                child, i = parse_subtree(i)
                children.append(child)
                i = skip_ws(i)
                if i >= len(newick):
                    raise ValueError("Unexpected end of Newick inside internal node.")
                if newick[i] == ",":
                    i += 1
                    continue
                if newick[i] == ")":
                    i += 1
                    break
                raise ValueError(f"Unexpected character '{newick[i]}' in internal node.")

            name, i = read_token(i, {":", ",", ")", ";"})
            length = None
            i = skip_ws(i)
            if i < len(newick) and newick[i] == ":":
                length, i = read_token(i + 1, {",", ")", ";"})
            return Node(name=name, length=length, children=children), i

        name, i = read_token(i, {":", ",", ")", ";"})
        if not name:
            raise ValueError("Encountered leaf node without a name.")
        length = None
        i = skip_ws(i)
        if i < len(newick) and newick[i] == ":":
            length, i = read_token(i + 1, {",", ")", ";"})
        return Node(name=name, length=length), i

    root, i = parse_subtree(0)
    i = skip_ws(i)
    if i >= len(newick) or newick[i] != ";":
        raise ValueError("Newick must end with a semicolon ';'.")
    return root


def _sum_lengths(length_a: str | None, length_b: str | None) -> str | None:
    if length_a is None:
        return length_b
    if length_b is None:
        return length_a
    try:
        return format(float(length_a) + float(length_b), ".10g")
    except ValueError:
        return length_a


def prune_tree(node: Node, keep_leaves: set[str]) -> Node | None:
    """Prune tree to only leaves in keep_leaves and collapse unary internal nodes."""
    if node.is_leaf:
        return node if node.name in keep_leaves else None

    pruned_children: list[Node] = []
    for child in node.children:
        pruned_child = prune_tree(child, keep_leaves)
        if pruned_child is not None:
            pruned_children.append(pruned_child)

    if not pruned_children:
        return None

    node.children = pruned_children

    if len(node.children) == 1:
        only_child = node.children[0]
        only_child.length = _sum_lengths(only_child.length, node.length)
        return only_child

    return node


def node_to_newick(node: Node) -> str:
    """Serialize Node tree back to Newick string."""
    if node.children:
        subtree = f"({','.join(node_to_newick(child) for child in node.children)}){node.name}"
    else:
        subtree = node.name

    if node.length is not None:
        subtree += f":{node.length}"
    return subtree


def extract_target_phyla(csv_path: Path) -> set[str]:
    with csv_path.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        if not reader.fieldnames or "phylum" not in reader.fieldnames:
            raise ValueError("Filtered CSV must contain a 'phylum' column.")
        return {row["phylum"].strip() for row in reader if row.get("phylum", "").strip()}


def extract_keep_genomes(taxonomy_path: Path, target_phyla: set[str]) -> set[str]:
    keep_genomes: set[str] = set()

    with taxonomy_path.open("r", encoding="utf-8") as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue

            parts = line.split(None, 1)
            if len(parts) != 2:
                continue

            genome, taxonomy = parts
            phylum = None
            for rank in taxonomy.split(";"):
                if rank.startswith("p__"):
                    phylum = rank[3:]
                    break

            if phylum and phylum in target_phyla:
                keep_genomes.add(genome)

    return keep_genomes


def count_leaves(node: Node) -> int:
    if node.is_leaf:
        return 1
    return sum(count_leaves(child) for child in node.children)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prune bac120.tree_stripped to only genomes whose phylum appears in "
            "a filtered CSV file."
        )
    )
    parser.add_argument(
        "filtered_csv",
        type=Path,
        help="Path to filtered CSV (must include a 'phylum' column).",
    )
    parser.add_argument(
        "--tree",
        type=Path,
        default=DEFAULT_TREE,
        help=f"Input Newick tree path (default: {DEFAULT_TREE})",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        default=DEFAULT_TAXONOMY,
        help=f"GTDB taxonomy TSV path (default: {DEFAULT_TAXONOMY})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output pruned Newick path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def main() -> None:
    sys.setrecursionlimit(1_000_000)
    args = parse_args()

    target_phyla = extract_target_phyla(args.filtered_csv)
    if not target_phyla:
        raise ValueError("No phyla found in filtered CSV; nothing to keep.")

    keep_genomes = extract_keep_genomes(args.taxonomy, target_phyla)
    if not keep_genomes:
        raise ValueError("No genomes from taxonomy matched the filtered phyla.")

    newick = args.tree.read_text(encoding="utf-8").strip()
    root = parse_newick(newick)
    original_leaf_count = count_leaves(root)

    pruned_root = prune_tree(root, keep_genomes)
    if pruned_root is None:
        raise ValueError("Pruning removed all leaves; check your inputs.")

    pruned_leaf_count = count_leaves(pruned_root)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(node_to_newick(pruned_root) + ";\n", encoding="utf-8")

    print(f"Target phyla: {len(target_phyla)}")
    print(f"Genomes kept (from taxonomy): {len(keep_genomes)}")
    print(f"Leaves before pruning: {original_leaf_count}")
    print(f"Leaves after pruning: {pruned_leaf_count}")
    print(f"Wrote pruned tree to {args.output}")


if __name__ == "__main__":
    main()
