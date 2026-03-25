#!/usr/bin/env python3
"""Generate iTOL phylum color files from a summary CSV of CO2/C1 fractions."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DEFAULT_SUMMARY = Path("output/c_by_e_donor_phylum_formate_threshold_1000_summary.csv")
DEFAULT_TEMPLATE_STRIP = Path("output/itol_dataset_strip_phylum.txt")
DEFAULT_OUTPUT_STRIP = Path("output/itol_dataset_strip_phylum_co2_c1_ratio.txt")
DEFAULT_OUTPUT_TREE = Path("output/itol_tree_colours_co2_c1_ratio.txt")
DEFAULT_OUTPUT_STRIP_COLLAPSED = Path("output/itol_dataset_strip_phylum_co2_c1_ratio_collapsed_tree.txt")
DEFAULT_OUTPUT_TREE_COLLAPSED = Path("output/itol_tree_colours_co2_c1_ratio_collapsed_tree.txt")

# CO2 -> green, C1 -> purple
GREEN_RGB = (44, 162, 95)
PURPLE_RGB = (122, 59, 184)


def rgb_to_hex(rgb: tuple[int, int, int]) -> str:
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def interpolate_color(purple_weight: float) -> str:
    """Return a blended hex color from green (0) to purple (1)."""
    w = min(1.0, max(0.0, purple_weight))
    r = round(GREEN_RGB[0] * (1 - w) + PURPLE_RGB[0] * w)
    g = round(GREEN_RGB[1] * (1 - w) + PURPLE_RGB[1] * w)
    b = round(GREEN_RGB[2] * (1 - w) + PURPLE_RGB[2] * w)
    return rgb_to_hex((r, g, b))


def read_summary_colors(summary_csv: Path) -> dict[str, str]:
    """Map phylum name to computed ratio color."""
    colors: dict[str, str] = {}

    with summary_csv.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        required = {"phylum", "co2_fraction", "c1_fraction"}

        if not reader.fieldnames:
            raise ValueError("Summary CSV is missing a header row.")

        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise ValueError(f"Summary CSV is missing required columns: {', '.join(missing)}")

        for row in reader:
            phylum = row["phylum"].strip()
            if not phylum:
                continue

            co2 = float(row["co2_fraction"])
            c1 = float(row["c1_fraction"])
            total = co2 + c1
            if total <= 0:
                purple_weight = 0.5
            else:
                purple_weight = c1 / total

            colors[phylum] = interpolate_color(purple_weight)

    return colors


def parse_strip_template(strip_template: Path) -> tuple[list[str], list[tuple[str, str, str]]]:
    """Return header lines and DATA rows as (node_id, color, label)."""
    lines = strip_template.read_text(encoding="utf-8").splitlines()

    try:
        data_idx = lines.index("DATA")
    except ValueError as exc:
        raise ValueError("Template strip file is missing a DATA section.") from exc

    header = lines[: data_idx + 1]
    rows: list[tuple[str, str, str]] = []

    for line in lines[data_idx + 1 :]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        node_id, color, label = parts[0], parts[1], parts[2]
        rows.append((node_id, color, label))

    return header, rows


def build_outputs(
    summary_colors: dict[str, str],
    template_rows: list[tuple[str, str, str]],
) -> tuple[list[str], list[str]]:
    """Build DATA lines for strip and range-tree color files."""
    strip_data: list[str] = []
    tree_data: list[str] = []

    for node_id, _old_color, label in template_rows:
        if not label.startswith("p__"):
            continue

        phylum = label[3:]
        color = summary_colors.get(phylum)
        if color is None:
            continue

        strip_data.append(f"{node_id}\t{color}\t{phylum}")
        tree_data.append(f"{node_id}\trange\t{color}\t{phylum}")

    return strip_data, tree_data


def build_collapsed_tree_outputs(summary_colors: dict[str, str]) -> tuple[list[str], list[str]]:
    """Build strip/range DATA lines for a phylum-collapsed tree where leaf ids are phylum names."""
    strip_data: list[str] = []
    tree_data: list[str] = []

    for phylum in sorted(summary_colors):
        color = summary_colors[phylum]
        strip_data.append(f"{phylum}\t{color}\t{phylum}")
        tree_data.append(f"{phylum}\trange\t{color}\t{phylum}")

    return strip_data, tree_data


def write_strip_file(header: list[str], data_lines: list[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    content = "\n".join(header + data_lines) + "\n"
    output_path.write_text(content, encoding="utf-8")


def write_tree_file(data_lines: list[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    header = ["TREE_COLORS", "SEPARATOR TAB", "DATA"]
    content = "\n".join(header + data_lines) + "\n"
    output_path.write_text(content, encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate iTOL strip/tree color files for only the phyla present in "
            "a summary CSV, using a green(CO2)-purple(C1) ratio color blend."
        )
    )
    parser.add_argument(
        "summary_csv",
        nargs="?",
        type=Path,
        default=DEFAULT_SUMMARY,
        help=f"Summary CSV path (default: {DEFAULT_SUMMARY})",
    )
    parser.add_argument(
        "--template-strip",
        type=Path,
        default=DEFAULT_TEMPLATE_STRIP,
        help=f"Template strip file path (default: {DEFAULT_TEMPLATE_STRIP})",
    )
    parser.add_argument(
        "--output-strip",
        type=Path,
        default=DEFAULT_OUTPUT_STRIP,
        help=f"Output strip file path (default: {DEFAULT_OUTPUT_STRIP})",
    )
    parser.add_argument(
        "--output-tree",
        type=Path,
        default=DEFAULT_OUTPUT_TREE,
        help=f"Output tree colors file path (default: {DEFAULT_OUTPUT_TREE})",
    )
    parser.add_argument(
        "--output-strip-collapsed",
        type=Path,
        default=DEFAULT_OUTPUT_STRIP_COLLAPSED,
        help=(
            "Output strip file for phylum-collapsed tree (default: "
            f"{DEFAULT_OUTPUT_STRIP_COLLAPSED})"
        ),
    )
    parser.add_argument(
        "--output-tree-collapsed",
        type=Path,
        default=DEFAULT_OUTPUT_TREE_COLLAPSED,
        help=(
            "Output tree colors file for phylum-collapsed tree (default: "
            f"{DEFAULT_OUTPUT_TREE_COLLAPSED})"
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    summary_colors = read_summary_colors(args.summary_csv)
    header, template_rows = parse_strip_template(args.template_strip)
    strip_data, tree_data = build_outputs(summary_colors, template_rows)
    strip_collapsed_data, tree_collapsed_data = build_collapsed_tree_outputs(summary_colors)

    if not strip_data:
        raise ValueError("No matching phyla found between summary and template strip file.")

    write_strip_file(header, strip_data, args.output_strip)
    write_tree_file(tree_data, args.output_tree)
    write_strip_file(header, strip_collapsed_data, args.output_strip_collapsed)
    write_tree_file(tree_collapsed_data, args.output_tree_collapsed)

    print(f"Summary phyla: {len(summary_colors)}")
    print(f"Matched phyla in template: {len(strip_data)}")
    print(f"Wrote strip colors to {args.output_strip}")
    print(f"Wrote tree colors to {args.output_tree}")
    print(f"Wrote collapsed-tree strip colors to {args.output_strip_collapsed}")
    print(f"Wrote collapsed-tree tree colors to {args.output_tree_collapsed}")


if __name__ == "__main__":
    main()
