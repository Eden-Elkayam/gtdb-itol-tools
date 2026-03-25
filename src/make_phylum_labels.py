#!/usr/bin/env python3
"""Generate a phylum-only iTOL LABELS dataset from a mixed-rank labels file."""

from __future__ import annotations

import argparse
from pathlib import Path


DEFAULT_INPUT = Path("output/itol_labels.txt")
DEFAULT_OUTPUT = Path("phylum_labels_only.txt")


def generate_phylum_labels_only(input_path: Path, output_path: Path) -> tuple[int, int]:
    """Return (total_data_rows_seen, phylum_rows_written)."""
    text = input_path.read_text(encoding="utf-8")
    lines = text.splitlines()

    if len(lines) < 3 or lines[0].strip() != "LABELS":
        raise ValueError("Input file does not look like an iTOL LABELS dataset.")

    data_started = False
    total_rows = 0
    written_rows = 0
    out_lines = ["LABELS", "SEPARATOR TAB", "DATA"]

    for line in lines:
        if not data_started:
            if line.strip() == "DATA":
                data_started = True
            continue

        if not line.strip():
            continue

        parts = line.split("\t")
        if len(parts) < 2:
            continue

        total_rows += 1
        node_id = parts[0].strip()
        label = parts[1].strip()

        if not label.startswith("p__"):
            continue

        phylum_only = label.split(";", 1)[0]
        phylum_name = phylum_only[3:] if phylum_only.startswith("p__") else phylum_only
        out_lines.append(f"{node_id}\t{phylum_name}")
        written_rows += 1

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    return total_rows, written_rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a phylum-only LABELS file from output/itol_labels.txt."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input LABELS file (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output phylum-only LABELS file (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    total_rows, written_rows = generate_phylum_labels_only(args.input, args.output)
    print(f"Read {total_rows} label rows from {args.input}")
    print(f"Wrote {written_rows} phylum-only rows to {args.output}")


if __name__ == "__main__":
    main()
