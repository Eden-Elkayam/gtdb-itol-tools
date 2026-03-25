#!/usr/bin/env python3
"""Summarize filtered phylum rows into CO2/C1 totals and fractions."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


DEFAULT_OUTPUT = Path("output/phylum_co2_c1_summary.csv")
CO2_COLUMNS = ["CBB", "rTCA", "WL", "3HP", "3HP/4HB", "DC/4HB"]
C1_COLUMNS = ["Serine Cycle", "Reductive Glycine Pathway", "RuMP"]


def _to_float(value: str, column: str, phylum: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(
            f"Non-numeric value '{value}' in column '{column}' for phylum '{phylum}'."
        ) from exc


def _to_int(value: str, column: str, phylum: str) -> int:
    try:
        return int(float(value))
    except ValueError as exc:
        raise ValueError(
            f"Non-numeric value '{value}' in column '{column}' for phylum '{phylum}'."
        ) from exc


def summarize(input_csv: Path, output_csv: Path) -> int:
    required = {"phylum", "n_genomes", *CO2_COLUMNS, *C1_COLUMNS}
    aggregates: dict[str, dict[str, float | int | None]] = defaultdict(
        lambda: {"n_genomes": None, "co2_sum": 0.0, "c1_sum": 0.0}
    )
    warnings: list[str] = []

    with input_csv.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        if not reader.fieldnames:
            raise ValueError("Input CSV is missing a header row.")

        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise ValueError(f"Input CSV is missing required columns: {', '.join(missing)}")

        for row in reader:
            phylum = row["phylum"].strip()
            if not phylum:
                continue

            n_genomes = _to_int(row["n_genomes"], "n_genomes", phylum)
            co2_sum = sum(_to_float(row[col], col, phylum) for col in CO2_COLUMNS)
            c1_sum = sum(_to_float(row[col], col, phylum) for col in C1_COLUMNS)

            existing_n = aggregates[phylum]["n_genomes"]
            if existing_n is None:
                aggregates[phylum]["n_genomes"] = n_genomes
            elif int(existing_n) != n_genomes:
                warnings.append(
                    f"Warning: n_genomes inconsistency for phylum '{phylum}' "
                    f"({int(existing_n)} vs {n_genomes}). Keeping first value {int(existing_n)}."
                )

            aggregates[phylum]["co2_sum"] += co2_sum
            aggregates[phylum]["c1_sum"] += c1_sum

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", newline="", encoding="utf-8") as outfile:
        fieldnames = [
            "phylum",
            "n_genomes",
            "co2_sum",
            "c1_sum",
            "co2_fraction",
            "c1_fraction",
            "overall_fraction",
        ]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        count = 0
        for phylum in sorted(aggregates):
            n_genomes_raw = aggregates[phylum]["n_genomes"]
            n_genomes = float(n_genomes_raw if n_genomes_raw is not None else 0)
            co2_sum = float(aggregates[phylum]["co2_sum"])
            c1_sum = float(aggregates[phylum]["c1_sum"])

            if n_genomes <= 0:
                co2_fraction = 0.0
                c1_fraction = 0.0
                overall_fraction = 0.0
            else:
                co2_fraction = co2_sum / n_genomes
                c1_fraction = c1_sum / n_genomes
                overall_fraction = (co2_sum + c1_sum) / n_genomes

            writer.writerow(
                {
                    "phylum": phylum,
                    "n_genomes": int(n_genomes),
                    "co2_sum": format(co2_sum, ".10g"),
                    "c1_sum": format(c1_sum, ".10g"),
                    "co2_fraction": format(co2_fraction, ".10g"),
                    "c1_fraction": format(c1_fraction, ".10g"),
                    "overall_fraction": format(overall_fraction, ".10g"),
                }
            )
            count += 1

    if warnings:
        print(f"Warnings ({len(warnings)}):")
        for warning in warnings[:20]:
            print(f"- {warning}")
        if len(warnings) > 20:
            print(f"- ... and {len(warnings) - 20} more")

    return count


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create per-phylum summary with n_genomes, co2_sum, c1_sum, and "
            "fractions out of n_genomes from a filtered CSV file."
        )
    )
    parser.add_argument("input_csv", type=Path, help="Path to filtered input CSV.")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output CSV path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    n_phyla = summarize(args.input_csv, args.output)
    print(f"Wrote {n_phyla} phylum summaries to {args.output}")


if __name__ == "__main__":
    main()
