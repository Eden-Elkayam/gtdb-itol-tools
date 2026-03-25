#!/usr/bin/env python3
"""Run the donor-filtered GTDB/iTOL pipeline in one command."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path


# ========================= USER INPUT (EDIT THIS BLOCK) =========================
# Run `python3 run_donor_threshold_pipeline.py` after setting these values.
# You can still override any of them via CLI flags if you prefer.
USER_INPUT = {
    "donor_csv": "c_by_e_donor_phylum.csv",
    "tree": "output/bac120.tree_stripped",
    "taxonomy": "bac120_taxonomy.tsv",
    "template_strip": "output/itol_dataset_strip_phylum.txt",
    # Set to "all" for all donors, or a list like ["formate oxidation", "H2 oxidation"].
    "donors": "all",
    # Legacy single-donor key (still supported if "donors" is omitted).
    "donor": "formate oxidation",
    "overall_fraction_threshold": 0.2,
    "output_dir": "pipeline_0.3_threshold",
}
# ================================================================================


PATHWAY_COLUMNS = [
    "CBB",
    "rTCA",
    "WL",
    "3HP",
    "3HP/4HB",
    "DC/4HB",
    "Serine Cycle",
    "Reductive Glycine Pathway",
    "RuMP",
]


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-zA-Z0-9]+", "_", value.strip()).strip("_").lower()
    return slug or "value"


def threshold_tag(value: float) -> str:
    text = f"{value:g}"
    return text.replace("-", "neg").replace(".", "p")


def filter_and_aggregate_by_phylum(
    input_csv: Path,
    output_csv: Path,
    donors: set[str] | None,
    warnings: list[str],
) -> tuple[int, int]:
    """Filter by donor and aggregate to one row per phylum.

    Pathway columns are summed across selected donor rows.
    n_genomes is kept as-is per phylum (must be consistent across rows).
    Returns (kept_rows, unique_phyla).
    """
    required_columns = {"phylum", "domain", "n_genomes", "electron_donor", *PATHWAY_COLUMNS}

    aggregates: dict[str, dict[str, object]] = {}

    with input_csv.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        if not reader.fieldnames:
            raise ValueError("Input CSV is missing a header row.")

        missing = sorted(required_columns.difference(reader.fieldnames))
        if missing:
            raise ValueError(f"Input CSV is missing required columns: {', '.join(missing)}")

        output_csv.parent.mkdir(parents=True, exist_ok=True)

        kept_rows = 0
        for row in reader:
            donor_value = row.get("electron_donor")
            if donors is not None and donor_value not in donors:
                continue

            phylum = row["phylum"].strip()
            if not phylum:
                continue

            domain = row["domain"].strip()

            try:
                n_genomes = int(float(row["n_genomes"]))
            except ValueError as exc:
                raise ValueError(
                    f"Non-numeric n_genomes '{row['n_genomes']}' for phylum '{phylum}'."
                ) from exc

            if phylum not in aggregates:
                aggregates[phylum] = {
                    "phylum": phylum,
                    "domain": domain,
                    "n_genomes": n_genomes,
                    "pathways": {col: 0.0 for col in PATHWAY_COLUMNS},
                }
            else:
                existing_n = int(aggregates[phylum]["n_genomes"])
                if n_genomes != existing_n:
                    warnings.append(
                        "Warning: n_genomes inconsistency for phylum "
                        f"'{phylum}' ({existing_n} vs {n_genomes}). Keeping first value {existing_n}."
                    )

                existing_domain = str(aggregates[phylum]["domain"])
                if domain != existing_domain:
                    warnings.append(
                        "Warning: domain inconsistency for phylum "
                        f"'{phylum}' ({existing_domain} vs {domain}). Keeping first value {existing_domain}."
                    )

            pathway_acc = aggregates[phylum]["pathways"]
            assert isinstance(pathway_acc, dict)
            for col in PATHWAY_COLUMNS:
                raw = row[col]
                try:
                    pathway_acc[col] += float(raw)
                except ValueError as exc:
                    raise ValueError(
                        f"Non-numeric value '{raw}' in pathway column '{col}' for phylum '{phylum}'."
                    ) from exc

            kept_rows += 1

    output_columns = ["phylum", "domain", "n_genomes", *PATHWAY_COLUMNS]
    with output_csv.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=output_columns)
        writer.writeheader()

        for phylum in sorted(aggregates):
            row = aggregates[phylum]
            pathway_acc = row["pathways"]
            assert isinstance(pathway_acc, dict)

            out = {
                "phylum": str(row["phylum"]),
                "domain": str(row["domain"]),
                "n_genomes": int(row["n_genomes"]),
            }
            for col in PATHWAY_COLUMNS:
                out[col] = format(float(pathway_acc[col]), ".10g")
            writer.writerow(out)

    return kept_rows, len(aggregates)


def filter_summary_by_overall_fraction(
    summary_csv: Path,
    output_csv: Path,
    overall_fraction_threshold: float,
) -> int:
    with summary_csv.open("r", newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile)
        if not reader.fieldnames:
            raise ValueError("Summary CSV is missing a header row.")
        if "overall_fraction" not in reader.fieldnames:
            raise ValueError("Summary CSV is missing required column: overall_fraction")
        if "phylum" not in reader.fieldnames:
            raise ValueError("Summary CSV is missing required column: phylum")

        output_csv.parent.mkdir(parents=True, exist_ok=True)
        with output_csv.open("w", newline="", encoding="utf-8") as outfile:
            writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
            writer.writeheader()

            kept = 0
            for row in reader:
                raw_fraction = row.get("overall_fraction", "")
                try:
                    overall_fraction = float(raw_fraction)
                except ValueError as exc:
                    phylum = row.get("phylum", "")
                    raise ValueError(
                        f"Non-numeric overall_fraction '{raw_fraction}' for phylum '{phylum}'."
                    ) from exc

                if overall_fraction <= overall_fraction_threshold:
                    continue

                writer.writerow(row)
                kept += 1

    return kept


def run_cmd(args: list[str]) -> None:
    print("Running:", " ".join(args))
    subprocess.run(args, check=True)


def donors_tag(donors: set[str] | None) -> str:
    if donors is None:
        return "all_donors"
    joined = "_or_".join(sorted(donors))
    return slugify(joined)


def normalize_donors_input(raw: object) -> set[str] | None:
    """Convert donor input into a donor set, or None to represent all donors."""
    if raw is None:
        return None

    if isinstance(raw, str):
        value = raw.strip()
        if not value:
            return set()
        if value.lower() == "all":
            return None
        return {value}

    if isinstance(raw, list):
        donors = {str(item).strip() for item in raw if str(item).strip()}
        if len(donors) == 1 and "all" in {d.lower() for d in donors}:
            return None
        return donors

    raise ValueError("Donor input must be a string, list, or None.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter by donor(s), aggregate by phylum, summarize fractions, "
            "and generate iTOL color files in one output folder."
        )
    )
    parser.add_argument("--donor-csv", type=Path, required=False, help="Input donor CSV.")
    parser.add_argument("--tree", type=Path, required=False, help="Input tree Newick file.")
    parser.add_argument("--taxonomy", type=Path, required=False, help="Input taxonomy TSV file.")
    parser.add_argument(
        "--template-strip",
        type=Path,
        required=False,
        help="Input iTOL strip template file.",
    )
    parser.add_argument("--donor", required=False, help="Single target electron donor string.")
    parser.add_argument(
        "--donors",
        nargs="+",
        required=False,
        help="One or more donor names. Example: --donors 'formate oxidation' 'H2 oxidation'",
    )
    parser.add_argument(
        "--all-donors",
        action="store_true",
        help="Use all donor categories present in the input CSV.",
    )
    parser.add_argument(
        "--overall-fraction-threshold",
        type=float,
        required=False,
        help="Keep phyla only when overall_fraction is strictly greater than this threshold.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=False,
        help="Output folder path (created if it does not exist).",
    )
    return parser.parse_args()


def resolve_params(args: argparse.Namespace) -> argparse.Namespace:
    """Use CLI values when provided; otherwise fallback to USER_INPUT block."""
    script_root = Path(__file__).resolve().parent

    donor_csv = args.donor_csv or (script_root / USER_INPUT["donor_csv"])
    tree = args.tree or (script_root / USER_INPUT["tree"])
    taxonomy = args.taxonomy or (script_root / USER_INPUT["taxonomy"])
    template_strip = args.template_strip or (script_root / USER_INPUT["template_strip"])

    if args.all_donors:
        donors = None
    elif args.donors:
        donors = normalize_donors_input(args.donors)
    elif args.donor:
        donors = normalize_donors_input(args.donor)
    elif "donors" in USER_INPUT:
        donors = normalize_donors_input(USER_INPUT["donors"])
    else:
        donors = normalize_donors_input(USER_INPUT.get("donor"))

    overall_fraction_threshold = (
        args.overall_fraction_threshold
        if args.overall_fraction_threshold is not None
        else float(USER_INPUT["overall_fraction_threshold"])
    )
    output_dir = args.output_dir or (script_root / USER_INPUT["output_dir"])

    resolved = argparse.Namespace(
        donor_csv=Path(donor_csv),
        tree=Path(tree),
        taxonomy=Path(taxonomy),
        template_strip=Path(template_strip),
        donors=donors,
        overall_fraction_threshold=float(overall_fraction_threshold),
        output_dir=Path(output_dir),
    )

    missing: list[str] = []
    if not resolved.donor_csv.exists():
        missing.append(f"--donor-csv / USER_INPUT['donor_csv'] -> {resolved.donor_csv}")
    if not resolved.tree.exists():
        missing.append(f"--tree / USER_INPUT['tree'] -> {resolved.tree}")
    if not resolved.taxonomy.exists():
        missing.append(f"--taxonomy / USER_INPUT['taxonomy'] -> {resolved.taxonomy}")
    if not resolved.template_strip.exists():
        missing.append(
            f"--template-strip / USER_INPUT['template_strip'] -> {resolved.template_strip}"
        )
    if resolved.donors is not None and not resolved.donors:
        missing.append("--donor/--donors or USER_INPUT['donors'] must include at least one donor")

    if missing:
        raise ValueError(
            "Missing or invalid required inputs:\n- " + "\n- ".join(missing)
        )

    return resolved


def main() -> None:
    args = resolve_params(parse_args())

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    donor_slug = donors_tag(args.donors)
    overall_frac_tag = threshold_tag(args.overall_fraction_threshold)

    filtered_csv = output_dir / f"thresholded_{donor_slug}.csv"
    summary_csv = output_dir / f"summary_{donor_slug}.csv"
    summary_filtered_csv = (
        output_dir / f"summary_{donor_slug}_overall_fraction_gt_{overall_frac_tag}.csv"
    )
    pruned_tree = (
        output_dir / f"bac120_tree_pruned_{donor_slug}_overall_fraction_gt_{overall_frac_tag}.nwk"
    )
    strip_file = (
        output_dir
        / f"itol_dataset_strip_phylum_co2_c1_ratio_{donor_slug}_overall_fraction_gt_{overall_frac_tag}.txt"
    )
    tree_colors_file = (
        output_dir
        / f"itol_tree_colours_co2_c1_ratio_{donor_slug}_overall_fraction_gt_{overall_frac_tag}.txt"
    )
    strip_file_collapsed = (
        output_dir
        / f"itol_dataset_strip_phylum_co2_c1_ratio_{donor_slug}_overall_fraction_gt_{overall_frac_tag}_collapsed_tree.txt"
    )
    tree_colors_file_collapsed = (
        output_dir
        / f"itol_tree_colours_co2_c1_ratio_{donor_slug}_overall_fraction_gt_{overall_frac_tag}_collapsed_tree.txt"
    )

    warnings: list[str] = []
    kept_rows, unique_phyla = filter_and_aggregate_by_phylum(
        input_csv=args.donor_csv,
        output_csv=filtered_csv,
        donors=args.donors,
        warnings=warnings,
    )
    if kept_rows == 0:
        raise ValueError("No rows matched donor filtering. Nothing else to run.")

    root = Path(__file__).resolve().parent
    prune_script = root / "src/prune_tree_by_phyla.py"
    summarize_script = root / "src/summarize_phyla_fractions.py"
    colors_script = root / "src/generate_itol_phylum_ratio_colors.py"

    run_cmd([sys.executable, str(summarize_script), str(filtered_csv), "--output", str(summary_csv)])

    kept_phyla = filter_summary_by_overall_fraction(
        summary_csv=summary_csv,
        output_csv=summary_filtered_csv,
        overall_fraction_threshold=args.overall_fraction_threshold,
    )
    if kept_phyla == 0:
        raise ValueError(
            "No phyla passed overall_fraction threshold. "
            "Lower --overall-fraction-threshold or adjust upstream filters."
        )

    run_cmd(
        [
            sys.executable,
            str(prune_script),
            str(summary_filtered_csv),
            "--tree",
            str(args.tree),
            "--taxonomy",
            str(args.taxonomy),
            "--output",
            str(pruned_tree),
        ]
    )

    run_cmd(
        [
            sys.executable,
            str(colors_script),
            str(summary_filtered_csv),
            "--template-strip",
            str(args.template_strip),
            "--output-strip",
            str(strip_file),
            "--output-tree",
            str(tree_colors_file),
            "--output-strip-collapsed",
            str(strip_file_collapsed),
            "--output-tree-collapsed",
            str(tree_colors_file_collapsed),
        ]
    )

    print("Pipeline complete.")
    if args.donors is None:
        print("Donors used: all donors")
    else:
        print(f"Donors used ({len(args.donors)}): {', '.join(sorted(args.donors))}")
    print(f"Rows matched donor filter: {kept_rows}")
    print(f"Unique phyla after donor aggregation: {unique_phyla}")
    if warnings:
        print(f"Warnings ({len(warnings)}):")
        for warning in warnings[:20]:
            print(f"- {warning}")
        if len(warnings) > 20:
            print(f"- ... and {len(warnings) - 20} more")
    print(f"Overall-fraction threshold: > {args.overall_fraction_threshold}")
    print(f"Filtered CSV: {filtered_csv}")
    print(f"Summary CSV: {summary_csv}")
    print(f"Summary filtered by overall_fraction: {summary_filtered_csv}")
    print(f"Pruned tree: {pruned_tree}")
    print(f"iTOL strip colors: {strip_file}")
    print(f"iTOL tree colors: {tree_colors_file}")
    print(f"iTOL strip colors (phylum-collapsed tree): {strip_file_collapsed}")
    print(f"iTOL tree colors (phylum-collapsed tree): {tree_colors_file_collapsed}")


if __name__ == "__main__":
    main()
