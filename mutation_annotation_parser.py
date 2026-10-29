"""Parse free-text mutation annotations into a gene x category count matrix.

Some report/spreadsheet exports summarize mutations per gene as free text
in a single cell, e.g. ``"C>T(3), C>A"`` or an indel count folded into a
"Mutation Type" column such as ``"frameshift deletion(2)"``. This script
parses that text into structured (gene, mutation category, count) records,
aggregates them into a gene x category matrix, cross-checks the total
against a reported "Number of mutations" column (if present), and plots
mutation frequency by gene.

Input format
------------
A CSV with, at minimum, the columns:

    Gene name, Sample, Mutation Type, <mutation detail column>

``Gene name`` and ``Sample`` may be sparsely filled (blank cells are
forward-filled from the row above, as is typical when these columns come
from merged cells in a spreadsheet export). The mutation detail column
holds substitution notation like ``"C>T"`` or ``"C>T(3)"``; if no
substitution notation is found in a row, the row is counted as one
"Insertion/Deletion" event (with a count parsed from a parenthesized
number in ``Mutation Type`` if present). An optional ``Number of
mutations`` column, if present, is used only to sanity-check that parsed
counts sum to the reported total.

Usage
-----
    python mutation_annotation_parser.py \\
        --input-csv data/mutation_report.csv \\
        --output-dir results/mutation_annotation \\
        --mutation-detail-column "Unnamed: 4"
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

SNV_PATTERN = re.compile(r"([ACGT]\s*>\s*[ACGT])\s*(?:\(\s*(\d+)\s*\))?")
COUNT_IN_PARENS_PATTERN = re.compile(r"\(\s*(\d+)\s*\)")

PREFERRED_CATEGORY_ORDER = [
    "C>T", "C>A", "C>G", "T>C", "T>A", "T>G",
    "A>T", "A>C", "A>G", "G>T", "G>C", "G>A",
    "Insertion/Deletion",
]


def parse_mutation_records(
    df: pd.DataFrame, gene_col: str, sample_col: str, type_col: str, detail_col: str
) -> pd.DataFrame:
    df = df.copy()
    df[gene_col] = df[gene_col].ffill()
    df[sample_col] = df[sample_col].ffill()

    records = []
    for _, row in df.iterrows():
        gene, value = row[gene_col], row[detail_col]
        if pd.isna(gene) or pd.isna(value):
            continue
        value = str(value).strip()

        snv_matches = SNV_PATTERN.findall(value)
        if snv_matches:
            for mutation, count_str in snv_matches:
                records.append({
                    "Gene name": gene,
                    "Mutation category": mutation.replace(" ", ""),
                    "Count": int(count_str) if count_str else 1,
                })
        else:
            count_match = COUNT_IN_PARENS_PATTERN.search(str(row.get(type_col, "")))
            records.append({
                "Gene name": gene,
                "Mutation category": "Insertion/Deletion",
                "Count": int(count_match.group(1)) if count_match else 1,
            })

    return pd.DataFrame(records)


def build_pivot_table(records: pd.DataFrame, gene_order: list[str]) -> pd.DataFrame:
    summary = records.groupby(["Gene name", "Mutation category"], as_index=False)["Count"].sum()
    pivot = summary.pivot(index="Gene name", columns="Mutation category", values="Count").fillna(0)
    pivot = pivot.reindex(gene_order)

    existing = [c for c in PREFERRED_CATEGORY_ORDER if c in pivot.columns]
    other = [c for c in pivot.columns if c not in existing]
    return pivot[existing + other]


def plot_mutation_frequency(pivot: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(pivot.index)), 7))
    pivot.plot(kind="bar", ax=ax, width=0.85)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Mutation frequency")
    ax.set_title("Mutation Frequency by Gene and Category")
    ax.tick_params(axis="x", rotation=90)
    ax.legend(title="Mutation category", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/mutation_annotation"))
    parser.add_argument("--gene-column", default="Gene name")
    parser.add_argument("--sample-column", default="Sample")
    parser.add_argument("--mutation-type-column", default="Mutation Type")
    parser.add_argument(
        "--mutation-detail-column", default="Unnamed: 4",
        help=(
            "Column holding free-text substitution notation such as 'C>T(3)'. "
            "Spreadsheet-to-CSV exports often leave this column unnamed; if the "
            "given name is not found, the 5th column (index 4) is used instead. "
            "(default: %(default)r)"
        ),
    )
    parser.add_argument("--total-mutations-column", default="Number of mutations")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input_csv)

    detail_col = args.mutation_detail_column
    if detail_col not in df.columns:
        detail_col = df.columns[4]
        print(f"Column {args.mutation_detail_column!r} not found; using column index 4 ({detail_col!r}) instead.")

    records = parse_mutation_records(df, args.gene_column, args.sample_column, args.mutation_type_column, detail_col)
    records.to_csv(args.output_dir / "mutation_records.csv", index=False)

    # Forward-fill first, same as parse_mutation_records does internally: a
    # gene whose *first* row happens to be a blank merged-cell (rather than a
    # later row) would otherwise be silently dropped from gene_order -- and
    # therefore from the reindexed pivot table below -- despite its mutations
    # being present in records.
    gene_order = df[args.gene_column].ffill().dropna().drop_duplicates().tolist()
    pivot = build_pivot_table(records, gene_order)
    pivot.to_csv(args.output_dir / "mutation_frequency_by_gene.csv")

    if args.total_mutations_column in df.columns:
        expected_total = int(df[args.total_mutations_column].dropna().sum())
        observed_total = int(records["Count"].sum())
        status = "OK: totals match." if expected_total == observed_total else (
            f"WARNING: totals do not match (difference = {expected_total - observed_total})."
        )
        print(f"Expected total: {expected_total}, observed total: {observed_total}. {status}")

    plot_mutation_frequency(pivot, args.output_dir / "mutation_frequency_by_gene.png")
    print(f"Wrote mutation annotation results to: {args.output_dir}")


if __name__ == "__main__":
    main()
