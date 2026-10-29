"""Shared enrichment terms between two functional enrichment result sets.

Functional enrichment tools (GO, KEGG, Reactome, ...) are typically run
separately per sample group and per ontology/category, producing one
result table per (category, group) pair. This script finds the terms
that came up as enriched in *both* groups for each category, and writes
the matching rows from both sides to a single CSV per category.

Input format
------------
For each category name in ``--categories`` (default: BP, CC, KEGG, MF, R
-- GO Biological Process / Cellular Component / Molecular Function, KEGG
pathways, and Reactome pathways), the input directory must contain two
CSV files:

    {category}_{group_a_suffix}.csv
    {category}_{group_b_suffix}.csv

each with at least a ``Term`` column. Output: ``{category}_shared_terms.csv``
with all rows (from both files) whose term is present in both.

Usage
-----
    python functional_enrichment_overlap.py \\
        --input-dir data/enrichment \\
        --output-dir results/enrichment_overlap
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

DEFAULT_CATEGORIES = ["BP", "CC", "KEGG", "MF", "R"]


def shared_terms_for_category(
    input_dir: Path, category: str, group_a_suffix: str, group_b_suffix: str
) -> tuple[pd.DataFrame, int] | None:
    a_path = input_dir / f"{category}_{group_a_suffix}.csv"
    b_path = input_dir / f"{category}_{group_b_suffix}.csv"
    if not a_path.exists() or not b_path.exists():
        return None

    a_df = pd.read_csv(a_path)
    b_df = pd.read_csv(b_path)
    a_df["_term_clean"] = a_df["Term"].astype(str).str.strip()
    b_df["_term_clean"] = b_df["Term"].astype(str).str.strip()

    shared = set(a_df["_term_clean"]) & set(b_df["_term_clean"])

    a_shared = a_df[a_df["_term_clean"].isin(shared)].copy()
    b_shared = b_df[b_df["_term_clean"].isin(shared)].copy()
    a_shared.insert(0, "Source", f"{category}_{group_a_suffix}")
    b_shared.insert(0, "Source", f"{category}_{group_b_suffix}")

    result = pd.concat([a_shared, b_shared], ignore_index=True).drop(columns=["_term_clean"])
    return result, len(shared)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/enrichment_overlap"))
    parser.add_argument("--categories", nargs="+", default=DEFAULT_CATEGORIES)
    parser.add_argument("--group-a-suffix", default="A")
    parser.add_argument("--group-b-suffix", default="B")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for category in args.categories:
        outcome = shared_terms_for_category(args.input_dir, category, args.group_a_suffix, args.group_b_suffix)
        if outcome is None:
            print(f"{category}: input files not found, skipped")
            continue
        result, n_shared = outcome
        output_path = args.output_dir / f"{category}_shared_terms.csv"
        result.to_csv(output_path, index=False)
        print(f"{category}: {n_shared} shared terms saved to {output_path}")


if __name__ == "__main__":
    main()
