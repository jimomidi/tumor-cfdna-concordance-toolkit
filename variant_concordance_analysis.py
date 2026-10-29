"""Paired-sample variant concordance from VCFs (e.g. tumor tissue vs cfDNA).

For each paired sample (e.g. a tumor biopsy and a matched cell-free DNA
draw from the same case), this script compares the two variant call sets
and reports:

  * the number of variants shared, and unique to each side,
  * concordance, expressed either as a fraction of the reference
    sample's variants (``--metric overlap``, the default) or as the
    Jaccard index of the two variant sets (``--metric jaccard``),
  * three summary plots: concordance per pair, a stacked bar of
    shared / reference-only / other-only variant counts, and a clean
    dot plot of concordance per pair.

A variant is defined as the (chromosome, position, reference allele,
alternate allele) tuple, so this is a comparison of the raw variant
calls rather than of any downstream annotation.

Input format
------------
A pairing sheet (CSV) with one row per case and the columns:

    pair_id, reference_vcf, other_vcf

``reference_vcf`` and ``other_vcf`` may be absolute paths or paths
relative to the pairing sheet's location. Which physical sample type
(tumor, plasma, ...) plays "reference" vs "other" is up to the sample
sheet -- the script itself is agnostic to sample type labels.

Usage
-----
    python variant_concordance_analysis.py \\
        --pairing-sheet pairs.csv \\
        --output-dir results/concordance \\
        --metric overlap \\
        --normalize-chromosomes
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from vcf_utils import read_variants


def compare_pair(
    reference_vcf: Path, other_vcf: Path, normalize_chrom: bool
) -> dict:
    """Compute shared/unique variant counts and concordance for one pair."""
    reference_variants = set(read_variants(reference_vcf, normalize_chrom))
    other_variants = set(read_variants(other_vcf, normalize_chrom))

    shared = reference_variants & other_variants
    reference_only = reference_variants - other_variants
    other_only = other_variants - reference_variants
    union = reference_variants | other_variants

    overlap_pct = (len(shared) / len(reference_variants) * 100) if reference_variants else 0.0
    jaccard_pct = (len(shared) / len(union) * 100) if union else 0.0

    return {
        "Shared": len(shared),
        "Reference_only": len(reference_only),
        "Other_only": len(other_only),
        "Union": len(union),
        "Concordance_overlap_pct": overlap_pct,
        "Concordance_jaccard_pct": jaccard_pct,
    }


def build_concordance_table(
    pairing_sheet: pd.DataFrame, base_dir: Path, normalize_chrom: bool
) -> pd.DataFrame:
    rows = []
    for _, row in pairing_sheet.iterrows():
        stats = compare_pair(
            base_dir / row["reference_vcf"],
            base_dir / row["other_vcf"],
            normalize_chrom=normalize_chrom,
        )
        rows.append({"Pair": row["pair_id"], **stats})
    return pd.DataFrame(rows)


def plot_concordance_bar(table: pd.DataFrame, metric_col: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(table["Pair"], table[metric_col])
    ax.set_title("Variant Concordance per Pair")
    ax.set_ylabel("Concordance (%)")
    ax.set_xlabel("Pair")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_variant_distribution(table: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(table["Pair"], table["Shared"], label="Shared")
    ax.bar(table["Pair"], table["Reference_only"], bottom=table["Shared"], label="Reference only")
    ax.bar(
        table["Pair"],
        table["Other_only"],
        bottom=table["Shared"] + table["Reference_only"],
        label="Other only",
    )
    ax.set_title("Variant Distribution per Pair")
    ax.set_ylabel("Number of Variants")
    ax.set_xlabel("Pair")
    ax.legend()
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_concordance_dotplot(table: pd.DataFrame, metric_col: str, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.5, 3))
    ax.scatter(table["Pair"], table[metric_col], s=180, edgecolor="black", zorder=3)

    y_max = table[metric_col].max()
    offset = y_max * 0.07 if y_max else 1.0
    for x, y in zip(table["Pair"], table[metric_col]):
        ax.text(x, y + offset, f"{y:.1f}", ha="center", va="bottom", fontsize=10)

    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    ax.set_ylim(0, y_max * 1.25 if y_max else 1.0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("Variant Concordance per Pair")
    ax.set_ylabel("Concordance (%)")
    ax.set_xlabel("Pair")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pairing-sheet", type=Path, required=True,
        help="CSV with columns: pair_id, reference_vcf, other_vcf",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/concordance"),
        help="Directory for output tables and figures (default: %(default)s)",
    )
    parser.add_argument(
        "--metric", choices=["overlap", "jaccard"], default="overlap",
        help=(
            "Concordance metric to plot: 'overlap' = shared / reference "
            "variants; 'jaccard' = shared / union. Both are always written "
            "to the output table. (default: %(default)s)"
        ),
    )
    parser.add_argument(
        "--normalize-chromosomes", action="store_true",
        help="Strip 'chr' prefixes and unify 'MT'/'M' before comparing variants.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pairing_sheet = pd.read_csv(args.pairing_sheet)
    base_dir = args.pairing_sheet.parent

    table = build_concordance_table(pairing_sheet, base_dir, args.normalize_chromosomes)
    table.to_csv(args.output_dir / "variant_concordance.csv", index=False)

    metric_col = f"Concordance_{args.metric}_pct"
    plot_concordance_bar(table, metric_col, args.output_dir / "concordance_bar.png")
    plot_variant_distribution(table, args.output_dir / "variant_distribution.png")
    plot_concordance_dotplot(table, metric_col, args.output_dir / "concordance_dotplot.png")

    print(f"Wrote variant concordance results to: {args.output_dir}")


if __name__ == "__main__":
    main()
