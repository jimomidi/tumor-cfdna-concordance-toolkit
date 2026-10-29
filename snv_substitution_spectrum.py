"""Single-nucleotide substitution spectrum (SBS6-style) from paired VCFs.

Given a set of per-sample VCF files, each labelled with a sample group
(e.g. "tumor" vs "cfdna"), this script counts single-base substitutions
and collapses them onto the 6 canonical, strand-symmetric substitution
classes (C>A, C>G, C>T, T>A, T>C, T>G) -- the same classes used by the
COSMIC single-base-substitution (SBS) signatures. It reports:

  1. The overall substitution spectrum across all samples.
  2. The per-sample substitution counts (stacked contribution plot).
  3. The spectrum broken down by sample group, with group totals
     compared side by side.

Input format
------------
A sample sheet (CSV) with one row per VCF file and at least the columns:

    sample_id, group, vcf_path

``group`` is a free-text label (e.g. "tumor", "cfdna") used to split
samples into two comparison groups; ``vcf_path`` may be absolute or
relative to the sample sheet's location.

Usage
-----
    python snv_substitution_spectrum.py \\
        --sample-sheet samples.csv \\
        --output-dir results/snv_spectrum
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

from vcf_utils import SNV_CLASS_ORDER, classify_snv, read_snv_variants


def count_substitution_classes(vcf_path: Path) -> pd.Series:
    """Count the 6 canonical substitution classes observed in one VCF."""
    counts = pd.Series(0, index=SNV_CLASS_ORDER, dtype=int)
    for variant in read_snv_variants(vcf_path):
        snv_class = classify_snv(variant.ref, variant.alt)
        if snv_class is not None:
            counts[snv_class] += 1
    return counts


def build_spectrum_table(sample_sheet: pd.DataFrame, base_dir: Path) -> pd.DataFrame:
    """Return a (substitution class x sample) count matrix.

    Column order follows the sample sheet; row order follows
    ``SNV_CLASS_ORDER``.
    """
    columns = {}
    for _, row in sample_sheet.iterrows():
        vcf_path = base_dir / row["vcf_path"]
        columns[row["sample_id"]] = count_substitution_classes(vcf_path)
    return pd.DataFrame(columns).reindex(SNV_CLASS_ORDER)


def plot_overall_spectrum(spectrum: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    spectrum.sum(axis=1).plot(kind="bar", ax=ax)
    ax.set_title("Substitution Spectrum (SBS6)")
    ax.set_ylabel("Mutation Counts")
    ax.set_xlabel("Substitution Class")
    ax.set_xticklabels(spectrum.index, rotation=0)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_per_sample_contribution(spectrum: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    spectrum.T.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title("Substitution Class Contribution per Sample")
    ax.set_ylabel("Counts")
    ax.set_xlabel("Samples")
    ax.legend(title="Substitution Class", bbox_to_anchor=(1.05, 1))
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_group_comparison(
    spectrum: pd.DataFrame, sample_sheet: pd.DataFrame, output_path: Path
) -> None:
    """Side-by-side bar chart comparing the spectrum summed within each group."""
    groups = sample_sheet.set_index("sample_id")["group"]
    group_totals = spectrum.T.groupby(groups).sum().T

    x = range(len(spectrum.index))
    width = 0.8 / max(len(group_totals.columns), 1)

    fig, ax = plt.subplots(figsize=(8, 5))
    for i, group_name in enumerate(group_totals.columns):
        offsets = [xi + i * width for xi in x]
        ax.bar(offsets, group_totals[group_name], width=width, label=group_name)

    ax.set_xticks([xi + width * (len(group_totals.columns) - 1) / 2 for xi in x])
    ax.set_xticklabels(spectrum.index)
    ax.set_title("Substitution Spectrum by Sample Group")
    ax.set_ylabel("Mutation Counts")
    ax.set_xlabel("Substitution Class")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample-sheet", type=Path, required=True,
        help="CSV with columns: sample_id, group, vcf_path",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("results/snv_spectrum"),
        help="Directory for output tables and figures (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sample_sheet = pd.read_csv(args.sample_sheet)
    base_dir = args.sample_sheet.parent

    spectrum = build_spectrum_table(sample_sheet, base_dir)
    spectrum.to_csv(args.output_dir / "substitution_spectrum_counts.csv")

    plot_overall_spectrum(spectrum, args.output_dir / "overall_spectrum.png")
    plot_per_sample_contribution(spectrum, args.output_dir / "per_sample_contribution.png")
    if "group" in sample_sheet.columns and sample_sheet["group"].nunique() > 1:
        plot_group_comparison(spectrum, sample_sheet, args.output_dir / "group_comparison.png")

    print(f"Wrote substitution spectrum results to: {args.output_dir}")


if __name__ == "__main__":
    main()
