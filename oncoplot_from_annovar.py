"""Oncoprint-style mutation matrix from ANNOVAR ``exonic_variant_function`` output.

Reads one ANNOVAR-style ``*.exonic_variant_function`` annotation file per
sample, restricts to a configurable set of protein-changing mutation
categories, selects the most frequently mutated genes, and renders a
gene x sample oncoprint where each cell is coloured proportionally by the
mix of mutation categories observed for that gene in that sample. Marginal
bar charts show the total (category-weighted) mutation burden per sample
(top) and per gene (left).

Input format
------------
A directory of ANNOVAR-style ``exonic_variant_function`` files (12
whitespace/tab-separated columns: record index, mutation type, gene
annotation, chromosome, start, end, ref, alt, zygosity, and two further
genotype-related columns -- the standard shape produced by ANNOVAR's
``annotate_variation.pl --geneanno`` combined with per-sample genotype
columns). One file per sample; the sample ID is taken from the file name
with its extension removed unless a sample sheet is supplied with
``--sample-sheet`` (CSV with columns: sample_id, file_name).

Usage
-----
    python oncoplot_from_annovar.py \\
        --input-dir data/annotated_variants \\
        --output-dir results/oncoplot \\
        --top-n-genes 20
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, Patch
import numpy as np
import pandas as pd

ANNOVAR_COLUMNS = [
    "record_id", "mutation_type", "gene_annotation", "chromosome",
    "start", "end", "ref", "alt", "zygosity", "filter_flag",
    "value_1", "value_2",
]

DEFAULT_MUTATION_TYPES = [
    "frameshift insertion",
    "frameshift deletion",
    "nonsynonymous SNV",
    "stopgain",
    "stoploss",
]

# Palette covering the common ANNOVAR exonic categories; only the
# categories selected via --mutation-types are actually used.
MUTATION_TYPE_COLORS = {
    "nonsynonymous SNV": "#89CFF0",
    "synonymous SNV": "#CFA4F7",
    "stopgain": "#FF7A85",
    "stoploss": "#FFD54F",
    "frameshift insertion": "#8DE4A0",
    "frameshift deletion": "#FFB37E",
    "nonframeshift deletion": "#FFF07C",
    "nonframeshift insertion": "#C3A68F",
}
DEFAULT_CELL_COLOR = "#CCCCCC"


def load_annotation_file(path: Path, sample_id: str, exclude_gene_prefixes: list[str]) -> pd.DataFrame:
    """Parse one ANNOVAR-style exonic_variant_function file into a tidy frame.

    One output row per (variant, gene) pair -- variants annotated against
    multiple overlapping genes are exploded into multiple rows.
    """
    data = pd.read_csv(
        path, sep="\t", header=None, na_values=".", dtype=str,
        engine="python", on_bad_lines="skip",
    )

    n_expected = len(ANNOVAR_COLUMNS)
    if data.shape[1] < n_expected:
        for i in range(data.shape[1], n_expected):
            data[i] = np.nan
    elif data.shape[1] > n_expected:
        data = data.iloc[:, :n_expected]
    data.columns = ANNOVAR_COLUMNS

    data["sample_id"] = sample_id

    # A cell can list multiple gene annotations separated by commas.
    data = data.assign(
        gene_annotation=data["gene_annotation"].str.split(",")
    ).explode("gene_annotation")

    # Each annotation looks like "GENE:NM_000000:exon4:c.123A>G:p.T41A" or
    # "GENE(...)"; the gene symbol is whatever precedes the first ":", "(",
    # or ",".
    first_token = data["gene_annotation"].str.split(":", n=1).str[0]
    data["gene"] = first_token.str.split(r"[(),]", n=1).str[0].str.strip()

    if exclude_gene_prefixes:
        pattern = "^(?:" + "|".join(re.escape(p) for p in exclude_gene_prefixes) + ")"
        data = data[~data["gene"].str.match(pattern, case=False, na=False)]

    return data[["sample_id", "gene", "mutation_type"]]


def load_all_samples(
    input_dir: Path,
    glob_pattern: str,
    sample_sheet: pd.DataFrame | None,
    exclude_gene_prefixes: list[str],
) -> pd.DataFrame:
    if sample_sheet is not None:
        file_specs = list(zip(sample_sheet["sample_id"], sample_sheet["file_name"]))
    else:
        file_specs = [(p.stem, p.name) for p in sorted(input_dir.glob(glob_pattern))]

    if not file_specs:
        raise FileNotFoundError(f"No annotation files found in {input_dir} matching {glob_pattern!r}")

    frames = [
        load_annotation_file(input_dir / file_name, sample_id, exclude_gene_prefixes)
        for sample_id, file_name in file_specs
    ]
    return pd.concat(frames, ignore_index=True), [s for s, _ in file_specs]


def select_top_genes(
    annotations: pd.DataFrame, mutation_types: list[str], top_n: int
) -> tuple[pd.DataFrame, list[str]]:
    filtered = annotations[annotations["mutation_type"].isin(mutation_types)]
    gene_totals = (
        filtered.groupby("gene").size().reset_index(name="total")
        .sort_values("total", ascending=False)
    )
    top_genes = gene_totals.head(top_n)["gene"].tolist()
    return filtered[filtered["gene"].isin(top_genes)], top_genes


def build_count_matrix(
    filtered: pd.DataFrame, genes: list[str], samples: list[str]
) -> tuple[np.ndarray, dict]:
    """Return a (gene x sample) total-count matrix and a nested count dict.

    The nested dict has the shape ``counts[gene][sample][mutation_type] = n``
    and is what the proportional oncoprint cells are drawn from.
    """
    per_cell = (
        filtered.groupby(["gene", "sample_id", "mutation_type"])
        .size().reset_index(name="count")
    )

    counts: dict = {gene: {} for gene in genes}
    for row in per_cell.itertuples(index=False):
        counts[row.gene].setdefault(row.sample_id, {})[row.mutation_type] = row.count

    count_matrix = np.zeros((len(genes), len(samples)))
    for i, gene in enumerate(genes):
        for j, sample in enumerate(samples):
            count_matrix[i, j] = sum(counts.get(gene, {}).get(sample, {}).values())

    return count_matrix, counts


def plot_oncoprint(
    genes: list[str],
    samples: list[str],
    counts: dict,
    mutation_types: list[str],
    color_map: dict,
    output_path: Path,
) -> None:
    """Draw the oncoprint: proportional stacked cells plus marginal bar charts."""
    count_matrix = np.zeros((len(genes), len(samples)))
    type_counts_by_sample = {m: np.zeros(len(samples)) for m in mutation_types}
    type_counts_by_gene = {m: np.zeros(len(genes)) for m in mutation_types}

    for i, gene in enumerate(genes):
        for j, sample in enumerate(samples):
            cell = counts.get(gene, {}).get(sample, {})
            cell_total = sum(cell.values())
            count_matrix[i, j] = cell_total
            if cell_total == 0:
                continue
            for mtype, value in cell.items():
                fraction = value / cell_total
                if mtype in type_counts_by_sample:
                    type_counts_by_sample[mtype][j] += fraction
                    type_counts_by_gene[mtype][i] += fraction

    presence = (count_matrix > 0).astype(int)
    gene_counts = presence.sum(axis=1)
    sample_counts = presence.sum(axis=0)

    fig = plt.figure(figsize=(max(4.2, 0.4 * len(samples)), max(5.5, 0.3 * len(genes))))
    gs = gridspec.GridSpec(
        2, 2, width_ratios=[1.2, 8], height_ratios=[1.0, 8], wspace=0.05, hspace=0.035
    )
    ax_top = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax_main = fig.add_subplot(gs[1, 1])

    # Main grid: proportionally stacked cell fills, then a uniform border on top.
    for i, gene in enumerate(genes):
        for j, sample in enumerate(samples):
            cell = counts.get(gene, {}).get(sample, {})
            cell_total = sum(cell.values())
            if cell_total > 0:
                y = i
                for mtype in mutation_types:
                    value = cell.get(mtype, 0)
                    if value > 0:
                        height = value / cell_total
                        ax_main.add_patch(Rectangle(
                            (j, y), 1, height,
                            facecolor=color_map.get(mtype, DEFAULT_CELL_COLOR),
                            edgecolor="none",
                        ))
                        y += height
            ax_main.add_patch(Rectangle((j, i), 1, 1, facecolor="none", edgecolor="black", linewidth=0.4))

    ax_main.set_xlim(0, len(samples))
    ax_main.set_ylim(0, len(genes))
    ax_main.invert_yaxis()
    ax_main.set_xticks(np.arange(len(samples)) + 0.5)
    ax_main.set_xticklabels(samples, rotation=90, ha="center", fontsize=11)
    ax_main.set_yticks([])
    for i, gene in enumerate(genes):
        ax_main.text(len(samples) + 0.2, i + 0.5, gene, va="center", fontsize=10)
    for spine in ax_main.spines.values():
        spine.set_visible(False)

    # Top margin: stacked bar of mutation-category burden per sample.
    bottom = np.zeros(len(samples))
    for mtype in mutation_types:
        ax_top.bar(
            np.arange(len(samples)) + 0.5, type_counts_by_sample[mtype],
            bottom=bottom, color=color_map.get(mtype, DEFAULT_CELL_COLOR),
        )
        bottom += type_counts_by_sample[mtype]
    ax_top.set_xlim(0, len(samples))
    ax_top.set_xticks([])
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)

    # Left margin: stacked horizontal bar of mutation-category burden per gene.
    y_pos = np.arange(len(genes)) + 0.5
    left = np.zeros(len(genes))
    for mtype in mutation_types:
        ax_left.barh(
            y_pos, type_counts_by_gene[mtype], left=left,
            color=color_map.get(mtype, DEFAULT_CELL_COLOR),
        )
        left += type_counts_by_gene[mtype]
    ax_left.set_ylim(0, len(genes))
    ax_left.invert_yaxis()
    max_gene_count = max(gene_counts.max(), 1)
    ax_left.set_xlim(0, max_gene_count)
    ax_left.invert_xaxis()
    ax_left.set_yticks([])
    ax_left.spines["left"].set_visible(False)
    ax_left.spines["top"].set_visible(False)

    legend_elements = [
        Patch(facecolor=color_map.get(m, DEFAULT_CELL_COLOR), edgecolor="black", label=m)
        for m in mutation_types
    ]
    ax_main.legend(handles=legend_elements, bbox_to_anchor=(1.4, 0.5), loc="center left", frameon=False)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--glob-pattern", default="*.exonic_variant_function")
    parser.add_argument(
        "--sample-sheet", type=Path, default=None,
        help="Optional CSV with columns: sample_id, file_name (overrides --glob-pattern discovery)",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/oncoplot"))
    parser.add_argument("--top-n-genes", type=int, default=20)
    parser.add_argument(
        "--mutation-types", nargs="+", default=DEFAULT_MUTATION_TYPES,
        help="ANNOVAR mutation categories to include (default: %(default)s)",
    )
    parser.add_argument(
        "--exclude-gene-prefixes", nargs="*", default=["LOC"],
        help="Gene symbol prefixes to drop, e.g. placeholder/uncharacterized loci (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sample_sheet = pd.read_csv(args.sample_sheet) if args.sample_sheet else None
    annotations, samples = load_all_samples(
        args.input_dir, args.glob_pattern, sample_sheet, args.exclude_gene_prefixes
    )

    filtered, top_genes = select_top_genes(annotations, args.mutation_types, args.top_n_genes)
    count_matrix, counts = build_count_matrix(filtered, top_genes, samples)

    matrix_df = pd.DataFrame(count_matrix, index=top_genes, columns=samples)
    matrix_df.to_csv(args.output_dir / "gene_by_sample_mutation_counts.csv")

    long_form = (
        filtered.groupby(["gene", "sample_id", "mutation_type"])
        .size().reset_index(name="count")
    )
    long_form.to_csv(args.output_dir / "gene_sample_mutation_type_counts.csv", index=False)

    plot_oncoprint(
        top_genes, samples, counts, args.mutation_types, MUTATION_TYPE_COLORS,
        args.output_dir / "oncoprint.png",
    )

    print(f"Wrote oncoprint results for {len(top_genes)} genes x {len(samples)} samples to: {args.output_dir}")


if __name__ == "__main__":
    main()
