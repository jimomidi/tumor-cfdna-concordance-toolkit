"""Cross-sample gene-list overlap analysis for two sample groups.

Starting from one gene-list file per sample, this script builds a
combined gene matrix and then characterizes how genes are shared within
and between two sample groups (e.g. "cfdna" vs "tumor"), producing:

  1. ``gene_matrix.csv`` -- one column per sample, the deduplicated gene
     list read from that sample's input file.
  2. ``gene_frequency.csv`` -- every gene observed, which samples it was
     observed in, and how many samples that is.
  3. ``group_overlap.csv`` -- genes observed in at least one sample of
     *both* groups.
  4. ``group_unique_genes.csv`` -- genes unique to each group (i.e. not
     observed in the other group), side by side.
  5. ``pairwise_multiway_overlap.csv`` -- for every group, every 2-way up
     to N-way combination of its samples, the genes shared by exactly
     that combination or more (i.e. the intersection of that combination).
  6. Printed (not saved) exclusive N-way reports: for every combination
     of samples within a group, how many genes are observed in exactly
     that combination and no others.

``--exclude-gene-prefixes`` is applied consistently to every output
above, including the pairwise/multiway and exclusive-N-way reports.

Input format
------------
A directory containing one CSV file per sample. Each file is expected to
have the gene symbols in its first column (extra columns, if any, are
ignored); the sample name is taken from the file name (without
extension). Sample names are also used to decide group membership via
``--group-a-prefix`` / ``--group-b-prefix`` (case-insensitive "starts
with" match).

Usage
-----
    python gene_set_overlap_analysis.py \\
        --input-dir data/gene_lists \\
        --output-dir results/gene_overlap \\
        --group-a-prefix cfdna --group-b-prefix tumor \\
        --exclude-gene-prefixes LOC
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import pandas as pd


def read_gene_list(path: Path) -> list[str]:
    """Read the deduplicated (case-insensitive), stripped gene list from a CSV."""
    first_column = pd.read_csv(path).iloc[:, 0]
    genes = first_column.dropna().astype(str).str.strip()
    genes = genes[genes != ""]

    seen, unique_genes = set(), []
    for gene in genes:
        key = gene.lower()
        if key not in seen:
            seen.add(key)
            unique_genes.append(gene)
    return unique_genes


def build_gene_matrix(input_dir: Path, glob_pattern: str = "*.csv") -> pd.DataFrame:
    files = sorted(input_dir.glob(glob_pattern))
    if not files:
        raise FileNotFoundError(f"No gene-list files found in {input_dir} matching {glob_pattern!r}")
    columns = {path.stem: pd.Series(read_gene_list(path)) for path in files}
    return pd.DataFrame(columns)


def apply_prefix_exclusion(genes: list[str], exclude_prefixes: list[str]) -> list[str]:
    if not exclude_prefixes:
        return genes
    prefixes = tuple(p.upper() for p in exclude_prefixes)
    return [g for g in genes if not g.strip().upper().startswith(prefixes)]


def gene_frequency_table(matrix: pd.DataFrame, exclude_prefixes: list[str]) -> pd.DataFrame:
    gene_to_samples: dict[str, dict] = {}
    for sample in matrix.columns:
        genes = matrix[sample].dropna().astype(str).str.strip()
        genes = genes[genes != ""]
        for gene in genes:
            key = gene.lower()
            entry = gene_to_samples.setdefault(key, {"Gene": gene, "Samples": []})
            if sample not in entry["Samples"]:
                entry["Samples"].append(sample)

    rows = []
    for entry in gene_to_samples.values():
        gene = entry["Gene"].strip()
        if apply_prefix_exclusion([gene], exclude_prefixes) == []:
            continue
        rows.append({"Gene": gene, "Samples": ", ".join(entry["Samples"]), "Frequency": len(entry["Samples"])})

    return pd.DataFrame(rows).sort_values("Frequency", ascending=False).reset_index(drop=True)


def _gene_set(matrix: pd.DataFrame, columns: list[str], exclude_prefixes: list[str]) -> dict[str, str]:
    """Map lower-cased gene -> original-cased gene, pooled across ``columns``."""
    gene_dict: dict[str, str] = {}
    for col in columns:
        genes = matrix[col].dropna().astype(str).str.strip()
        genes = apply_prefix_exclusion(list(genes[genes != ""]), exclude_prefixes)
        for gene in genes:
            gene_dict.setdefault(gene.lower(), gene)
    return gene_dict


def group_columns(matrix: pd.DataFrame, prefix: str) -> list[str]:
    return [c for c in matrix.columns if c.lower().startswith(prefix.lower())]


def group_vs_group_overlap(
    matrix: pd.DataFrame, group_a_cols: list[str], group_b_cols: list[str], exclude_prefixes: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    a_genes = _gene_set(matrix, group_a_cols, exclude_prefixes)
    b_genes = _gene_set(matrix, group_b_cols, exclude_prefixes)

    shared_keys = sorted(set(a_genes) & set(b_genes))
    overlap_df = pd.DataFrame({"Shared_Genes": [a_genes[k] for k in shared_keys]})

    a_key_set, b_key_set = set(a_genes), set(b_genes)
    a_unique = sorted(g for k, g in a_genes.items() if k not in b_key_set)
    b_unique = sorted(g for k, g in b_genes.items() if k not in a_key_set)
    max_len = max(len(a_unique), len(b_unique))
    a_list = a_unique + [None] * (max_len - len(a_unique))
    b_list = b_unique + [None] * (max_len - len(b_unique))
    unique_df = pd.DataFrame({"Group_A_Unique_Genes": a_list, "Group_B_Unique_Genes": b_list})

    return overlap_df, unique_df


def pairwise_multiway_overlap(
    matrix: pd.DataFrame, columns: list[str], group_name: str, exclude_prefixes: list[str]
) -> list[dict]:
    gene_sets = {col: _gene_set(matrix, [col], exclude_prefixes) for col in columns}
    rows = []
    for n in range(2, len(columns) + 1):
        for combo in itertools.combinations(columns, n):
            overlap_keys = set(gene_sets[combo[0]])
            for col in combo[1:]:
                overlap_keys &= set(gene_sets[col])
            overlap_genes = sorted(gene_sets[combo[0]][k] for k in overlap_keys)
            rows.append({
                "Group": group_name,
                "Overlap_Level": f"{n}-way",
                "Samples": ", ".join(combo),
                "Number_of_Genes": len(overlap_genes),
                "Genes": ", ".join(overlap_genes),
            })
    return rows


def exclusive_nway_report(matrix: pd.DataFrame, columns: list[str], group_name: str, exclude_prefixes: list[str]) -> str:
    """Genes observed in *exactly* each combination of samples, and no others."""
    gene_sets = {col: set(_gene_set(matrix, [col], exclude_prefixes)) for col in columns}
    all_genes = set().union(*gene_sets.values()) if gene_sets else set()
    presence = {gene: [c for c in columns if gene in gene_sets[c]] for gene in all_genes}

    lines = [f"\n{'=' * 80}", group_name, "=" * 80]
    for n in range(len(columns), 0, -1):
        lines.append(f"\nEXCLUSIVE {n}-WAY")
        lines.append("-" * 80)
        for combo in itertools.combinations(columns, n):
            combo_set = set(combo)
            count = sum(1 for present_in in presence.values() if set(present_in) == combo_set)
            lines.append(f"{' + '.join(combo)} : {count} genes")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/gene_overlap"))
    parser.add_argument("--group-a-prefix", default="cfdna")
    parser.add_argument("--group-b-prefix", default="tumor")
    parser.add_argument(
        "--exclude-gene-prefixes", nargs="*", default=[],
        help="Gene symbol prefixes to drop from all outputs, e.g. placeholder loci (default: none)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matrix = build_gene_matrix(args.input_dir)
    matrix.to_csv(args.output_dir / "gene_matrix.csv", index=False)

    frequency_table = gene_frequency_table(matrix, args.exclude_gene_prefixes)
    frequency_table.to_csv(args.output_dir / "gene_frequency.csv", index=False)

    group_a_cols = group_columns(matrix, args.group_a_prefix)
    group_b_cols = group_columns(matrix, args.group_b_prefix)
    if not group_a_cols or not group_b_cols:
        raise ValueError(
            f"Could not find sample columns for both groups "
            f"(group A prefix={args.group_a_prefix!r}: {group_a_cols}, "
            f"group B prefix={args.group_b_prefix!r}: {group_b_cols})"
        )

    overlap_df, unique_df = group_vs_group_overlap(
        matrix, group_a_cols, group_b_cols, args.exclude_gene_prefixes
    )
    overlap_df.to_csv(args.output_dir / "group_overlap.csv", index=False)
    unique_df.to_csv(args.output_dir / "group_unique_genes.csv", index=False)

    multiway_rows = pairwise_multiway_overlap(matrix, group_a_cols, "Group_A", args.exclude_gene_prefixes)
    multiway_rows += pairwise_multiway_overlap(matrix, group_b_cols, "Group_B", args.exclude_gene_prefixes)
    pd.DataFrame(multiway_rows).to_csv(args.output_dir / "pairwise_multiway_overlap.csv", index=False)

    report_path = args.output_dir / "exclusive_nway_report.txt"
    with open(report_path, "w") as handle:
        handle.write(exclusive_nway_report(matrix, group_a_cols, "Group_A", args.exclude_gene_prefixes))
        handle.write(exclusive_nway_report(matrix, group_b_cols, "Group_B", args.exclude_gene_prefixes))

    print(f"Wrote gene overlap results to: {args.output_dir}")


if __name__ == "__main__":
    main()
