"""Genome-wide copy-number variation (CNV) profiling across multiple samples.

Consumes per-sample copy-number "ratio" tracks and discrete CNV segment
calls (e.g. Control-FREEC-style ``*_ratio.txt`` and ``*_CNVs`` output),
lays every sample out along a shared, concatenated genome-position axis,
and produces:

  * a genome-wide log2-ratio scatter plot per sample,
  * a GISTIC-style gain/loss segment plot per sample,
  * a combined multi-sample overview figure (one profile track per sample
    plus a stacked "CNV landscape" row),
  * recurrent CNV regions, both as exact coordinate matches and as
    reciprocal-overlap-merged clusters across samples,
  * summary statistics (event counts, lengths, per-sample burden, top
    recurrent regions),
  * a multi-sheet Excel workbook collecting all of the above tables.

Input format
------------
A sample sheet (CSV) with columns:

    sample_id, group, sample_type, ratio_file, cnv_file

``ratio_file`` is a tab-separated file with at least the columns
``Chromosome, Start, Ratio, MedianRatio, CopyNumber`` (as produced by
Control-FREEC); ``cnv_file`` is a headerless, tab-separated file with
columns ``Chromosome, Start, End, CopyNumber, CNV_Type`` (``CNV_Type`` is
expected to be "gain" or "loss"). Paths may be absolute or relative to the
sample sheet's location.

Chromosome names are used as-is after stripping a leading "chr"/"Chr". If
your reference uses accession-style contig names (e.g. RefSeq
``NC_XXXXXX.X``) instead of chromosome numbers, supply
``--chrom-map`` -- a two-column TSV of ``accession<TAB>chromosome_label`` --
to translate them before plotting.

Usage
-----
    python cnv_profiling_pipeline.py \\
        --sample-sheet samples.csv \\
        --output-dir results/cnv \\
        --chrom-order chrom_order.txt
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GAIN_COLOR = "#E91E63"
LOSS_COLOR = "#1E88E5"


# --------------------------------------------------------------------------
# Loading and chromosome-axis setup
# --------------------------------------------------------------------------

def _natural_sort_key(label: str):
    """Sort '1', '2', ..., '10', 'X', 'Y', 'M' in the conventional order."""
    return (0, int(label)) if label.isdigit() else (1, label)


def load_chrom_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}
    mapping = {}
    with open(path) as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                mapping[parts[0]] = parts[1]
    return mapping


def normalize_chromosome(value: str, chrom_map: dict[str, str]) -> str:
    value = str(value).replace("chr", "").replace("Chr", "").strip()
    return chrom_map.get(value, value)


def read_ratio_file(path: Path, sample_id: str, group: str, sample_type: str, chrom_map: dict) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.replace("chrChromosome", "Chromosome") for c in df.columns]

    df["Chromosome"] = df["Chromosome"].apply(lambda x: normalize_chromosome(x, chrom_map))
    for col in ["Start", "Ratio", "MedianRatio", "CopyNumber"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df[df["Ratio"] > 0].copy()
    df["Log2Ratio"] = np.log2(df["Ratio"])
    df["SampleID"], df["Group"], df["SampleType"] = sample_id, group, sample_type
    df["SourceFile"] = Path(path).name
    return df


def read_cnv_file(path: Path, sample_id: str, group: str, sample_type: str, chrom_map: dict) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None)
    df.columns = ["Chromosome", "Start", "End", "CopyNumber", "CNV_Type"]

    df["Chromosome"] = df["Chromosome"].apply(lambda x: normalize_chromosome(x, chrom_map))
    for col in ["Start", "End", "CopyNumber"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["CNV_Type"] = df["CNV_Type"].astype(str).str.lower()

    df["SampleID"], df["Group"], df["SampleType"] = sample_id, group, sample_type
    df["SourceFile"] = Path(path).name
    return df


def load_samples(sample_sheet: pd.DataFrame, base_dir: Path, chrom_map: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    ratio_frames, cnv_frames = [], []
    for row in sample_sheet.itertuples(index=False):
        ratio_frames.append(
            read_ratio_file(base_dir / row.ratio_file, row.sample_id, row.group, row.sample_type, chrom_map)
        )
        cnv_frames.append(
            read_cnv_file(base_dir / row.cnv_file, row.sample_id, row.group, row.sample_type, chrom_map)
        )
    return pd.concat(ratio_frames, ignore_index=True), pd.concat(cnv_frames, ignore_index=True)


def compute_genome_axis(
    ratio_all: pd.DataFrame, cnv_all: pd.DataFrame, chrom_order: list[str] | None
) -> tuple[pd.DataFrame, pd.DataFrame, list[float], list[str]]:
    """Concatenate per-chromosome coordinates into one genome-wide axis.

    Returns the ratio/CNV tables with added genome-position columns, plus
    tick positions and labels for plotting.
    """
    if chrom_order is None:
        chrom_order = sorted(ratio_all["Chromosome"].unique(), key=_natural_sort_key)

    chrom_lengths = (
        ratio_all.groupby("Chromosome")["Start"].max().reindex(chrom_order).dropna().to_dict()
    )

    offsets, tick_pos, tick_lab = {}, [], []
    current = 0
    for chrom in chrom_order:
        if chrom in chrom_lengths:
            offsets[chrom] = current
            tick_pos.append(current + chrom_lengths[chrom] / 2)
            tick_lab.append(chrom)
            current += chrom_lengths[chrom]

    ratio_all = ratio_all.copy()
    cnv_all = cnv_all.copy()
    ratio_all["GenomePosition"] = ratio_all["Start"] + ratio_all["Chromosome"].map(offsets)
    cnv_all["GenomeStart"] = cnv_all["Start"] + cnv_all["Chromosome"].map(offsets)
    cnv_all["GenomeEnd"] = cnv_all["End"] + cnv_all["Chromosome"].map(offsets)

    return ratio_all, cnv_all, tick_pos, tick_lab


# --------------------------------------------------------------------------
# Per-sample plots
# --------------------------------------------------------------------------

def plot_genome_wide_profile(
    sample_df: pd.DataFrame, title: str, tick_pos, tick_lab, output_path: Path,
    gain_threshold: float, loss_threshold: float,
) -> None:
    sdf = sample_df.dropna(subset=["GenomePosition", "Log2Ratio"])
    fig, ax = plt.subplots(figsize=(16, 4.5))
    ax.scatter(sdf["GenomePosition"], sdf["Log2Ratio"], s=4, alpha=0.7)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.axhline(gain_threshold, color="red", linestyle="--", linewidth=0.8)
    ax.axhline(loss_threshold, color="blue", linestyle="--", linewidth=0.8)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lab, rotation=90)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Log2 ratio")
    ax.set_title(f"Genome-wide CNV profile: {title}")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_gistic_style(
    sample_df: pd.DataFrame, title: str, tick_pos, tick_lab, output_path: Path
) -> None:
    sdf = sample_df.dropna(subset=["GenomeStart", "GenomeEnd"])
    fig, ax = plt.subplots(figsize=(16, 3.8))
    ax.axhline(0, color="black", linewidth=0.8)
    for row in sdf.itertuples(index=False):
        y = 1 if row.CNV_Type == "gain" else -1
        color = "red" if row.CNV_Type == "gain" else "blue"
        ax.hlines(y, row.GenomeStart, row.GenomeEnd, colors=color, linewidth=4)
    ax.set_yticks([-1, 0, 1])
    ax.set_yticklabels(["Loss", "Neutral", "Gain"])
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lab, rotation=90)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("CNV event")
    ax.set_title(f"GISTIC-style CNV segment plot: {title}")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_combined_overview(
    ratio_all: pd.DataFrame, cnv_all: pd.DataFrame, sample_order: list[str],
    tick_pos, tick_lab, output_path: Path,
    guide_threshold: float = 0.25,
) -> None:
    """One profile row per sample plus a shared CNV-landscape summary row."""
    fig, axes = plt.subplots(
        nrows=len(sample_order) + 1, ncols=1, figsize=(18, 1.8 * len(sample_order) + 3),
        sharex=True,
        gridspec_kw={"height_ratios": [1] * len(sample_order) + [2.4]},
    )

    for ax, sample_id in zip(axes[:-1], sample_order):
        sdf = ratio_all[ratio_all["SampleID"] == sample_id].dropna(subset=["GenomePosition", "Log2Ratio"])
        sdf = sdf.sort_values("GenomePosition")
        ax.plot(sdf["GenomePosition"], sdf["Log2Ratio"], linewidth=1.0, alpha=0.9, color="#9400D4")
        ax.axhline(0, color="black", linewidth=0.6)
        ax.axhline(guide_threshold, color=GAIN_COLOR, linewidth=1.2)
        ax.axhline(-guide_threshold, color=LOSS_COLOR, linewidth=1.2)
        ax.set_ylim(-3.5, 3.0)
        ax.set_ylabel(sample_id, fontsize=14, fontweight="bold", labelpad=14)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    landscape_ax = axes[-1]
    y_positions = {sample_id: i for i, sample_id in enumerate(sample_order)}
    for sample_id in sample_order:
        sdf = cnv_all[cnv_all["SampleID"] == sample_id].dropna(subset=["GenomeStart", "GenomeEnd"])
        y = y_positions[sample_id]
        for row in sdf.itertuples(index=False):
            color = GAIN_COLOR if row.CNV_Type == "gain" else (LOSS_COLOR if row.CNV_Type == "loss" else "gray")
            landscape_ax.hlines(y=y, xmin=row.GenomeStart, xmax=row.GenomeEnd, colors=color, linewidth=6, alpha=0.9)

    for x in tick_pos:
        landscape_ax.axvline(x, color="lightgray", linewidth=0.4, alpha=0.5)
    landscape_ax.set_yticks(list(y_positions.values()))
    landscape_ax.set_yticklabels(sample_order, fontsize=12)
    landscape_ax.set_ylabel("CNV\nlandscape", fontsize=14, fontweight="bold")
    landscape_ax.set_ylim(-1, len(sample_order))
    landscape_ax.invert_yaxis()
    landscape_ax.spines["top"].set_visible(False)
    landscape_ax.spines["right"].set_visible(False)

    axes[-1].set_xticks(tick_pos)
    axes[-1].set_xticklabels(tick_lab, rotation=90, fontsize=11)
    axes[-1].set_xlabel("Chromosome", fontsize=14, fontweight="bold")

    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
# Recurrent regions and summary statistics
# --------------------------------------------------------------------------

def recurrent_regions_exact(cnv_all: pd.DataFrame) -> pd.DataFrame:
    """Regions recurring at exactly the same coordinates across samples."""
    n_samples = cnv_all["SampleID"].nunique()
    table = (
        cnv_all.groupby(["Chromosome", "Start", "End", "CNV_Type"])
        .agg(
            SampleCount=("SampleID", "nunique"),
            Samples=("SampleID", lambda x: ", ".join(sorted(set(x)))),
            Groups=("Group", lambda x: ", ".join(sorted(set(x)))),
        )
        .reset_index()
    )
    table["Frequency_All_Samples"] = table["SampleCount"] / n_samples
    return table


def recurrent_regions_by_overlap(cnv_all: pd.DataFrame, overlap_threshold: float = 0.50) -> pd.DataFrame:
    """Cluster CNV events by reciprocal overlap, independent of exact breakpoints.

    Two events are merged into the same region if the length of their
    overlap is at least ``overlap_threshold`` of *both* event lengths
    (reciprocal overlap), following the convention commonly used to call
    "the same" CNV across imprecisely-breakpointed samples.
    """
    n_samples = cnv_all["SampleID"].nunique()
    events = cnv_all.dropna(subset=["Chromosome", "Start", "End", "CNV_Type", "SampleID"]).copy()
    events["Start"] = events["Start"].astype(int)
    events["End"] = events["End"].astype(int)
    events["CNV_Type"] = events["CNV_Type"].astype(str).str.lower()

    merged_regions = []
    for (chrom, cnv_type), group in events.groupby(["Chromosome", "CNV_Type"]):
        group = group.sort_values(["Start", "End"]).reset_index(drop=True)
        used = set()

        for i, anchor in group.iterrows():
            if i in used:
                continue
            cluster_indices = [i]
            used.add(i)

            for j, candidate in group.iterrows():
                if j in used:
                    continue
                overlap_len = max(0, min(anchor["End"], candidate["End"]) - max(anchor["Start"], candidate["Start"]))
                len_anchor = anchor["End"] - anchor["Start"]
                len_candidate = candidate["End"] - candidate["Start"]
                if len_anchor <= 0 or len_candidate <= 0:
                    continue
                reciprocal_overlap = min(overlap_len / len_anchor, overlap_len / len_candidate)
                if reciprocal_overlap >= overlap_threshold:
                    cluster_indices.append(j)
                    used.add(j)

            cluster = group.loc[cluster_indices]
            merged_regions.append({
                "Chromosome": chrom,
                "Start": int(cluster["Start"].min()),
                "End": int(cluster["End"].max()),
                "CNV_Type": cnv_type,
                "SampleCount": cluster["SampleID"].nunique(),
                "Samples": ", ".join(sorted(cluster["SampleID"].dropna().unique())),
                "Groups": ", ".join(sorted(cluster["Group"].dropna().unique())),
                "SampleTypes": ", ".join(sorted(cluster["SampleType"].dropna().unique())),
                "Frequency_All_Samples": cluster["SampleID"].nunique() / n_samples,
            })

    return pd.DataFrame(merged_regions)


def format_recurrent_region_table(overlap_regions: pd.DataFrame) -> pd.DataFrame:
    table = overlap_regions.copy()
    table["Frequency"] = (table["Frequency_All_Samples"] * 100).round(1).astype(str) + "%"
    table["Type"] = table["CNV_Type"].replace({"gain": "Gain", "loss": "Loss"})
    table["Region"] = (
        "chr" + table["Chromosome"].astype(str)
        + ":" + table["Start"].astype(int).astype(str)
        + "-" + table["End"].astype(int).astype(str)
    )
    table["Annotated_Genes"] = ""  # placeholder for a downstream gene-annotation step
    columns = ["Region", "Type", "Frequency", "SampleCount", "Samples", "Groups", "SampleTypes", "Annotated_Genes"]
    return table[columns].sort_values(["SampleCount", "Region"], ascending=[False, True])


def summarize_cnv_burden(cnv_all: pd.DataFrame) -> dict[str, pd.DataFrame]:
    stats = cnv_all.dropna(subset=["SampleID", "Group", "SampleType", "Chromosome", "Start", "End", "CNV_Type"]).copy()
    stats["Start"] = stats["Start"].astype(int)
    stats["End"] = stats["End"].astype(int)
    stats["Length_Mb"] = (stats["End"] - stats["Start"]) / 1e6
    stats["CNV_Type"] = stats["CNV_Type"].astype(str).str.lower()

    overall = stats.groupby("CNV_Type").agg(
        EventCount=("CNV_Type", "size"),
        UniqueSamples=("SampleID", "nunique"),
        TotalLength_Mb=("Length_Mb", "sum"),
        MedianLength_Mb=("Length_Mb", "median"),
    ).reset_index()

    per_sample = (
        stats.groupby(["SampleID", "Group", "SampleType", "CNV_Type"]).size()
        .reset_index(name="EventCount")
        .pivot_table(index=["SampleID", "Group", "SampleType"], columns="CNV_Type", values="EventCount", fill_value=0)
        .reset_index()
    )
    per_sample.columns.name = None
    for col in ["gain", "loss"]:
        if col not in per_sample.columns:
            per_sample[col] = 0
    length_per_sample = stats.groupby(["SampleID", "Group", "SampleType"]).agg(
        Total_CNV_length_Mb=("Length_Mb", "sum")
    ).reset_index()
    per_sample = per_sample.merge(length_per_sample, on=["SampleID", "Group", "SampleType"], how="left")
    per_sample["Total_CNV_events"] = per_sample["gain"] + per_sample["loss"]

    by_sample_type = stats.groupby(["SampleType", "CNV_Type"]).agg(
        EventCount=("CNV_Type", "size"),
        TotalLength_Mb=("Length_Mb", "sum"),
        MedianLength_Mb=("Length_Mb", "median"),
    ).reset_index()

    top_chromosomes = (
        stats.groupby(["Chromosome", "CNV_Type"]).agg(
            EventCount=("CNV_Type", "size"),
            SampleCount=("SampleID", "nunique"),
            TotalLength_Mb=("Length_Mb", "sum"),
        ).reset_index()
        .sort_values(["SampleCount", "EventCount"], ascending=[False, False]).head(10)
    )

    return {
        "overall_by_cnv_type": overall,
        "per_sample_burden": per_sample,
        "by_sample_type": by_sample_type,
        "top_chromosomes": top_chromosomes,
    }


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample-sheet", type=Path, required=True,
        help="CSV with columns: sample_id, group, sample_type, ratio_file, cnv_file",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/cnv"))
    parser.add_argument(
        "--chrom-map", type=Path, default=None,
        help="Optional TSV of accession<TAB>chromosome_label, for accession-style contig names",
    )
    parser.add_argument(
        "--chrom-order", type=Path, default=None,
        help="Optional text file, one chromosome label per line, giving the plotting order "
             "(default: natural-sorted order of the chromosome labels found in the data)",
    )
    parser.add_argument("--overlap-threshold", type=float, default=0.50,
                         help="Reciprocal overlap fraction for merging recurrent CNV regions (default: %(default)s)")
    parser.add_argument("--min-recurrent-sample-count", type=int, default=3)
    parser.add_argument("--top-n-recurrent-regions", type=int, default=20)
    parser.add_argument("--gain-threshold", type=float, default=np.log2(1.25),
                         help="Log2-ratio guide line for gains in per-sample profile plots (default: log2(1.25))")
    parser.add_argument("--loss-threshold", type=float, default=np.log2(0.75),
                         help="Log2-ratio guide line for losses in per-sample profile plots (default: log2(0.75))")
    parser.add_argument("--overview-guide-threshold", type=float, default=0.25,
                         help="Flat log2-ratio guide line used in the combined overview figure (default: %(default)s)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    profile_dir = args.output_dir / "genome_wide_profiles"
    gistic_dir = args.output_dir / "gistic_style_plots"
    profile_dir.mkdir(exist_ok=True)
    gistic_dir.mkdir(exist_ok=True)

    sample_sheet = pd.read_csv(args.sample_sheet)
    base_dir = args.sample_sheet.parent
    chrom_map = load_chrom_map(args.chrom_map)

    ratio_all, cnv_all = load_samples(sample_sheet, base_dir, chrom_map)

    chrom_order = None
    if args.chrom_order:
        chrom_order = [line.strip() for line in open(args.chrom_order) if line.strip()]
    ratio_all, cnv_all, tick_pos, tick_lab = compute_genome_axis(ratio_all, cnv_all, chrom_order)

    sample_order = list(sample_sheet["sample_id"])

    for sample_id in sample_order:
        info = sample_sheet.set_index("sample_id").loc[sample_id]
        title = f"{sample_id} | {info['group']} | {info['sample_type']}"
        plot_genome_wide_profile(
            ratio_all[ratio_all["SampleID"] == sample_id], title, tick_pos, tick_lab,
            profile_dir / f"{sample_id}_genome_wide_profile.png",
            args.gain_threshold, args.loss_threshold,
        )
        plot_gistic_style(
            cnv_all[cnv_all["SampleID"] == sample_id], title, tick_pos, tick_lab,
            gistic_dir / f"{sample_id}_gistic_style.png",
        )

    plot_combined_overview(
        ratio_all, cnv_all, sample_order, tick_pos, tick_lab,
        args.output_dir / "combined_overview.png", args.overview_guide_threshold,
    )

    exact_regions = recurrent_regions_exact(cnv_all)
    overlap_regions = recurrent_regions_by_overlap(cnv_all, args.overlap_threshold)
    formatted_table = format_recurrent_region_table(overlap_regions)

    top_recurrent = exact_regions[
        exact_regions["SampleCount"] >= args.min_recurrent_sample_count
    ].sort_values(["SampleCount", "Chromosome", "Start"], ascending=[False, True, True]).head(args.top_n_recurrent_regions)

    summary = summarize_cnv_burden(cnv_all)

    excel_path = args.output_dir / "cnv_analysis_summary.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        sample_sheet.to_excel(writer, sheet_name="sample_annotation", index=False)
        ratio_all.to_excel(writer, sheet_name="all_ratio_data", index=False)
        cnv_all.to_excel(writer, sheet_name="all_cnv_events", index=False)
        exact_regions.to_excel(writer, sheet_name="recurrent_regions_exact", index=False)
        overlap_regions.to_excel(writer, sheet_name="recurrent_regions_overlap", index=False)
        formatted_table.to_excel(writer, sheet_name="recurrent_regions_formatted", index=False)
        top_recurrent.to_excel(writer, sheet_name="top_recurrent_regions", index=False)
        for name, table in summary.items():
            table.to_excel(writer, sheet_name=name[:31], index=False)

    print(f"Wrote CNV profiling results to: {args.output_dir}")
    print(f"  Samples: {len(sample_order)}")
    print(f"  Total CNV events: {len(cnv_all)}")
    print(f"  Excel summary: {excel_path}")


if __name__ == "__main__":
    main()
