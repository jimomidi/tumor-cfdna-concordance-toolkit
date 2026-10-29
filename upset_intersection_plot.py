"""Two-group UpSet-style intersection plot from a precomputed intersection table.

Draws a custom UpSet-style figure for two sample groups (e.g. "plasma"
samples vs "tumor" samples), each with any number of member samples:

  * a top bar chart with one pair of bars per intersection combination
    (group A size, group B size, side by side),
  * a bottom dot matrix showing, for each combination, which samples
    from each group participate in it (connected by a vertical line),
  * a left horizontal bar chart of each sample's total set size.

This script does not compute intersections itself -- it visualizes a
table of intersection sizes that you have already computed (e.g. from
variant sets, gene sets, or any other per-sample collection). This keeps
the plot reusable regardless of what kind of set is being compared.

Input format
------------
A CSV with one row per intersection combination and the columns:

    <group-a-size-column>, <group-a-prefix>_<sample>, ...
    <group-b-size-column>, <group-b-prefix>_<sample>, ...

The size columns hold the intersection's total size for that group (e.g.
number of shared variants); the per-sample columns hold 0/1 membership
flags indicating whether that sample participates in the combination.
Column names starting with ``--group-a-prefix`` / ``--group-b-prefix``
are treated as that group's per-sample membership columns.

Usage
-----
    python upset_intersection_plot.py \\
        --intersection-table intersections.csv \\
        --output results/upset_plot.png \\
        --group-a-prefix plasma --group-a-size-column plasma_size \\
        --group-b-prefix tumor --group-b-size-column tumor_size
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

GROUP_A_COLOR = "#D6A6FF"
GROUP_B_COLOR = "#90EE90"


def load_intersection_table(
    path: Path, group_a_size_col: str, group_b_size_col: str, group_a_prefix: str, group_b_prefix: str
) -> tuple[pd.DataFrame, list[str], list[str]]:
    df = pd.read_csv(path)
    group_a_samples = [c for c in df.columns if c.startswith(group_a_prefix) and c != group_a_size_col]
    group_b_samples = [c for c in df.columns if c.startswith(group_b_prefix) and c != group_b_size_col]

    df = df.copy()
    df["_sort_total"] = df[group_a_size_col] + df[group_b_size_col]
    df = df.sort_values(by=["_sort_total", group_a_size_col, group_b_size_col], ascending=False).reset_index(drop=True)

    return df, group_a_samples, group_b_samples


def plot_upset(
    df: pd.DataFrame,
    group_a_size_col: str, group_b_size_col: str,
    group_a_samples: list[str], group_b_samples: list[str],
    group_a_label: str, group_b_label: str,
    output_path: Path,
) -> None:
    all_samples = group_a_samples + group_b_samples
    n = len(df)
    x = list(range(n))
    bar_width = 0.46

    set_sizes = {s: df.loc[df[s] == 1, group_a_size_col].sum() for s in group_a_samples}
    set_sizes.update({s: df.loc[df[s] == 1, group_b_size_col].sum() for s in group_b_samples})

    fig = plt.figure(figsize=(max(10, 0.8 * n), 8))
    gs = GridSpec(2, 2, width_ratios=[1.2, 4], height_ratios=[3.6, 1.8], wspace=0.05, hspace=0.05, figure=fig)
    ax_set = fig.add_subplot(gs[1, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_bar)

    a_x = [i - bar_width / 2 for i in x]
    b_x = [i + bar_width / 2 for i in x]
    ax_bar.bar(a_x, df[group_a_size_col], width=bar_width, color=GROUP_A_COLOR, edgecolor="black", label=group_a_label)
    ax_bar.bar(b_x, df[group_b_size_col], width=bar_width, color=GROUP_B_COLOR, edgecolor="black", label=group_b_label)
    ax_bar.set_ylabel("Intersection size")
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_bar.legend(handles=[
        Patch(facecolor=GROUP_A_COLOR, edgecolor="black", label=group_a_label),
        Patch(facecolor=GROUP_B_COLOR, edgecolor="black", label=group_b_label),
    ], frameon=False, loc="upper right")

    y_map = {sample: i for i, sample in enumerate(all_samples)}
    for i in x:
        for sample in all_samples:
            ax_mat.scatter(i, y_map[sample], s=90, color="#d9d9d9", zorder=1)

    for group_samples, color in [(group_a_samples, GROUP_A_COLOR), (group_b_samples, GROUP_B_COLOR)]:
        for i, row in df.iterrows():
            active = [s for s in group_samples if row[s] == 1]
            if active:
                ys = [y_map[s] for s in active]
                ax_mat.scatter([i] * len(ys), ys, s=130, color=color, edgecolor="black", linewidth=0.8, zorder=3)
                if len(ys) > 1:
                    ax_mat.plot([i, i], [min(ys), max(ys)], color=color, linewidth=1.6, zorder=2)

    ax_mat.set_yticks(list(y_map.values()))
    ax_mat.set_yticklabels(all_samples)
    ax_mat.invert_yaxis()
    ax_mat.set_xlabel("Intersection combinations")
    for side in ("top", "right", "bottom"):
        ax_mat.spines[side].set_visible(False)
    ax_mat.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_mat.axhline(len(group_a_samples) - 0.5, color="black", linewidth=0.8, linestyle="--", alpha=0.7)

    colors = [GROUP_A_COLOR] * len(group_a_samples) + [GROUP_B_COLOR] * len(group_b_samples)
    ax_set.barh(list(y_map.values()), [set_sizes[s] for s in all_samples], color=colors, edgecolor="black")
    ax_set.set_yticks(list(y_map.values()))
    ax_set.set_yticklabels(all_samples)
    ax_set.invert_yaxis()
    ax_set.set_xlabel("Set size")
    ax_set.invert_xaxis()
    ax_set.spines["top"].set_visible(False)
    ax_set.spines["right"].set_visible(False)
    ax_set.set_ylim(ax_mat.get_ylim())

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--intersection-table", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("results/upset_plot.png"))
    parser.add_argument("--group-a-size-column", default="group_a_size")
    parser.add_argument("--group-b-size-column", default="group_b_size")
    parser.add_argument("--group-a-prefix", default="A_")
    parser.add_argument("--group-b-prefix", default="B_")
    parser.add_argument("--group-a-label", default="Group A")
    parser.add_argument("--group-b-label", default="Group B")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    df, group_a_samples, group_b_samples = load_intersection_table(
        args.intersection_table, args.group_a_size_column, args.group_b_size_column,
        args.group_a_prefix, args.group_b_prefix,
    )
    plot_upset(
        df, args.group_a_size_column, args.group_b_size_column,
        group_a_samples, group_b_samples, args.group_a_label, args.group_b_label,
        args.output,
    )
    print(f"Wrote UpSet-style plot to: {args.output}")


if __name__ == "__main__":
    main()
