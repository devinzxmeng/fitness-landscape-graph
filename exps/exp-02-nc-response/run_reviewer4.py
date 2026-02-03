"""Reviewer #4 analysis: fitness advantage of aztreonam peak group across concentrations.

Takes all genotypes in the global peak supernode (from the azt_108 graph) and
computes their fitness advantage over external 1-mutation and 2-mutation
neighbors at each drug concentration. Produces a box plot showing how the
peak's fitness advantage changes with concentration.
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

sys.path.insert(0, os.path.dirname(__file__))
from graph_analysis import (
    compute_group_fitness_advantage,
    get_peak_genotypes,
    load_graph,
)


def make_boxplot(
    result_df: pl.DataFrame,
    output_dir: str,
    group_size: int,
) -> None:
    """Create box plot of fitness advantage vs concentration.

    Args:
        result_df: Output of compute_group_fitness_advantage.
        output_dir: Directory to save figures.
        group_size: Number of genotypes in the peak group.
    """
    concentrations = sorted(result_df["concentration"].unique().to_list())
    distances = sorted(result_df["distance"].unique().to_list())

    fig, ax = plt.subplots(figsize=(10, 5))

    n_conc = len(concentrations)
    n_dist = len(distances)
    width = 0.35
    x_base = np.arange(n_conc)

    colors = {1: "#4C72B0", 2: "#DD8452"}
    labels_added = set()

    for i, dist in enumerate(distances):
        dist_df = result_df.filter(pl.col("distance") == dist)
        box_data = []
        positions = []
        for j, conc in enumerate(concentrations):
            conc_data = dist_df.filter(pl.col("concentration") == conc)[
                "fitness_diff"
            ].to_list()
            if conc_data:
                box_data.append(conc_data)
                positions.append(x_base[j] + (i - (n_dist - 1) / 2) * width)

        if not box_data:
            continue

        bp = ax.boxplot(
            box_data,
            positions=positions,
            widths=width * 0.8,
            patch_artist=True,
            showfliers=True,
            flierprops={"markersize": 3, "alpha": 0.5},
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(colors.get(dist, "grey"))
            patch.set_alpha(0.7)
        for median_line in bp["medians"]:
            median_line.set_color("black")

        if dist not in labels_added:
            bp["boxes"][0].set_label(f"{dist}-mutation neighbors")
            labels_added.add(dist)

    ax.axhline(y=0, color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_xticks(x_base)
    ax.set_xticklabels([str(c) for c in concentrations])
    ax.set_xlabel("Aztreonam concentration (µg/mL)")
    ax.set_ylabel("Fitness advantage (group member − external neighbor)")
    ax.set_title(
        f"Fitness advantage of azt_108 global peak group ({group_size} genotypes) "
        f"over external neighbors"
    )
    ax.legend(loc="upper left")

    plt.tight_layout()
    for ext in ("pdf", "png"):
        path = os.path.join(output_dir, f"fitness_advantage_boxplot.{ext}")
        plt.savefig(path, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved box plot to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Reviewer #4: aztreonam peak group fitness advantage analysis"
    )
    parser.add_argument(
        "--pairs-path",
        type=str,
        default="data/processed/azt_pairs.csv",
        help="Path to aztreonam pairwise table CSV",
    )
    parser.add_argument(
        "--graph-path",
        type=str,
        default="exps/exp-01-reproduce-results/outputs/run1-reproduce-results/azt_108_0.graphml",
        help="Path to GraphML file to extract the peak group from",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="exps/exp-02-nc-response/outputs",
        help="Directory to save output files",
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=2,
        help="Maximum mutation distance to consider (1 or 2)",
    )
    parser.add_argument(
        "--peak-rank",
        type=int,
        default=0,
        help="Which peak to analyze (0 = highest fitness, 1 = second, etc.)",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load peak group from graph
    print(f"Loading graph from {args.graph_path}")
    graph = load_graph(args.graph_path)
    group_genotypes = get_peak_genotypes(graph, rank=args.peak_rank)
    print(f"Peak group has {len(group_genotypes)} genotypes")

    # Load pair table
    print(f"Loading pair table from {args.pairs_path}")
    pairs_df = pl.read_csv(args.pairs_path)
    concentrations = sorted(pairs_df["concentration"].unique().to_list())
    print(f"Concentrations: {concentrations}")

    # Compute fitness advantages
    print("Computing fitness advantages...")
    result_df = compute_group_fitness_advantage(
        pairs_df=pairs_df,
        group_genotypes=group_genotypes,
        concentrations=concentrations,
        max_distance=args.max_distance,
    )

    # Save CSV
    csv_path = os.path.join(args.output_dir, "fitness_advantage_data.csv")
    result_df.write_csv(csv_path)
    print(f"Saved data to {csv_path} ({result_df.height} rows)")

    # Print summary
    for conc in concentrations:
        conc_data = result_df.filter(pl.col("concentration") == conc)
        for dist in sorted(conc_data["distance"].unique().to_list()):
            dist_data = conc_data.filter(pl.col("distance") == dist)
            diffs = dist_data["fitness_diff"].to_list()
            if diffs:
                print(
                    f"  conc={conc:>6.1f}  dist={dist}  "
                    f"n={len(diffs):>5d}  "
                    f"median_adv={np.median(diffs):>+.3f}  "
                    f"mean_adv={np.mean(diffs):>+.3f}"
                )

    # Make box plot
    make_boxplot(result_df, args.output_dir, len(group_genotypes))


if __name__ == "__main__":
    main()
