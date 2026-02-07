"""Analyze fitness advantage of genotype groups across drug concentrations.

Takes a group of genotypes (from peak supernode, custom node selection, or file)
and computes their fitness advantage over external 1-mutation and 2-mutation
neighbors at each drug concentration. Produces a box plot showing how the
group's fitness advantage changes with concentration.

Supports three input modes:
1. --node-ids: Analyze genotypes from specific graph nodes (comma-separated)
2. --genotypes-file: Analyze genotypes from a text file (one per line)
3. --peak-rank: Convenience shortcut to analyze the Nth peak (0=highest fitness)
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from fitness_landscape_graph.fitness_advantage import FitnessAdvantageAnalyzer
from fitness_landscape_graph.graph_analyzer import GraphAnalyzer


def make_boxplot(
    result_df: pl.DataFrame,
    output_dir: str,
    group_size: int,
    antibiotic: str = "aztreonam",
) -> None:
    """Create box plot of fitness advantage vs concentration.

    Args:
        result_df: Output of FitnessAdvantageAnalyzer.compute_fitness_advantage.
        output_dir: Directory to save figures.
        group_size: Number of genotypes in the peak group.
        antibiotic: Antibiotic name for plot labels.
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
    ax.set_xlabel(f"{antibiotic.capitalize()} concentration (µg/mL)")
    ax.set_ylabel("Fitness advantage (group member − external neighbor)")
    ax.set_title(
        f"Fitness advantage of genotype group ({group_size} genotypes) "
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
        description="Analyze fitness advantage of genotype groups across concentrations"
    )

    # Mutually exclusive group for genotype selection
    selection_group = parser.add_mutually_exclusive_group(required=True)

    selection_group.add_argument(
        "--node-ids",
        type=str,
        help="Comma-separated list of node IDs to analyze (e.g., 'node1,node2,node3')",
    )

    selection_group.add_argument(
        "--genotypes-file",
        type=str,
        help="Path to file with one genotype per line",
    )

    selection_group.add_argument(
        "--peak-rank",
        type=int,
        help="Convenience: analyze the Nth peak (0=highest fitness, 1=second, etc.)",
    )

    # Required arguments
    parser.add_argument(
        "--pairs-path",
        type=str,
        required=True,
        help="Path to pairwise table CSV",
    )
    parser.add_argument(
        "--graph-path",
        type=str,
        required=True,
        help="Path to GraphML file",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory to save output files",
    )

    # Optional arguments
    parser.add_argument(
        "--max-distance",
        type=int,
        default=2,
        help="Maximum mutation distance to consider (1 or 2)",
    )
    parser.add_argument(
        "--antibiotic",
        type=str,
        default="aztreonam",
        help="Antibiotic name for plot labels",
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load graph
    print(f"Loading graph from {args.graph_path}")
    graph_analyzer = GraphAnalyzer(args.graph_path)

    # Resolve genotypes based on selection method
    if args.node_ids:
        node_list = [n.strip() for n in args.node_ids.split(",")]
        group_genotypes = graph_analyzer.get_nodes_genotypes(node_list)
        print(f"Selected {len(node_list)} nodes: {node_list}")
    elif args.genotypes_file:
        with open(args.genotypes_file) as f:
            group_genotypes = set(line.strip() for line in f if line.strip())
        print(f"Loaded {len(group_genotypes)} genotypes from file")
    elif args.peak_rank is not None:
        group_genotypes = graph_analyzer.get_peak_genotypes(rank=args.peak_rank)
        print(f"Selected peak rank {args.peak_rank}")

    print(f"Analyzing {len(group_genotypes)} genotypes")

    if not group_genotypes:
        print("ERROR: No genotypes found. Exiting.")
        sys.exit(1)

    # Load pair table
    print(f"Loading pair table from {args.pairs_path}")
    pairs_df = pl.read_csv(args.pairs_path)
    concentrations = sorted(pairs_df["concentration"].unique().to_list())
    print(f"Concentrations: {concentrations}")

    # Create fitness advantage analyzer
    fitness_analyzer = FitnessAdvantageAnalyzer(pairs_df, graph_analyzer)

    # Compute fitness advantages
    print("Computing fitness advantages...")
    result_df = fitness_analyzer.compute_fitness_advantage(
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
    make_boxplot(result_df, args.output_dir, len(group_genotypes), args.antibiotic)


if __name__ == "__main__":
    main()
