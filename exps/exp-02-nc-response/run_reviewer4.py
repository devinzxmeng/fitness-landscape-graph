"""Reviewer #4 analysis: fitness advantage of aztreonam peak across concentrations.

Produces a box plot showing how the fitness advantage of the global peak
genotype over its 1-mutation and 2-mutation neighbors changes with drug
concentration. This addresses the reviewer's request for a quantitative
diagnostic of the peak "absorption" at intermediate concentrations.
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import polars as pl

# Add project root to path so we can import the helper module
sys.path.insert(0, os.path.dirname(__file__))
from graph_analysis import compute_fitness_advantage_across_concentrations


def make_boxplot(
    result_df: pl.DataFrame,
    target: str,
    output_dir: str,
) -> None:
    """Create box plot of fitness advantage vs concentration.

    Args:
        result_df: Output of compute_fitness_advantage_across_concentrations.
        target: Target genotype string (for title).
        output_dir: Directory to save figures.
    """
    concentrations = sorted(result_df["concentration"].unique().to_list())
    distances = sorted(result_df["distance"].unique().to_list())

    fig, ax = plt.subplots(figsize=(10, 5))

    # Positions for grouped box plots
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

        label = f"{dist}-mutation" if dist not in labels_added else None
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

        # Legend entry
        if dist not in labels_added:
            bp["boxes"][0].set_label(f"{dist}-mutation neighbors")
            labels_added.add(dist)

    ax.axhline(y=0, color="grey", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.set_xticks(x_base)
    ax.set_xticklabels([str(c) for c in concentrations])
    ax.set_xlabel("Aztreonam concentration (µg/mL)")
    ax.set_ylabel("Fitness advantage (target − neighbor)")
    ax.set_title(f"Fitness advantage of peak {target}")
    ax.legend(loc="upper left")

    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, "fitness_advantage_boxplot.pdf"),
        bbox_inches="tight",
        dpi=300,
    )
    plt.savefig(
        os.path.join(output_dir, "fitness_advantage_boxplot.png"),
        bbox_inches="tight",
        dpi=300,
    )
    plt.close()
    print(f"Saved box plot to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Reviewer #4: aztreonam peak fitness advantage analysis"
    )
    parser.add_argument(
        "--pairs-path",
        type=str,
        default="data/processed/azt_pairs.csv",
        help="Path to aztreonam pairwise table CSV",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="exps/exp-02-nc-response/outputs",
        help="Directory to save output files",
    )
    parser.add_argument(
        "--target-genotype",
        type=str,
        default="P.LKN...K....",
        help="13-char mutant profile of the peak to track",
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=2,
        help="Maximum mutation distance to consider (1 or 2)",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Loading pair table from {args.pairs_path}")
    pairs_df = pl.read_csv(args.pairs_path)
    concentrations = sorted(pairs_df["concentration"].unique().to_list())
    print(f"Concentrations: {concentrations}")
    print(f"Target genotype: {args.target_genotype}")

    print("Computing fitness advantages...")
    result_df = compute_fitness_advantage_across_concentrations(
        pairs_df=pairs_df,
        target=args.target_genotype,
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
                    f"n={len(diffs):>3d}  "
                    f"median_adv={np.median(diffs):>+.3f}  "
                    f"mean_adv={np.mean(diffs):>+.3f}"
                )

    # Make box plot
    make_boxplot(result_df, args.target_genotype, args.output_dir)


if __name__ == "__main__":
    main()
