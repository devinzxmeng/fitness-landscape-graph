"""Visualize global peak existence as a heatmap.

Creates a heatmap showing which threshold/concentration combinations
have a global peak (highest-fitness node that is also a large peak).

Input:
    CSV file from analyze_global_peaks.py

Output:
    PNG heatmap with:
    - X-axis: Neutral threshold (0.15 to 0.45, raw values)
    - Y-axis: Concentration (12.0, 36.0, 108.0, 324.0, log scale)
    - Color: Binary indicator (has global peak: True/False)
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import polars as pl


def plot_heatmap(df: pl.DataFrame, output_path: Path) -> None:
    """Create heatmap visualization of global peak existence.

    Args:
        df: Results DataFrame with columns: concentration, neutral_threshold,
            has_global_peak.
        output_path: Path to save the output PNG file.
    """
    # Get unique values
    concentrations = sorted(df["concentration"].unique().to_list())
    thresholds = sorted(df["neutral_threshold"].unique().to_list())

    print(f"Concentrations: {concentrations}")
    print(
        f"Thresholds: {len(thresholds)} values from {min(thresholds):.2f} "
        f"to {max(thresholds):.2f}"
    )

    # Create 2D array for heatmap [concentrations × thresholds]
    heatmap_data = np.zeros((len(concentrations), len(thresholds)))

    for row in df.iter_rows(named=True):
        conc = row["concentration"]
        thresh = row["neutral_threshold"]
        has_peak = row["has_global_peak"]

        # Find indices
        try:
            i = concentrations.index(conc)
            j = thresholds.index(thresh)
            heatmap_data[i, j] = 1 if has_peak else 0
        except ValueError:
            print(f"⚠ Warning: Could not find index for c={conc}, t={thresh}")
            continue

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 6))

    # Plot heatmap
    im = ax.imshow(
        heatmap_data,
        aspect="auto",
        cmap="RdYlGn",  # Red (no peak) to Green (has peak)
        vmin=0,
        vmax=1,
        interpolation="nearest",
    )

    # Set ticks and labels
    ax.set_xticks(np.arange(len(thresholds)))
    ax.set_yticks(np.arange(len(concentrations)))

    # X-axis: show every 3rd threshold to avoid crowding
    x_labels = [f"{t:.2f}" if i % 3 == 0 else "" for i, t in enumerate(thresholds)]
    ax.set_xticklabels(x_labels, rotation=45, ha="right")

    # Y-axis: show all concentrations with log scale indication
    y_labels = [f"{c:.1f}" for c in concentrations]
    ax.set_yticklabels(y_labels)

    # Labels
    ax.set_xlabel("Neutral Threshold", fontsize=12, fontweight="bold")
    ax.set_ylabel("Concentration (µg/mL)", fontsize=12, fontweight="bold")
    ax.set_title(
        "Global Peak Existence in Aztreonam Fitness Landscapes\n"
        "(Global peak = highest-fitness node with group_size > 32)",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, ticks=[0, 1])
    cbar.ax.set_yticklabels(["No", "Yes"])
    cbar.set_label("Has Global Peak", rotation=270, labelpad=20, fontweight="bold")

    # Add grid for better readability
    ax.set_xticks(np.arange(len(thresholds)) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(concentrations)) - 0.5, minor=True)
    ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.5, alpha=0.3)

    # Add text annotations in each cell
    for i in range(len(concentrations)):
        for j in range(len(thresholds)):
            text = "✓" if heatmap_data[i, j] == 1 else "✗"
            color = "white" if heatmap_data[i, j] == 1 else "black"
            ax.text(
                j,
                i,
                text,
                ha="center",
                va="center",
                color=color,
                fontsize=8,
                fontweight="bold",
            )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✓ Heatmap saved to: {output_path}")

    plt.close()


def plot_summary_by_concentration(df: pl.DataFrame, output_path: Path) -> None:
    """Create line plot showing fraction of thresholds with global peak.

    Args:
        df: Results DataFrame.
        output_path: Path to save the output PNG file.
    """
    # Calculate fraction of thresholds with global peak for each concentration
    summary = (
        df.group_by("concentration")
        .agg(
            [
                pl.col("has_global_peak").mean().alias("fraction_with_peak"),
                pl.col("has_global_peak").sum().alias("count_with_peak"),
                pl.col("has_global_peak").len().alias("total"),
            ]
        )
        .sort("concentration")
    )

    concentrations = summary["concentration"].to_list()
    fractions = summary["fraction_with_peak"].to_list()
    counts = summary["count_with_peak"].to_list()
    totals = summary["total"].to_list()

    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot line
    ax.plot(
        concentrations,
        fractions,
        marker="o",
        markersize=10,
        linewidth=2,
        color="#2E86AB",
    )

    # Add value labels
    for c, f, count, total in zip(concentrations, fractions, counts, totals):
        ax.annotate(
            f"{count}/{total}\n({f * 100:.0f}%)",
            xy=(c, f),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            fontsize=10,
            fontweight="bold",
        )

    ax.set_xlabel("Concentration (µg/mL)", fontsize=12, fontweight="bold")
    ax.set_ylabel(
        "Fraction of Thresholds\nwith Global Peak", fontsize=12, fontweight="bold"
    )
    ax.set_title(
        "Global Peak Prevalence Across Concentrations\n(Aztreonam)",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)
    ax.set_xscale("log")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✓ Summary plot saved to: {output_path}")

    plt.close()


def plot_summary_by_threshold(df: pl.DataFrame, output_path: Path) -> None:
    """Create line plot showing fraction of concentrations with global peak.

    Args:
        df: Results DataFrame.
        output_path: Path to save the output PNG file.
    """
    # Calculate fraction of concentrations with global peak for each threshold
    summary = (
        df.group_by("neutral_threshold")
        .agg(
            [
                pl.col("has_global_peak").mean().alias("fraction_with_peak"),
                pl.col("has_global_peak").sum().alias("count_with_peak"),
                pl.col("has_global_peak").len().alias("total"),
            ]
        )
        .sort("neutral_threshold")
    )

    thresholds = summary["neutral_threshold"].to_list()
    fractions = summary["fraction_with_peak"].to_list()

    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot line
    ax.plot(
        thresholds, fractions, marker="o", markersize=6, linewidth=2, color="#A23B72"
    )

    # Add shaded regions for threshold ranges
    ax.axvspan(0.15, 0.24, alpha=0.1, color="red", label="Low (0.15-0.24)")
    ax.axvspan(0.25, 0.34, alpha=0.1, color="yellow", label="Medium (0.25-0.34)")
    ax.axvspan(0.35, 0.45, alpha=0.1, color="green", label="High (0.35-0.45)")

    ax.set_xlabel("Neutral Threshold", fontsize=12, fontweight="bold")
    ax.set_ylabel(
        "Fraction of Concentrations\nwith Global Peak", fontsize=12, fontweight="bold"
    )
    ax.set_title(
        "Global Peak Prevalence Across Thresholds\n(Aztreonam)",
        fontsize=14,
        fontweight="bold",
        pad=20,
    )

    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(0.14, 0.46)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"✓ Threshold summary plot saved to: {output_path}")

    plt.close()


def main():
    """Main visualization workflow."""
    # Determine input/output paths
    script_dir = Path(__file__).parent
    results_dir = script_dir / "outputs" / "global-peak-robustness"
    input_csv = results_dir / "global_peak_analysis.csv"

    if not input_csv.exists():
        print(f"❌ Input CSV not found: {input_csv}")
        print("   Please run analyze_global_peaks.py first")
        sys.exit(1)

    print("=" * 80)
    print("GLOBAL PEAK HEATMAP VISUALIZATION")
    print("=" * 80)
    print(f"Input: {input_csv}")
    print("")

    # Load data
    df = pl.read_csv(input_csv)
    print(f"Loaded {len(df)} records")
    print("")

    # Create visualizations
    heatmap_path = results_dir / "global_peak_heatmap.png"
    plot_heatmap(df, heatmap_path)

    summary_conc_path = results_dir / "global_peak_by_concentration.png"
    plot_summary_by_concentration(df, summary_conc_path)

    summary_thresh_path = results_dir / "global_peak_by_threshold.png"
    plot_summary_by_threshold(df, summary_thresh_path)

    print("\n" + "=" * 80)
    print("✓ All visualizations complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
