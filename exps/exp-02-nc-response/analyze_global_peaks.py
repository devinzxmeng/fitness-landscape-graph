"""Analyze global peak existence across threshold/concentration combinations.

Loads all Aztreonam graphs from the global-peak-robustness analysis and
checks each for the presence of a global peak (highest-fitness node that
is also a large peak).

Output:
    CSV file with columns:
    - concentration: Drug concentration (12.0, 36.0, 108.0, 324.0)
    - neutral_threshold: Neutral threshold value (0.15 to 0.45)
    - has_global_peak: Boolean indicator
    - peak_node: Node ID (if exists)
    - peak_fitness: Fitness value (if exists)
    - peak_group_size: Group size (if exists)
"""

import sys
from pathlib import Path

import polars as pl

from fitness_landscape_graph.graph_analyzer import GraphAnalyzer


def analyze_all_graphs(
    output_dir: Path,
    min_group_size: int = 32,
) -> pl.DataFrame:
    """Analyze all graphs in the output directory for global peaks.

    Args:
        output_dir: Directory containing graph files (azt_c*_t*.graphml).
        min_group_size: Minimum group size for a peak to be considered global.

    Returns:
        Polars DataFrame with analysis results.
    """
    # Find all Aztreonam graph files
    graph_files = sorted(output_dir.glob("azt_c*_t*.graphml"))

    if not graph_files:
        print(f"❌ No graph files found in {output_dir}")
        print("   Looking for pattern: azt_c*_t*.graphml")
        return pl.DataFrame()

    print(f"Found {len(graph_files)} graph files to analyze")
    print(f"Using min_group_size = {min_group_size}")
    print("")

    results = []

    for i, graph_path in enumerate(graph_files, 1):
        # Parse filename: azt_c{conc}_t{threshold}.graphml
        # Dots are replaced with underscores in the filename
        # Example: azt_c12_0_t0_15.graphml -> concentration=12.0, threshold=0.15
        filename = graph_path.stem  # Remove .graphml

        try:
            # Split on "_t" to separate concentration and threshold parts
            # azt_c12_0_t0_15 -> ["azt_c12_0", "0_15"]
            if "_t" not in filename:
                print(
                    f"⚠ Warning: Could not parse filename {filename}: no '_t' separator"
                )
                continue

            conc_part, thresh_part = filename.rsplit("_t", 1)

            # Extract concentration: "azt_c12_0" -> "c12_0" -> "12_0" -> 12.0
            conc_str = conc_part.split("_", 1)[1]  # Remove "azt_" prefix -> "c12_0"
            conc_str = conc_str[1:]  # Remove "c" prefix -> "12_0"
            concentration = float(conc_str.replace("_", "."))  # 12.0

            # Extract threshold: "0_15" -> 0.15
            neutral_threshold = float(thresh_part.replace("_", "."))

        except (IndexError, ValueError) as e:
            print(f"⚠ Warning: Could not parse filename {filename}: {e}")
            continue

        # Load graph and analyze
        try:
            analyzer = GraphAnalyzer(str(graph_path))
            has_peak, node_id, info = analyzer.has_global_peak(
                min_group_size=min_group_size
            )

            result = {
                "concentration": concentration,
                "neutral_threshold": neutral_threshold,
                "has_global_peak": has_peak,
                "peak_node": node_id if has_peak else None,
                "peak_fitness": info.get("fitness") if has_peak else None,
                "peak_group_size": info.get("group_size") if has_peak else None,
            }

            results.append(result)

            # Progress indicator
            if i % 10 == 0 or i == len(graph_files):
                print(f"  Processed {i}/{len(graph_files)} graphs...")

        except Exception as e:
            print(f"❌ Error analyzing {graph_path.name}: {e}")
            continue

    if not results:
        print("❌ No results collected")
        return pl.DataFrame()

    # Create DataFrame
    df = pl.DataFrame(results)

    # Sort by concentration, then threshold
    df = df.sort(["concentration", "neutral_threshold"])

    return df


def print_summary(df: pl.DataFrame) -> None:
    """Print summary statistics of the analysis.

    Args:
        df: Results DataFrame from analyze_all_graphs.
    """
    print("\n" + "=" * 80)
    print("ANALYSIS SUMMARY")
    print("=" * 80)

    total_graphs = len(df)
    total_with_peak = df["has_global_peak"].sum()

    print(f"\nTotal graphs analyzed: {total_graphs}")
    print(
        f"Graphs with global peak: {total_with_peak} "
        f"({100 * total_with_peak / total_graphs:.1f}%)"
    )

    # Summary by concentration
    print("\nBy concentration:")
    conc_summary = (
        df.group_by("concentration")
        .agg(
            [
                pl.col("has_global_peak").sum().alias("count_with_peak"),
                pl.col("has_global_peak").len().alias("total"),
            ]
        )
        .sort("concentration")
    )

    for row in conc_summary.iter_rows(named=True):
        conc = row["concentration"]
        count = row["count_with_peak"]
        total = row["total"]
        pct = 100 * count / total if total > 0 else 0
        print(f"  {conc:6.1f}: {count:2d}/{total:2d} ({pct:5.1f}%)")

    # Summary by threshold range
    print("\nBy threshold range:")
    threshold_ranges = [
        (0.15, 0.24, "0.15-0.24"),
        (0.25, 0.34, "0.25-0.34"),
        (0.35, 0.45, "0.35-0.45"),
    ]

    for t_min, t_max, label in threshold_ranges:
        subset = df.filter(
            (pl.col("neutral_threshold") >= t_min)
            & (pl.col("neutral_threshold") <= t_max)
        )
        count = subset["has_global_peak"].sum()
        total = len(subset)
        pct = 100 * count / total if total > 0 else 0
        print(f"  {label}: {count:2d}/{total:2d} ({pct:5.1f}%)")

    # Highest fitness peaks
    peaks_df = df.filter(pl.col("has_global_peak"))
    if len(peaks_df) > 0:
        print("\nTop 5 highest fitness peaks:")
        top_peaks = peaks_df.sort("peak_fitness", descending=True).head(5)
        for row in top_peaks.iter_rows(named=True):
            print(
                f"  c={row['concentration']:6.1f}, t={row['neutral_threshold']:.2f}: "
                f"fitness={row['peak_fitness']:.4f}, size={row['peak_group_size']}"
            )


def main():
    """Main analysis workflow."""
    # Determine output directory
    script_dir = Path(__file__).parent
    output_dir = script_dir / "outputs" / "global-peak-robustness"

    if not output_dir.exists():
        print(f"❌ Output directory not found: {output_dir}")
        print("   Please run build_all_threshold_graphs.sh first")
        sys.exit(1)

    print("=" * 80)
    print("GLOBAL PEAK ROBUSTNESS ANALYSIS")
    print("=" * 80)
    print(f"Input directory: {output_dir}")
    print("")

    # Run analysis
    df = analyze_all_graphs(output_dir, min_group_size=32)

    if len(df) == 0:
        print("❌ No data collected. Exiting.")
        sys.exit(1)

    # Print summary
    print_summary(df)

    # Save results
    output_csv = output_dir / "global_peak_analysis.csv"
    df.write_csv(output_csv)

    print("\n" + "=" * 80)
    print(f"✓ Results saved to: {output_csv}")
    print("=" * 80)


if __name__ == "__main__":
    main()
