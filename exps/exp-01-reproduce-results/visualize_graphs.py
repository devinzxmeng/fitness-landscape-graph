"""Generate PNG visualizations from GraphML files using GraphAnalyzer.

This script loads all GraphML files from the run1-reproduce-results directory
and generates high-quality PNG images using the GraphAnalyzer.plot() function.
"""

from pathlib import Path

from fitness_landscape_graph.graph_analyzer import GraphAnalyzer, VisConfig


def main():
    """Generate PNG visualizations for all GraphML files."""
    # Define paths
    graphml_dir = Path("outputs/run1-reproduce-results")
    output_dir = graphml_dir / "images"

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all GraphML files
    graphml_files = sorted(graphml_dir.glob("*.graphml"))

    if not graphml_files:
        print(f"No GraphML files found in {graphml_dir}")
        return

    print(f"Found {len(graphml_files)} GraphML files")
    print(f"Output directory: {output_dir}")
    print("=" * 80)

    # Configure visualization
    config = VisConfig(
        figure_width=1200,
        figure_height=800,
    )

    # Process each file
    for i, graphml_path in enumerate(graphml_files, 1):
        print(f"[{i}/{len(graphml_files)}] Processing {graphml_path.name}...", end=" ")

        try:
            # Load graph and generate figure
            analyzer = GraphAnalyzer(str(graphml_path))
            fig = analyzer.plot(config=config)

            # Save as PNG
            output_path = output_dir / f"{graphml_path.stem}.png"
            fig.write_image(str(output_path), format="png", scale=2)

            print(f"✓ Saved to {output_path.name}")

        except Exception as e:
            print(f"✗ Error: {e}")

    print("=" * 80)
    print(f"✓ Complete! Images saved to {output_dir}")


if __name__ == "__main__":
    main()
