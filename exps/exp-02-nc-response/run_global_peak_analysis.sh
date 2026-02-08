#!/bin/bash

# Global peak robustness analysis — graph building step
# Builds 124 graphs (31 thresholds × 4 concentrations) in parallel.
# After building, open and run neutral_threshold_robustness.ipynb for
# analysis and visualization.

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$( cd "$SCRIPT_DIR/../.." &> /dev/null && pwd )"
OUTPUT_DIR="$SCRIPT_DIR/outputs/global-peak-robustness"

echo "========================================================================"
echo "GLOBAL PEAK ROBUSTNESS ANALYSIS"
echo "========================================================================"
echo ""

# Step 1: Build all graphs in parallel
echo "Step 1: Building 124 graphs (31 thresholds × 4 concentrations)..."
echo "------------------------------------------------------------------------"

python -m fitness_landscape_graph.build_graphs_parallel \
    --base-path "$PROJECT_ROOT" \
    --output-dir "$OUTPUT_DIR" \
    --neutral-thresholds 0.15 0.45 0.01 \
    --concentrations 12.0 36.0 108.0 324.0

echo ""
echo "========================================================================"
echo "GRAPH BUILDING COMPLETE"
echo "========================================================================"
echo ""
echo "Graphs saved to: $OUTPUT_DIR/"
echo ""
echo "Next: open and run neutral_threshold_robustness.ipynb for analysis + visualization."
