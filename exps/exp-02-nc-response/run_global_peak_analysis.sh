#!/bin/bash

# Global peak robustness analysis workflow
# Deliverable: Robustness/sensitivity analysis showing global peak existence
# across different neutrality cutoff values and drug concentrations.
#
# Steps:
# 1. Build 124 graphs (31 thresholds × 4 concentrations) in parallel
# 2. Analyze graphs for global peak existence
# 3. Create heatmap and summary visualizations

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
echo "Step 2: Analyzing graphs for global peak existence..."
echo "------------------------------------------------------------------------"
python "$SCRIPT_DIR/analyze_global_peaks.py"

echo ""
echo "Step 3: Creating visualizations..."
echo "------------------------------------------------------------------------"
python "$SCRIPT_DIR/plot_global_peak_heatmap.py"

echo ""
echo "========================================================================"
echo "WORKFLOW COMPLETE"
echo "========================================================================"
echo ""
echo "Results in: $OUTPUT_DIR/"
echo "  - global_peak_analysis.csv"
echo "  - global_peak_heatmap.png"
echo "  - global_peak_by_concentration.png"
echo "  - global_peak_by_threshold.png"
