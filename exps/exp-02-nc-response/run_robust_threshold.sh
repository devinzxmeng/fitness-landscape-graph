#!/bin/bash

# Robust neutral threshold analysis script
# Tests multiple neutral threshold values to assess graph robustness

EXP_NAME="exp-02-nc-response"
RUN_NAME="robust-neutral-threshold"

PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." &> /dev/null && pwd )"

# Define neutral thresholds to test
THRESHOLDS=(0.08 0.14 0.22 0.27 0.40 0.45)

BASE_OUTPUT_DIR="$PROJECT_ROOT/exps/${EXP_NAME}/outputs/${RUN_NAME}"
mkdir -p "$BASE_OUTPUT_DIR"

echo "Starting robust neutral threshold analysis..."
echo "Output directory: $BASE_OUTPUT_DIR"
echo "Testing ${#THRESHOLDS[@]} threshold values"
echo ""

# Fixed parameters
TINY_INITIAL_THRESHOLD=0.02
LARGE_EDGE_THRESHOLD=5.5
NUM_FORBIDDEN_PAIRS=1

# Loop through each threshold
for threshold in "${THRESHOLDS[@]}"; do
    # Convert threshold to folder-safe format (0.08 -> 0-08)
    folder_suffix=$(echo "$threshold" | tr '.' '-')
    RUN_DIR="$BASE_OUTPUT_DIR/neutral-threshold-$folder_suffix"
    
    echo "Running with neutral-threshold = $threshold"
    echo "  Output: $RUN_DIR"
    
    mkdir -p "$RUN_DIR"
    
    python -u -m fitness_landscape_graph.build_graph \
        --base-path "$PROJECT_ROOT" \
        --output-dir "$RUN_DIR" \
        --neutral-threshold "$threshold" \
        --tiny-initial-threshold "$TINY_INITIAL_THRESHOLD" \
        --large-edge-threshold "$LARGE_EDGE_THRESHOLD" \
        --num-forbidden-pairs "$NUM_FORBIDDEN_PAIRS" \
        > "$RUN_DIR/build_graph.log" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Completed successfully"
    else
        echo "  ✗ Failed (check $RUN_DIR/build_graph.log)"
    fi
    echo ""
done

echo "All threshold runs completed!"
echo "Results saved in: $BASE_OUTPUT_DIR"
