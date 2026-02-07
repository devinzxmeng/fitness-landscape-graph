#!/bin/bash
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." &> /dev/null && pwd )"
OUTPUT_DIR="$PROJECT_ROOT/exps/exp-02-nc-response/outputs/fitness-advantage-azt-12"
mkdir -p "$OUTPUT_DIR"

python -u "$PROJECT_ROOT/exps/exp-02-nc-response/analyze_peak_fitness_advantage.py" \
    --pairs-path "$PROJECT_ROOT/data/processed/azt_pairs.csv" \
    --graph-path "$PROJECT_ROOT/exps/exp-01-reproduce-results/outputs/run1-reproduce-results/azt_12_0.graphml" \
    --output-dir "$OUTPUT_DIR" \
    --max-distance 2 \
    --peak-rank 0 \
    --antibiotic aztreonam \
    2>&1 | tee "$OUTPUT_DIR/peak_fitness_advantage.log"
