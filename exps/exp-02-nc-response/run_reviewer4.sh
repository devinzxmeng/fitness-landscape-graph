#!/bin/bash
PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../.." &> /dev/null && pwd )"
OUTPUT_DIR="$PROJECT_ROOT/exps/exp-02-nc-response/outputs"
mkdir -p "$OUTPUT_DIR"

python -u "$PROJECT_ROOT/exps/exp-02-nc-response/run_reviewer4.py" \
    --pairs-path "$PROJECT_ROOT/data/processed/azt_pairs.csv" \
    --output-dir "$OUTPUT_DIR" \
    --target-genotype "P.LKN...K...." \
    --max-distance 2 \
    2>&1 | tee "$OUTPUT_DIR/reviewer4.log"
