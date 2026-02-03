#!/bin/bash

EXP_NAME="exp-01-reproduce-results"
RUN_NAME="run1-reproduce-results"

PROJECT_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

RUN_DIR="$PROJECT_ROOT/exps/${EXP_NAME}/outputs/${RUN_NAME}"
mkdir -p "$RUN_DIR"

python -u -m fitness_landscape_graph.build_graph \
    --base-path "$PROJECT_ROOT" \
    --output-dir "$RUN_DIR" \
    --neutral-threshold 0.4 \
    --tiny-initial-threshold 0.02 \
    --large-edge-threshold 5.5 \
    --num-forbidden-pairs 1 \
    > "$RUN_DIR/build_graph.log" 2>&1
