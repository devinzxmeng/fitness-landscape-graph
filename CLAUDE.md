# CLAUDE.md

## Project Overview

Fitness Landscape Graph: builds and visualizes fitness landscape graphs for antibiotic resistance mutations in TEM-1 beta-lactamase genotypes. Analyzes mutation networks for two antibiotics (Ampicillin, Aztreonam) across multiple drug concentrations.

## Tech Stack

- Python 3.11+
- Polars (DataFrames), NetworkX (graphs), NumPy, SciPy (stats), Logomaker, Matplotlib
- Ruff for linting/formatting
- Gephi (external) for visualizing graphs

## Repository Structure

```
src/fitness_landscape_graph/
  build_graph.py       # Main entry point & CLI argument parsing
  graph_builder.py     # Core GraphBuilder class (network construction & simplification)
  preprocess.py        # Data cleaning & filtering pipeline
  pair_table.py        # Concentration-specific pairwise mutation analysis
  pair_table_global.py # Global fitness (AUC) calculations & pair generation
  make_logo.py         # Sequence logo visualization utilities
data/
  raw/combined-auc/    # Raw genotype AUC CSV files (~20-30MB each)
  processed/           # Preprocessed pair tables
exps/                  # Experiment outputs (graphs, logs)
build_graph.sh         # Main execution script
```

## Setup

```bash
mamba create -n fitness-landscape-graph python=3.11
mamba activate fitness-landscape-graph
pip install -e .
```

## Running

```bash
# Recommended: use the shell script
bash build_graph.sh

# Or run directly
python -m fitness_landscape_graph.build_graph \
    --base-path . \
    --output-dir ./outputs \
    --neutral-threshold 0.4 \
    --tiny-initial-threshold 0.02 \
    --large-edge-threshold 5.5 \
    --num-forbidden-pairs 1
```

Output: GraphML files (one per antibiotic per concentration) + global fitness graphs.

## Code Style
- Ruff for linting and formatting (configured in pyproject.toml)
- Google-style docstrings
- Type hints for all public functions/classes
- Always document array/tensor shapes at the point of creation/return (inline comments preferred). Shape comments use format: `# [B, R, M, 3]`. Add dtype or other info after shape comments when it improves clarity.
- Use assert only for internal invariants, not user validation

## Key Concepts

- **Mutant profile**: 13-char string (`.` = wildtype, letter = mutation, `X` = dead). E.g. `..VK.T..K...D`
- **Genotype**: 13-position mutation profile (e.g., `P.LS.T.SKM[LQ]D`)
- **Fitness**: Normalized AUC-derived resistance measurement (0 = dead, 1 = wildtype)
- **Pairs**: All mutant pairs differing by exactly one position, with fitness_diff
- **Peak node**: Local fitness maximum in the mutation network (out_degree=0, in_degree>0)
- **Forbidden pair**: Mutant pair with large fitness gap, protected from merging
- **Neutral merge**: Union-Find clustering of nodes with small fitness differences

## Pipeline

1. Preprocess raw AUC data → filter intended mutations, handle missing replicates
2. Compute normalized global fitness via trapezoid AUC integration
3. Build raw directed multigraph (edges = single-mutation transitions, weights = exp(|fitness_diff|))
4. Merge neutral clusters (Union-Find), respecting forbidden pairs
5. Detect fitness peaks, merge non-critical ancestors into peak clusters
6. Generate sequence logos, export as GraphML
