# AGENTS.md

Guidelines for agentic coding agents working on the fitness-landscape-graph codebase.

## Build/Lint/Test Commands

```bash
# Setup environment
mamba create -n fitness-landscape-graph python=3.11
mamba activate fitness-landscape-graph
pip install -e .

# Run main pipeline
bash build_graph.sh

# Or run directly
python -m fitness_landscape_graph.build_graph \
    --base-path . \
    --output-dir ./outputs \
    --neutral-threshold 0.4 \
    --tiny-initial-threshold 0.02 \
    --large-edge-threshold 5.5 \
    --num-forbidden-pairs 1

# Lint and format (required before commits)
ruff check src/          # Check for issues
ruff check --fix src/     # Auto-fix issues
ruff format src/          # Format code

# Single test (if tests exist - none currently)
pytest path/to/test_file.py::test_function -v

# All tests
pytest
```

## Tech Stack

- Python 3.11+
- Polars (DataFrames), NetworkX (graphs), NumPy, SciPy (stats), Logomaker, Matplotlib
- Ruff for linting/formatting (configured in pyproject.toml)
- Gephi (external) for graph visualization

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
```

## Code Style Guidelines

### Imports
- Order: stdlib → third-party → local (isort handles this)
- Use `from fitness_landscape_graph import X` for local imports
- Combine related imports with parentheses
- Example:
  ```python
  import argparse
  import json
  
  import networkx as nx
  import polars as pl
  
  from fitness_landscape_graph.graph_builder import GraphBuilder
  ```

### Formatting
- Line length: 88 characters (Black-compatible)
- 4 spaces for indentation
- Double quotes for strings
- Trailing commas for multi-line structures
- Run `ruff format` before committing

### Type Hints
- Use Python 3.10+ union syntax: `int | None`, `list[str]`
- Type all public functions and class methods
- Use `pl.DataFrame`, `nx.Graph` for domain types
- Example: `def build_graph(self, concentration: float | None = None) -> nx.MultiDiGraph`

### Docstrings
- Google-style docstrings when writing documentation
- Include Args and Returns sections
- Docstring formatting is optional (D10* rules ignored)

### Naming Conventions
- snake_case for functions and variables
- PascalCase for classes
- UPPER_CASE for module-level constants
- Descriptive names (avoid single letters except loop indices)

### Error Handling
- Use `assert` only for internal invariants, not user input validation
- Use `contextlib.suppress()` for expected exceptions when appropriate
- Log errors with the `logging` module (see `logger = logging.getLogger(__name__)` pattern)

### Shape Comments
Always document array shapes at creation/return points using format: `# [B, R, M, 3]`
Add dtype or other info after shape comments when it improves clarity.
Example:
```python
# [n_mutants, n_concs, n_replicates]
reshaped_data = data_array.reshape(n_mutants, n_concs, n_replicates)
```

## Key Concepts

- **Mutant profile**: 13-char string (`.` = wildtype, letter = mutation, `X` = dead). E.g. `..VK.T..K...D`
- **Genotype**: 13-position mutation profile notation like `P.LS.T.SKM[LQ]D`
- **Fitness**: Normalized AUC-derived resistance measurement (0 = dead, 1 = wildtype)
- **Pairs**: All mutant pairs differing by exactly one position, with fitness_diff
- **Peak node**: Local fitness maximum (out_degree=0, in_degree>0)
- **Forbidden pair**: Mutant pair with large fitness gap, protected from merging
- **Neutral merge**: Union-Find clustering of nodes with small fitness differences
- Graph files saved as GraphML format

## Pipeline Overview

1. Preprocess raw AUC data → filter intended mutations, handle missing replicates
2. Compute normalized global fitness via trapezoid AUC integration
3. Build raw directed multigraph (edges = single-mutation transitions, weights = exp(|fitness_diff|))
4. Merge neutral clusters (Union-Find), respecting forbidden pairs
5. Detect fitness peaks, merge non-critical ancestors into peak clusters
6. Generate sequence logos, export as GraphML

### Pre-commit Checklist
- [ ] `ruff check src/` passes with no errors
- [ ] `ruff format src/` applied
- [ ] Type hints present for public APIs
- [ ] Shape comments added for array operations
- [ ] No secrets or hardcoded paths committed
