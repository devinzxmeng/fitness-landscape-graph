# exp-02-nc-response

## Purpose

Addresses Reviewer #4's comments from the Nature Communications response. See `nc-response.md` for the full reviewer comments and interpretation.

## Deliverables

### 1. Fitness advantage over neighbors (Reviewer #4, point iii)

Shows how the peak's fitness advantage over its 1-mutation and 2-mutation neighbors changes with drug concentration. Distinguishes real biological flattening from threshold artifacts.

**Run:** Open and run `fitness_advantage_analysis.ipynb`

**Output:** `outputs/fitness-advantage-azt-12/`
- `fitness_advantage_boxplot.pdf` — box plot of fitness advantage vs concentration
- `fitness_advantage_data.csv` — raw data

### 2. Robustness/sensitivity of global peak existence (Reviewer #4, point iii)

Shows that the global peak disappearance at intermediate concentrations is robust across a range of neutrality cutoff values (0.15–0.45), not an artifact of a specific threshold choice.

**Run:**
```bash
bash run_global_peak_analysis.sh
```

**Output:** `outputs/global-peak-robustness/`
- `global_peak_heatmap.png` — threshold × concentration heatmap
- `global_peak_by_concentration.png` — summary by concentration
- `global_peak_by_threshold.png` — summary by threshold
- `global_peak_analysis.csv` — raw analysis data

**Key scripts:**
- `run_global_peak_analysis.sh` — master workflow (build → analyze → visualize)
- `analyze_global_peaks.py` — batch analysis of all graphs for global peak existence
- `plot_global_peak_heatmap.py` — visualization of results

## Code organization

Reusable analysis code has been moved to `src/fitness_landscape_graph/`:
- `graph_analyzer.py` — GraphML loading, visualization, global peak detection
- `fitness_advantage.py` — vectorized fitness advantage computation
- `build_graphs_parallel.py` — parallel graph building for threshold sweeps

## Dev notes (historical)

- Neutral threshold sweeping analysis required group_size > 32 for global peak detection
- Additional mode to detect big connection nodes was considered but not needed
