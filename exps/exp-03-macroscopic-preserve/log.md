check src/fitness_landscape_graph/graph_builder.py for the algorithm we are discussing here.

One of the limitation of current fitness landscape graph is that this neutral merging is too aggressive: two supernodes will be merged together is there is one edge falls below the neutral threshold. This limitation results in some macroscopic structure visible in lower neutral threshold disappear in high neutral threshold. However, we can't simple use neutral threshold because the whole point of this landscape graph representation is too coarse grain the landscpe such that the visualization is readable and interpretable for human.

My current solution to this landscape graph is forbidden pair. However, this solution introduce two new hyper parameter, and is not that clever.

Is it possible to design a better solution such that the landscape graph is more robust to the choice of neutral threshold, and the macroscopic structure is preserved in a automatic, intelligent way.

## Diagnosis: why the current approach fails

The current neutral merging uses Union-Find, which is equivalent to single-linkage clustering. The fundamental problem is that **Union-Find enforces transitivity**: if A merges with B and B merges with C, then A and C are in the same cluster — even if their fitness values are very different.

Neutrality is not transitive. A ≈ B and B ≈ C does not imply A ≈ C. A chain of individually small fitness differences can traverse a large fitness gap. This is why increasing the neutral threshold inevitably destroys macroscopic structure: longer chains form, bridging regions that should remain separate.

The forbidden pair mechanism is a manual patch for this — it identifies specific large edges and prevents merging across them. But it introduces two extra hyperparameters (`large_edge_threshold`, `num_forbidden_pairs`) and doesn't address the root cause.

Note that at the raw graph level, every node is a single genotype with at most 13 neighbors (one per mutation position), and each pair has at most one edge. Approaches based on aggregate inter-cluster edge statistics (average-linkage, edge density, etc.) are meaningless at this level — every pair has 0 or 1 edges. The problem is purely about transitive chaining along the graph.

## Idea 1: Fitness diameter constraint

**Core idea**: Within any merged cluster, the fitness diameter (max fitness - min fitness) should not exceed some bound.

During the Union-Find process, before merging two clusters, check whether the resulting cluster's fitness range would exceed a maximum allowed spread. If so, reject the merge even though the edge between them is small.

This is analogous to complete-linkage clustering in the fitness dimension, applied on the graph structure (only graph neighbors can initiate merges).

**Parameter**: one interpretable parameter — maximum allowed fitness spread within a cluster.

**Strength**: directly prevents the transitive chaining problem. A long chain of small steps that spans a large fitness range will be blocked once the cluster diameter exceeds the bound.

**Concern**: the result may still depend on merge order / starting point. Two different traversal orders through the sorted edge list could produce different clusterings, because which merges happen first determines which later merges are blocked by the diameter constraint. Need to investigate whether this order dependence is significant in practice.

## Idea 2: Basin-of-attraction decomposition (topological approach)

**Core idea**: instead of bottom-up clustering based on pairwise similarity, let the fitness function's own topology define the coarse-graining. The macroscopic structure of a fitness landscape is defined by:

1. **Peaks** — local fitness maxima (adaptive optima)
2. **Basins of attraction** — which peak does each genotype flow toward under selection?
3. **Barriers (saddle points)** — the fitness value where two adjacent basins first connect

These are properties of the fitness function's topology, not of pairwise neighbor similarity.

**Algorithm**:

- **Step 1: Assign every genotype to a basin.** From each node, follow the steepest uphill edge repeatedly until reaching a peak. This is deterministic and requires zero parameters. Every node ends up assigned to exactly one peak's basin.
- **Step 2: Compute barrier heights between adjacent basins.** For each pair of adjacent basins, the barrier height is determined by the saddle point — the fitness value at which the two basins first become connected. Concretely, this is related to the maximum fitness among boundary nodes (nodes in one basin that have a neighbor in the other basin). The barrier significance = how much the saddle fitness differs from the peak fitnesses.
- **Step 3: Merge basins separated by insignificant barriers.** If the barrier between two adjacent basins is below a threshold (meaning the saddle is close in fitness to both peaks), merge them — they are neutral variants of each other.

**Parameter**: one parameter — minimum significant barrier height. This replaces `neutral_threshold`, `large_edge_threshold`, and `num_forbidden_pairs`.

**Key reframing**: the old question was "is this edge neutral?" (local, susceptible to transitive chaining). The new question is "is this barrier significant?" (global, no chaining possible). A barrier height is a single number per pair of adjacent basins. There is no transitivity issue: if basin A and basin B have a low barrier, and basin B and basin C have a low barrier, it tells you nothing about the A–C barrier. The landscape itself determines that.

**Connection to existing pipeline**: this approach would unify the current neutral merge and peak merge steps. The current pipeline does neutral merge → peak detection → peak merge. The basin decomposition does peak detection first (via steepest ascent), then merges basins — achieving both goals in a principled order.

**Open question**: how to handle flat regions where steepest ascent is ambiguous (multiple neighbors with equal or near-equal fitness). Options include breaking ties consistently (e.g., lexicographic), or treating truly flat connected components as a single node before basin assignment.