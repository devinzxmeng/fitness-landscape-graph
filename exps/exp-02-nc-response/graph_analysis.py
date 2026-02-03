"""Helper functions for analyzing fitness landscape graph neighbors.

Provides utilities to extract 1-mutation and 2-mutation neighbors of a target
genotype from the pairwise table, and compute fitness advantages across
drug concentrations.
"""

import json
from statistics import median

import networkx as nx
import numpy as np
import polars as pl


def hamming_distance(a: str, b: str) -> int:
    """Compute Hamming distance between two equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def get_1mut_neighbors_fitness(
    pairs_df: pl.DataFrame,
    target: str,
    concentration: float,
) -> pl.DataFrame:
    """Get fitness differences for all 1-mutation neighbors at a concentration.

    Uses the pair table directly since it contains all 1-mutation pairs.

    Args:
        pairs_df: Pairwise table with columns mutant_profile1, mutant_profile2,
            concentration, median_diff, and replicate columns.
        target: 13-char mutant profile of the target genotype.
        concentration: Drug concentration to filter on.

    Returns:
        DataFrame with columns: neighbor, fitness_diff (target - neighbor),
            target_fitness, neighbor_fitness.
    """
    conc_df = pairs_df.filter(pl.col("concentration") == concentration)

    # Target as profile1: median_diff = profile1_median - profile2_median
    as_p1 = conc_df.filter(pl.col("mutant_profile1") == target).select(
        neighbor=pl.col("mutant_profile2"),
        target_rep1=pl.col("profile1_rep1"),
        target_rep2=pl.col("profile1_rep2"),
        target_rep3=pl.col("profile1_rep3"),
        neighbor_rep1=pl.col("profile2_rep1"),
        neighbor_rep2=pl.col("profile2_rep2"),
        neighbor_rep3=pl.col("profile2_rep3"),
    )

    # Target as profile2: median_diff = profile1_median - profile2_median
    as_p2 = conc_df.filter(pl.col("mutant_profile2") == target).select(
        neighbor=pl.col("mutant_profile1"),
        target_rep1=pl.col("profile2_rep1"),
        target_rep2=pl.col("profile2_rep2"),
        target_rep3=pl.col("profile2_rep3"),
        neighbor_rep1=pl.col("profile1_rep1"),
        neighbor_rep2=pl.col("profile1_rep2"),
        neighbor_rep3=pl.col("profile1_rep3"),
    )

    combined = pl.concat([as_p1, as_p2])

    # Compute median fitness for target and each neighbor
    rows = []
    for row in combined.iter_rows(named=True):
        t_fit = median(
            [row["target_rep1"], row["target_rep2"], row["target_rep3"]]
        )
        n_fit = median(
            [row["neighbor_rep1"], row["neighbor_rep2"], row["neighbor_rep3"]]
        )
        rows.append(
            {
                "neighbor": row["neighbor"],
                "target_fitness": t_fit,
                "neighbor_fitness": n_fit,
                "fitness_diff": t_fit - n_fit,
            }
        )

    return pl.DataFrame(rows)


def get_2mut_neighbors_fitness(
    pairs_df: pl.DataFrame,
    target: str,
    concentration: float,
) -> pl.DataFrame:
    """Get fitness differences for all 2-mutation neighbors at a concentration.

    Identifies 2-mutation neighbors by Hamming distance from the target,
    then extracts their fitness from the pair table's replicate columns.

    Args:
        pairs_df: Pairwise table.
        target: 13-char mutant profile of the target genotype.
        concentration: Drug concentration to filter on.

    Returns:
        DataFrame with columns: neighbor, fitness_diff (target - neighbor),
            target_fitness, neighbor_fitness.
    """
    conc_df = pairs_df.filter(pl.col("concentration") == concentration)

    # Get all unique genotypes at this concentration
    all_genotypes = set(conc_df["mutant_profile1"].unique().to_list()) | set(
        conc_df["mutant_profile2"].unique().to_list()
    )

    # Find 2-mutation neighbors
    neighbors_2mut = {g for g in all_genotypes if hamming_distance(g, target) == 2}

    if not neighbors_2mut:
        return pl.DataFrame(
            schema={
                "neighbor": pl.Utf8,
                "target_fitness": pl.Float64,
                "neighbor_fitness": pl.Float64,
                "fitness_diff": pl.Float64,
            }
        )

    # Extract target fitness from any row where target appears
    target_row = conc_df.filter(pl.col("mutant_profile1") == target).head(1)
    if target_row.height == 0:
        target_row = conc_df.filter(pl.col("mutant_profile2") == target).head(1)
        if target_row.height == 0:
            return pl.DataFrame(
                schema={
                    "neighbor": pl.Utf8,
                    "target_fitness": pl.Float64,
                    "neighbor_fitness": pl.Float64,
                    "fitness_diff": pl.Float64,
                }
            )
        t_fit = median(
            [
                target_row["profile2_rep1"][0],
                target_row["profile2_rep2"][0],
                target_row["profile2_rep3"][0],
            ]
        )
    else:
        t_fit = median(
            [
                target_row["profile1_rep1"][0],
                target_row["profile1_rep2"][0],
                target_row["profile1_rep3"][0],
            ]
        )

    # For each 2-mut neighbor, extract fitness from pair table
    rows = []
    for nb in neighbors_2mut:
        nb_row = conc_df.filter(pl.col("mutant_profile1") == nb).head(1)
        if nb_row.height > 0:
            n_fit = median(
                [
                    nb_row["profile1_rep1"][0],
                    nb_row["profile1_rep2"][0],
                    nb_row["profile1_rep3"][0],
                ]
            )
        else:
            nb_row = conc_df.filter(pl.col("mutant_profile2") == nb).head(1)
            if nb_row.height == 0:
                continue
            n_fit = median(
                [
                    nb_row["profile2_rep1"][0],
                    nb_row["profile2_rep2"][0],
                    nb_row["profile2_rep3"][0],
                ]
            )

        rows.append(
            {
                "neighbor": nb,
                "target_fitness": t_fit,
                "neighbor_fitness": n_fit,
                "fitness_diff": t_fit - n_fit,
            }
        )

    return pl.DataFrame(rows)


def compute_fitness_advantage_across_concentrations(
    pairs_df: pl.DataFrame,
    target: str,
    concentrations: list[float] | None = None,
    max_distance: int = 2,
) -> pl.DataFrame:
    """Compute fitness advantage of target over neighbors at all concentrations.

    Args:
        pairs_df: Pairwise table.
        target: 13-char mutant profile.
        concentrations: List of concentrations to analyze. If None, uses all.
        max_distance: Maximum mutation distance (1 or 2).

    Returns:
        DataFrame with columns: concentration, neighbor, distance,
            target_fitness, neighbor_fitness, fitness_diff.
    """
    if concentrations is None:
        concentrations = sorted(pairs_df["concentration"].unique().to_list())

    all_rows = []
    for conc in concentrations:
        # 1-mutation neighbors
        df_1mut = get_1mut_neighbors_fitness(pairs_df, target, conc)
        if df_1mut.height > 0:
            df_1mut = df_1mut.with_columns(
                pl.lit(conc).alias("concentration"),
                pl.lit(1).alias("distance"),
            )
            all_rows.append(df_1mut)

        # 2-mutation neighbors
        if max_distance >= 2:
            df_2mut = get_2mut_neighbors_fitness(pairs_df, target, conc)
            if df_2mut.height > 0:
                df_2mut = df_2mut.with_columns(
                    pl.lit(conc).alias("concentration"),
                    pl.lit(2).alias("distance"),
                )
                all_rows.append(df_2mut)

    if not all_rows:
        return pl.DataFrame(
            schema={
                "neighbor": pl.Utf8,
                "target_fitness": pl.Float64,
                "neighbor_fitness": pl.Float64,
                "fitness_diff": pl.Float64,
                "concentration": pl.Float64,
                "distance": pl.Int32,
            }
        )

    return pl.concat(all_rows)


def load_graph(graphml_path: str) -> nx.DiGraph:
    """Load a GraphML file and deserialize JSON-encoded node attributes.

    Args:
        graphml_path: Path to a .graphml file.

    Returns:
        nx.DiGraph with group_mutants as dict, numeric attributes cast properly.
    """
    graph = nx.read_graphml(graphml_path)
    for node in graph.nodes():
        data = graph.nodes[node]
        if isinstance(data.get("group_mutants"), str):
            data["group_mutants"] = json.loads(data["group_mutants"])
        for key in ("fitness", "group_size", "is_peak", "contain_wildtype"):
            if key in data and isinstance(data[key], str):
                try:
                    data[key] = float(data[key])
                except ValueError:
                    pass
        if "group_size" in data:
            data["group_size"] = int(data["group_size"])
        if "is_peak" in data:
            data["is_peak"] = int(data["is_peak"])
    return graph


def find_genotype_node(graph: nx.DiGraph, target: str) -> str | None:
    """Find which node contains a specific genotype in its group_mutants.

    Args:
        graph: A fitness landscape graph.
        target: A 13-character mutant profile string.

    Returns:
        The node label containing the target, or None.
    """
    for node, data in graph.nodes(data=True):
        if target in data.get("group_mutants", {}):
            return node
    return None
