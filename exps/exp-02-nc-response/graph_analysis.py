"""Helper functions for analyzing fitness landscape graph neighbors.

Provides utilities to extract 1-mutation and 2-mutation neighbors of a group
of genotypes (a peak supernode) from the pairwise table, and compute fitness
advantages across drug concentrations.
"""

import json
from statistics import median

import networkx as nx
import polars as pl


def hamming_distance(a: str, b: str) -> int:
    """Compute Hamming distance between two equal-length strings."""
    return sum(c1 != c2 for c1, c2 in zip(a, b))


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


def get_peak_genotypes(graph: nx.DiGraph, rank: int = 0) -> set[str]:
    """Get all genotype strings from a peak supernode.

    Args:
        graph: A fitness landscape graph.
        rank: 0 for the highest-fitness peak, 1 for second, etc.

    Returns:
        Set of genotype strings in the peak's group_mutants.
    """
    peaks = [
        (n, d) for n, d in graph.nodes(data=True) if d.get("is_peak") == 1
    ]
    peaks.sort(key=lambda x: x[1]["fitness"], reverse=True)
    if rank >= len(peaks):
        return set()
    _, data = peaks[rank]
    return set(data.get("group_mutants", {}).keys())


def _get_genotype_fitness(
    conc_df: pl.DataFrame, genotype: str
) -> float | None:
    """Extract median fitness of a genotype from the pair table at one concentration.

    Args:
        conc_df: Pair table already filtered to one concentration.
        genotype: 13-char mutant profile.

    Returns:
        Median fitness, or None if not found.
    """
    row = conc_df.filter(pl.col("mutant_profile1") == genotype).head(1)
    if row.height > 0:
        vals = [row["profile1_rep1"][0], row["profile1_rep2"][0], row["profile1_rep3"][0]]
        vals = [v for v in vals if v is not None and v == v]  # filter None and NaN
        return median(vals) if vals else None
    row = conc_df.filter(pl.col("mutant_profile2") == genotype).head(1)
    if row.height > 0:
        vals = [row["profile2_rep1"][0], row["profile2_rep2"][0], row["profile2_rep3"][0]]
        vals = [v for v in vals if v is not None and v == v]
        return median(vals) if vals else None
    return None


def get_group_external_neighbors(
    pairs_df: pl.DataFrame,
    group_genotypes: set[str],
    concentration: float,
    max_distance: int = 2,
) -> pl.DataFrame:
    """Get fitness differences between group members and external neighbors.

    For each genotype in the group, finds neighbors (by mutation distance)
    that are NOT in the group, and computes the fitness difference.

    Args:
        pairs_df: Pairwise table.
        group_genotypes: Set of genotype strings in the peak supernode.
        concentration: Drug concentration to filter on.
        max_distance: Maximum mutation distance (1 or 2).

    Returns:
        DataFrame with columns: group_member, neighbor, distance,
            group_member_fitness, neighbor_fitness, fitness_diff.
    """
    conc_df = pairs_df.filter(pl.col("concentration") == concentration)

    # For 1-mutation neighbors, use the pair table directly
    # Filter pairs where one member is in the group and the other is not
    in_group_p1 = conc_df.filter(
        pl.col("mutant_profile1").is_in(group_genotypes)
        & ~pl.col("mutant_profile2").is_in(group_genotypes)
    ).select(
        group_member=pl.col("mutant_profile1"),
        neighbor=pl.col("mutant_profile2"),
        gm_rep1=pl.col("profile1_rep1"),
        gm_rep2=pl.col("profile1_rep2"),
        gm_rep3=pl.col("profile1_rep3"),
        nb_rep1=pl.col("profile2_rep1"),
        nb_rep2=pl.col("profile2_rep2"),
        nb_rep3=pl.col("profile2_rep3"),
    )

    in_group_p2 = conc_df.filter(
        pl.col("mutant_profile2").is_in(group_genotypes)
        & ~pl.col("mutant_profile1").is_in(group_genotypes)
    ).select(
        group_member=pl.col("mutant_profile2"),
        neighbor=pl.col("mutant_profile1"),
        gm_rep1=pl.col("profile2_rep1"),
        gm_rep2=pl.col("profile2_rep2"),
        gm_rep3=pl.col("profile2_rep3"),
        nb_rep1=pl.col("profile1_rep1"),
        nb_rep2=pl.col("profile1_rep2"),
        nb_rep3=pl.col("profile1_rep3"),
    )

    combined_1mut = pl.concat([in_group_p1, in_group_p2])

    rows = []
    for row in combined_1mut.iter_rows(named=True):
        gm_vals = [v for v in [row["gm_rep1"], row["gm_rep2"], row["gm_rep3"]] if v is not None and v == v]
        nb_vals = [v for v in [row["nb_rep1"], row["nb_rep2"], row["nb_rep3"]] if v is not None and v == v]
        if not gm_vals or not nb_vals:
            continue
        gm_fit = median(gm_vals)
        nb_fit = median(nb_vals)
        rows.append(
            {
                "group_member": row["group_member"],
                "neighbor": row["neighbor"],
                "distance": 1,
                "group_member_fitness": gm_fit,
                "neighbor_fitness": nb_fit,
                "fitness_diff": gm_fit - nb_fit,
            }
        )

    # For 2-mutation neighbors, find external genotypes at Hamming distance 2
    if max_distance >= 2:
        all_genotypes = set(conc_df["mutant_profile1"].unique().to_list()) | set(
            conc_df["mutant_profile2"].unique().to_list()
        )
        external_genotypes = all_genotypes - group_genotypes

        # For each group member, find 2-mut external neighbors
        for gm in group_genotypes:
            gm_fit = _get_genotype_fitness(conc_df, gm)
            if gm_fit is None:
                continue
            for ext in external_genotypes:
                if hamming_distance(gm, ext) != 2:
                    continue
                nb_fit = _get_genotype_fitness(conc_df, ext)
                if nb_fit is None:
                    continue
                rows.append(
                    {
                        "group_member": gm,
                        "neighbor": ext,
                        "distance": 2,
                        "group_member_fitness": gm_fit,
                        "neighbor_fitness": nb_fit,
                        "fitness_diff": gm_fit - nb_fit,
                    }
                )

    if not rows:
        return pl.DataFrame(
            schema={
                "group_member": pl.Utf8,
                "neighbor": pl.Utf8,
                "distance": pl.Int32,
                "group_member_fitness": pl.Float64,
                "neighbor_fitness": pl.Float64,
                "fitness_diff": pl.Float64,
            }
        )

    return pl.DataFrame(rows)


def compute_group_fitness_advantage(
    pairs_df: pl.DataFrame,
    group_genotypes: set[str],
    concentrations: list[float] | None = None,
    max_distance: int = 2,
) -> pl.DataFrame:
    """Compute fitness advantage of group members over external neighbors.

    Args:
        pairs_df: Pairwise table.
        group_genotypes: Set of genotype strings in the peak supernode.
        concentrations: Concentrations to analyze. If None, uses all.
        max_distance: Maximum mutation distance (1 or 2).

    Returns:
        DataFrame with columns: concentration, group_member, neighbor,
            distance, group_member_fitness, neighbor_fitness, fitness_diff.
    """
    if concentrations is None:
        concentrations = sorted(pairs_df["concentration"].unique().to_list())

    all_dfs = []
    for conc in concentrations:
        print(f"  Processing concentration {conc}...")
        df = get_group_external_neighbors(
            pairs_df, group_genotypes, conc, max_distance
        )
        if df.height > 0:
            df = df.with_columns(pl.lit(conc).alias("concentration"))
            all_dfs.append(df)

    if not all_dfs:
        return pl.DataFrame(
            schema={
                "group_member": pl.Utf8,
                "neighbor": pl.Utf8,
                "distance": pl.Int32,
                "group_member_fitness": pl.Float64,
                "neighbor_fitness": pl.Float64,
                "fitness_diff": pl.Float64,
                "concentration": pl.Float64,
            }
        )

    return pl.concat(all_dfs)
