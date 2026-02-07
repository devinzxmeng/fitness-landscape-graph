"""Tests for FitnessAdvantageAnalyzer.

End-to-end tests with hand-crafted data to verify correct neighbor
detection and fitness difference calculations. Also includes performance
test for vectorized Hamming distance computation.

Test setup (used across tests):

    Genotypes (3-char for simplicity):
        AAA  fitness=0.8  (group member)
        AAB  fitness=0.6  (1-mut from AAA)
        ABA  fitness=0.4  (1-mut from AAA)
        ABB  fitness=0.3  (2-mut from AAA, 1-mut from AAB and ABA)
        BAA  fitness=0.5  (1-mut from AAA)
        BBA  fitness=0.7  (2-mut from AAA)

    1-mutation pairs in pair table:
        AAA-AAB, AAA-ABA, AAA-BAA, AAB-ABB, ABA-ABB, BAA-BBA

    Group = {AAA}
    Expected 1-mut external neighbors: AAB, ABA, BAA
    Expected 2-mut external neighbors: ABB, BBA
"""

import time

import polars as pl
import pytest

from fitness_landscape_graph.fitness_advantage import FitnessAdvantageAnalyzer


@pytest.fixture()
def sample_data():
    """Create hand-crafted pairs_df and long_df for testing.

    Returns (pairs_df, long_df) at concentration=1.0.
    """
    # Pair table: all 1-mutation pairs with pre-computed median_diff
    # median_diff = median(profile1) - median(profile2)
    pairs_df = pl.DataFrame(
        {
            "mutant_profile1": ["AAA", "AAA", "AAA", "AAB", "ABA", "BAA"],
            "mutant_profile2": ["AAB", "ABA", "BAA", "ABB", "ABB", "BBA"],
            "concentration": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            "median_diff": [
                0.8 - 0.6,   # AAA(0.8) - AAB(0.6) = 0.2
                0.8 - 0.4,   # AAA(0.8) - ABA(0.4) = 0.4
                0.8 - 0.5,   # AAA(0.8) - BAA(0.5) = 0.3
                0.6 - 0.3,   # AAB(0.6) - ABB(0.3) = 0.3
                0.4 - 0.3,   # ABA(0.4) - ABB(0.3) = 0.1
                0.5 - 0.7,   # BAA(0.5) - BBA(0.7) = -0.2
            ],
        }
    )

    # Long table: per-genotype fitness
    long_df = pl.DataFrame(
        {
            "mutant_profile": ["AAA", "AAB", "ABA", "ABB", "BAA", "BBA"],
            "concentration": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            "median": [0.8, 0.6, 0.4, 0.3, 0.5, 0.7],
        }
    )

    return pairs_df, long_df


@pytest.fixture()
def multi_conc_data():
    """Create data with two concentrations for testing.

    Concentration 1.0: same as sample_data
    Concentration 2.0: fitness values shifted (AAA still highest in group)
    """
    pairs_df = pl.DataFrame(
        {
            "mutant_profile1": [
                "AAA", "AAA", "AAA", "AAB", "ABA", "BAA",
                "AAA", "AAA", "AAA", "AAB", "ABA", "BAA",
            ],
            "mutant_profile2": [
                "AAB", "ABA", "BAA", "ABB", "ABB", "BBA",
                "AAB", "ABA", "BAA", "ABB", "ABB", "BBA",
            ],
            "concentration": [
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
            ],
            "median_diff": [
                0.2, 0.4, 0.3, 0.3, 0.1, -0.2,
                0.1, 0.2, 0.15, 0.15, 0.05, -0.1,
            ],
        }
    )

    long_df = pl.DataFrame(
        {
            "mutant_profile": [
                "AAA", "AAB", "ABA", "ABB", "BAA", "BBA",
                "AAA", "AAB", "ABA", "ABB", "BAA", "BBA",
            ],
            "concentration": [
                1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
            ],
            "median": [
                0.8, 0.6, 0.4, 0.3, 0.5, 0.7,
                0.5, 0.4, 0.3, 0.25, 0.35, 0.45,
            ],
        }
    )

    return pairs_df, long_df


class TestEndToEnd:
    """End-to-end tests verifying correct neighbors and fitness differences."""

    def test_finds_correct_1mut_neighbors(self, sample_data):
        """Group={AAA} should have 1-mut neighbors: AAB, ABA, BAA."""
        pairs_df, long_df = sample_data
        analyzer = FitnessAdvantageAnalyzer(pairs_df, long_df)

        result = analyzer.compute_fitness_advantage(
            group_genotypes={"AAA"},
            concentrations=[1.0],
            max_distance=1,
        )

        neighbors = set(result["neighbor"].to_list())
        assert neighbors == {"AAB", "ABA", "BAA"}
        assert all(d == 1 for d in result["distance"].to_list())

    def test_finds_correct_2mut_neighbors(self, sample_data):
        """Group={AAA} should have 2-mut neighbors: ABB, BBA."""
        pairs_df, long_df = sample_data
        analyzer = FitnessAdvantageAnalyzer(pairs_df, long_df)

        result = analyzer.compute_fitness_advantage(
            group_genotypes={"AAA"},
            concentrations=[1.0],
            max_distance=2,
        )

        dist2 = result.filter(pl.col("distance") == 2)
        neighbors_2mut = set(dist2["neighbor"].to_list())
        assert neighbors_2mut == {"ABB", "BBA"}

    def test_correct_fitness_diff_1mut(self, sample_data):
        """Verify fitness_diff = group_member_fitness - neighbor_fitness for 1-mut."""
        pairs_df, long_df = sample_data
        analyzer = FitnessAdvantageAnalyzer(pairs_df, long_df)

        result = analyzer.compute_fitness_advantage(
            group_genotypes={"AAA"},
            concentrations=[1.0],
            max_distance=1,
        )

        for row in result.iter_rows(named=True):
            if row["neighbor"] == "AAB":
                assert abs(row["fitness_diff"] - 0.2) < 1e-10  # 0.8 - 0.6
                assert abs(row["group_member_fitness"] - 0.8) < 1e-10
                assert abs(row["neighbor_fitness"] - 0.6) < 1e-10
            elif row["neighbor"] == "ABA":
                assert abs(row["fitness_diff"] - 0.4) < 1e-10  # 0.8 - 0.4
            elif row["neighbor"] == "BAA":
                assert abs(row["fitness_diff"] - 0.3) < 1e-10  # 0.8 - 0.5

    def test_correct_fitness_diff_2mut(self, sample_data):
        """Verify fitness_diff = group_member_fitness - neighbor_fitness for 2-mut."""
        pairs_df, long_df = sample_data
        analyzer = FitnessAdvantageAnalyzer(pairs_df, long_df)

        result = analyzer.compute_fitness_advantage(
            group_genotypes={"AAA"},
            concentrations=[1.0],
            max_distance=2,
        )

        dist2 = result.filter(pl.col("distance") == 2)
        for row in dist2.iter_rows(named=True):
            if row["neighbor"] == "ABB":
                assert abs(row["fitness_diff"] - 0.5) < 1e-10  # 0.8 - 0.3
            elif row["neighbor"] == "BBA":
                assert abs(row["fitness_diff"] - 0.1) < 1e-10  # 0.8 - 0.7

    def test_excludes_group_members_from_neighbors(self, sample_data):
        """When group has multiple members, they shouldn't appear as neighbors."""
        pairs_df, long_df = sample_data
        analyzer = FitnessAdvantageAnalyzer(pairs_df, long_df)

        # Group = {AAA, AAB}. AAB is 1-mut from AAA but should NOT be a neighbor.
        result = analyzer.compute_fitness_advantage(
            group_genotypes={"AAA", "AAB"},
            concentrations=[1.0],
            max_distance=2,
        )

        neighbors = set(result["neighbor"].to_list())
        assert "AAA" not in neighbors
        assert "AAB" not in neighbors

    def test_multiple_concentrations(self, multi_conc_data):
        """Verify results are produced for each concentration."""
        pairs_df, long_df = multi_conc_data
        analyzer = FitnessAdvantageAnalyzer(pairs_df, long_df)

        result = analyzer.compute_fitness_advantage(
            group_genotypes={"AAA"},
            concentrations=[1.0, 2.0],
            max_distance=1,
        )

        concs = set(result["concentration"].to_list())
        assert concs == {1.0, 2.0}

        # At conc 2.0, AAA-AAB fitness_diff = 0.1 (from median_diff in pairs_df)
        conc2 = result.filter(
            (pl.col("concentration") == 2.0) & (pl.col("neighbor") == "AAB")
        )
        assert abs(conc2["fitness_diff"][0] - 0.1) < 1e-10

    def test_output_columns(self, sample_data):
        """Verify output DataFrame has all expected columns."""
        pairs_df, long_df = sample_data
        analyzer = FitnessAdvantageAnalyzer(pairs_df, long_df)

        result = analyzer.compute_fitness_advantage(
            group_genotypes={"AAA"},
            concentrations=[1.0],
            max_distance=2,
        )

        expected_cols = {
            "group_member", "neighbor", "distance",
            "group_member_fitness", "neighbor_fitness",
            "fitness_diff", "concentration",
        }
        assert set(result.columns) == expected_cols

    def test_empty_result(self, sample_data):
        """Group with no external neighbors returns empty DataFrame."""
        pairs_df, long_df = sample_data
        analyzer = FitnessAdvantageAnalyzer(pairs_df, long_df)

        # All genotypes in the group → no external neighbors
        all_genotypes = set(long_df["mutant_profile"].to_list())
        result = analyzer.compute_fitness_advantage(
            group_genotypes=all_genotypes,
            concentrations=[1.0],
            max_distance=2,
        )

        assert result.height == 0


class TestHammingPerformance:
    """Performance test for vectorized Hamming distance computation."""

    def test_performance_benchmark(self):
        """Vectorized Hamming should handle 100×1000 genotypes under 1 second."""
        n_group = 100
        n_external = 1000

        group = [f"A{i:012d}" for i in range(n_group)]
        external = [f"B{i:012d}" for i in range(n_external)]

        start = time.time()
        result = FitnessAdvantageAnalyzer._compute_hamming_distances_vectorized(
            group, external, target_distance=2
        )
        elapsed = time.time() - start

        assert elapsed < 1.0
        assert isinstance(result, list)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
