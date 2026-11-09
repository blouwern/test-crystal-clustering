"""Truth-level matching for multi-particle samples.

The overlap matrix ``M[i, j]`` is the signal energy of primary ``j`` that lands
in reconstructed cluster ``i``.  Shared crystals are handled naturally: one
crystal can contribute to several primaries.

Matching uses the Hungarian algorithm to maximise the total matched signal
energy.  Merge/split diagnostics use each primary's *best* cluster (largest
overlap), because that is the quantity a merge or split manifests in.
"""

from __future__ import annotations

from itertools import combinations, permutations

import numpy as np
from scipy.optimize import linear_sum_assignment


def build_overlap_matrix(cluster_mods, signal_arrays):
    """Overlap matrix (n_clusters x n_primary) for one event."""
    n_clusters = len(cluster_mods)
    n_primary = len(signal_arrays)
    matrix = np.zeros((n_clusters, n_primary), dtype=np.float64)
    for i, mods in enumerate(cluster_mods):
        idx = np.asarray(mods, dtype=np.int64)
        for j in range(n_primary):
            matrix[i, j] = float(signal_arrays[j][idx].sum())
    return matrix


def assign(matrix):
    """Exact assignment, returned as ``cluster -> primary`` (-1 = none).

    Multi-particle events have at most a handful of clusters and primaries, so
    a brute-force search over permutations is faster than a SciPy call per
    event.  The result is the same maximum-total-overlap assignment as the
    Hungarian algorithm.  Larger matrices fall back to SciPy.
    """
    n_clusters, n_primary = matrix.shape
    assignment = np.full(n_clusters, -1, dtype=np.int64)
    if n_clusters == 0 or n_primary == 0:
        return assignment

    if n_clusters * n_primary <= 20:
        k = min(n_clusters, n_primary)
        best_total = -1.0
        best_pairs = []
        for cluster_subset in combinations(range(n_clusters), k):
            for primary_subset in permutations(range(n_primary), k):
                total = sum(matrix[c, p] for c, p in
                            zip(cluster_subset, primary_subset))
                if total > best_total:
                    best_total = total
                    best_pairs = list(zip(cluster_subset, primary_subset))
        for c, p in best_pairs:
            if matrix[c, p] > 0:
                assignment[c] = p
        return assignment

    size = max(n_clusters, n_primary)
    padded = np.zeros((size, size), dtype=np.float64)
    padded[:n_clusters, :n_primary] = matrix
    rows, cols = linear_sum_assignment(-padded)
    for r, c in zip(rows, cols):
        if r < n_clusters and c < n_primary and matrix[r, c] > 0:
            assignment[r] = c
    return assignment


def efficiency_purity(matrix, assignment, primary_totals):
    """Per-primary efficiency and per-cluster match purity."""
    n_clusters, n_primary = matrix.shape
    efficiency = np.zeros(n_primary, dtype=np.float64)
    for i, j in enumerate(assignment):
        if j >= 0 and primary_totals[j] > 0:
            efficiency[j] = matrix[i, j] / float(primary_totals[j])

    match_purity = np.zeros(n_clusters, dtype=np.float64)
    for i, j in enumerate(assignment):
        total = float(matrix[i].sum())
        if j >= 0 and total > 0:
            match_purity[i] = matrix[i, j] / total
    return efficiency, match_purity


def merge_split(matrix, primary_totals, split_threshold=0.9):
    """Merge and split diagnostics based on each primary's best cluster.

    Returns per-primary ``best_cluster`` / ``best_efficiency`` / ``split`` and
    the number of primaries involved in a merge (two primaries whose best
    cluster is the same).
    """
    n_clusters, n_primary = matrix.shape
    result = {
        "best_cluster": np.full(n_primary, -1, dtype=np.int64),
        "best_efficiency": np.zeros(n_primary, dtype=np.float64),
        "split": np.zeros(n_primary, dtype=bool),
        "merged_primaries": 0,
        "merge_event": False,
        "split_primaries": 0,
    }
    if n_clusters == 0 or n_primary == 0:
        return result

    best = np.argmax(matrix, axis=0)
    best_energy = matrix[best, np.arange(n_primary)]
    valid = best_energy > 0
    best_eff = np.zeros(n_primary, dtype=np.float64)
    best_eff[valid] = best_energy[valid] / np.maximum(primary_totals[valid], 1e-12)
    split = valid & (best_eff < split_threshold)

    groups = {}
    for j in range(n_primary):
        if valid[j]:
            groups.setdefault(int(best[j]), []).append(j)
    merged = set()
    for primaries in groups.values():
        if len(primaries) >= 2:
            merged.update(primaries)

    result.update(
        best_cluster=best,
        best_efficiency=best_eff,
        split=split,
        merged_primaries=len(merged),
        merge_event=len(merged) > 0,
        split_primaries=int(split.sum()),
    )
    return result
