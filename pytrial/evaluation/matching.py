"""Truth-level matching for multi-particle samples (Stage 2/3).

Planned API
-----------
``build_overlap_matrix(clusters, truth)``
    ``M[i, j]`` = energy of reconstructed cluster ``i`` that comes from truth
    primary ``j``.
``match_clusters(overlap)``
    One-to-one assignment (Hungarian) maximising the total matched energy.
``efficiency_purity(overlap, assignment)``
    Per-primary efficiency and per-cluster purity.
``merge_split_rates(overlap, assignment)``
    Merge rate, split rate and matched fraction.

The per-primary truth format (``TruthCrystal`` / ``TruthPrimary`` /
``TruthEnergy``) is specified in ``MULTIPARTICLE_STAGE_PLAN.md`` §5.3.
"""

from __future__ import annotations


def build_overlap_matrix(clusters, truth):
    raise NotImplementedError("Implemented in Stage 2 (Scheme B).")


def match_clusters(overlap):
    raise NotImplementedError("Implemented in Stage 2 (Scheme B).")


def efficiency_purity(overlap, assignment):
    raise NotImplementedError("Implemented in Stage 2 (Scheme B).")


def merge_split_rates(overlap, assignment):
    raise NotImplementedError("Implemented in Stage 2 (Scheme B).")
