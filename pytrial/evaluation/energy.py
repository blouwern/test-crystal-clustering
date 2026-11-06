"""Per-primary energy metrics for multi-particle samples (Stage 2).

Planned API
-----------
``energy_metrics(clusters, truth, overlap)``
    Total energy capture, per-primary energy recovery, leading-cluster
    containment and unclustered energy for a stitched sample.

The single-particle counterpart lives in ``pytrial/sweeps/run_task.py``
(``energy_metrics``); the multi-particle version needs the per-primary truth
to attribute shared crystals correctly, so it is kept separate.
"""

from __future__ import annotations


def energy_metrics(clusters, truth, overlap=None):
    raise NotImplementedError("Implemented in Stage 2 (Scheme B).")
