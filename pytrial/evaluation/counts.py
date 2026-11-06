"""Count-based metrics for multi-particle clustering.

These metrics only need the per-event cluster count and the true multiplicity
N.  They are intentionally independent of any truth-level matching, so they can
be used for the Stage 1 smoke test and for every later stage.

The historical single-particle scorer (``pytrial/scoring.py``) reads the true
number of events from ``max(EvtID)+1``; that silently drops trailing events
with zero clusters.  Here the true event count is always supplied explicitly.
"""

from __future__ import annotations

import numpy as np
import ROOT


def read_cluster_counts(cluster_path, n_events):
    """Return an array of length ``n_events`` with the cluster count per event.

    Events without a cluster entry stay at zero, so the denominator is always
    the true number of input events.
    """
    n_events = int(n_events)
    counts = np.zeros(n_events, dtype=np.int64)
    with ROOT.TFile.Open(str(cluster_path), "READ") as f:
        tree = f.Get("ECALClusters")
        if tree.GetEntries() == 0:
            return counts
        for entry in tree:
            evt = int(entry.EvtID)
            if 0 <= evt < n_events:
                counts[evt] += 1
    return counts


def count_metrics(counts, n_truth):
    """Count-based metrics for a fixed true multiplicity ``n_truth``.

    Returns a JSON-serialisable dict with the exact/under/over rates, MAE,
    RMSE, bias, distribution moments and a histogram.  ``frac_0`` ... are
    included as conveniences (``frac_3plus`` means ``>= 3``).
    """
    c = np.asarray(counts, dtype=np.int64)
    if c.size == 0:
        raise ValueError("counts is empty")
    n_truth = int(n_truth)
    diff = c.astype(np.int64) - n_truth
    abs_diff = np.abs(diff)

    values, frequencies = np.unique(c, return_counts=True)
    histogram = {int(v): int(f) for v, f in zip(values, frequencies)}

    result = {
        "n_events": int(c.size),
        "n_truth": n_truth,
        "exact_rate": float(np.mean(diff == 0)),
        "under_rate": float(np.mean(diff < 0)),
        "over_rate": float(np.mean(diff > 0)),
        "MAE": float(np.mean(abs_diff)),
        "RMSE": float(np.sqrt(np.mean(diff.astype(np.float64) ** 2))),
        "bias": float(np.mean(diff)),
        "mean_clusters": float(np.mean(c)),
        "std_clusters": float(np.std(c, ddof=1)) if c.size > 1 else 0.0,
        "min_clusters": int(c.min()),
        "max_clusters": int(c.max()),
        "histogram": histogram,
    }
    for k in (0, 1, 2):
        result[f"frac_{k}"] = float(np.mean(c == k))
    result["frac_3plus"] = float(np.mean(c >= 3))
    result["frac_lt"] = result["under_rate"]
    result["frac_eq"] = result["exact_rate"]
    result["frac_gt"] = result["over_rate"]
    return result
