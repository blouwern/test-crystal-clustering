"""Separation-bin helpers for Stage 3.

Distances are expressed in crystal pitches; the default binning follows the
Stage 3 plan: 0-0.5, 0.5-1, 1-2, 2-3, 3-5, >5 pitch.
"""

from __future__ import annotations

import numpy as np

DEFAULT_EDGES_PITCH = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 5.0, np.inf])


def assign_bins(sep_pitch, edges=DEFAULT_EDGES_PITCH):
    """Return the bin index (0..len(edges)-2) for each separation."""
    values = np.asarray(sep_pitch, dtype=np.float64)
    return np.digitize(values, np.asarray(edges, dtype=np.float64)[1:-1])


def bin_label(index, edges=DEFAULT_EDGES_PITCH):
    """Human-readable label for one bin index."""
    edges = np.asarray(edges, dtype=np.float64)
    lo, hi = edges[int(index)], edges[int(index) + 1]
    hi_s = "inf" if not np.isfinite(hi) else f"{hi:g}"
    return f"{lo:g}-{hi_s} pitch"


def crossing_radius(x, y, level=0.5):
    """First downward crossing of ``level``, linearly interpolated.

    ``x`` must be increasing (e.g. median separation per bin); ``y`` is the
    merge rate.  Returns NaN when the curve never drops below ``level``.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if len(x) < 2:
        return float("nan")
    for i in range(1, len(y)):
        if y[i - 1] >= level > y[i]:
            t = (level - y[i - 1]) / (y[i] - y[i - 1])
            return float(x[i - 1] + t * (x[i] - x[i - 1]))
    if y[0] < level:
        return float(x[0])
    return float("nan")
