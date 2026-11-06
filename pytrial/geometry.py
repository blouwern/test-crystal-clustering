"""Crystal geometry helpers for multi-particle stitching.

The detector topology file only stores the neighbour graph; for separation
studies (Stage 3) we also need the crystal centres.  The raw MC hit tree
``G4Run0/ECALSimHit`` stores the hit position ``x`` (3-vector) and ``ModID``, so
the crystal centre is the mean hit position per module - the same approach as
``data/utilities/verify_neighbor_graph.py``.

Only the standard library and numpy are needed for the geometry maths; ROOT is
imported lazily inside :func:`read_module_positions` so the pure functions can
be used without opening a file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

N_MOD = 568
DEFAULT_RAW = "data/raw/ep_iso.root"
DEFAULT_POSITIONS = "data/utilities/ecal_module_positions.npz"


def read_module_positions(raw_path, n_hits=500_000, n_mod=N_MOD):
    """Return ``(centre, counts)`` for every module.

    ``centre[i]`` is the mean MC hit position of module ``i`` (mm) and
    ``counts[i]`` the number of hits used.  Modules that were never hit have a
    zero centre and zero count; they can be masked with ``counts > 0``.
    """
    import ROOT

    raw_path = str(raw_path)
    tree = "G4Run0/ECALSimHit"
    df = ROOT.RDataFrame(tree, raw_path).Range(0, int(n_hits))
    data = df.AsNumpy(["ModID", "x"])
    mod = np.asarray(data["ModID"], dtype=np.int64)
    pts = np.stack([np.asarray(v, dtype=np.float64) for v in data["x"]])

    acc = np.zeros((int(n_mod), 3), dtype=np.float64)
    cnt = np.zeros(int(n_mod), dtype=np.int64)
    np.add.at(acc, mod, pts)
    np.add.at(cnt, mod, 1)
    centre = np.divide(acc, cnt[:, None], out=np.zeros_like(acc),
                       where=cnt[:, None] > 0)
    return centre, cnt


def save_positions(path, centre, counts):
    """Save a position map to a compressed ``.npz`` file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, centre=np.asarray(centre, dtype=np.float64),
                        counts=np.asarray(counts, dtype=np.int64))
    return path


def load_positions(path):
    """Load a position map saved by :func:`save_positions`."""
    data = np.load(path)
    return data["centre"], data["counts"]


def energy_weighted_centroid(edeps, centre):
    """Energy-weighted centroid of one event (mm)."""
    e = np.asarray(edeps, dtype=np.float64)
    total = float(e.sum())
    if total <= 0.0:
        return np.full(3, np.nan)
    return (e[:, None] * np.asarray(centre, dtype=np.float64)).sum(axis=0) / total


def pairwise_distances(centroids):
    """Symmetric pairwise Euclidean distances (mm) between centroids."""
    c = np.asarray(centroids, dtype=np.float64)
    n = len(c)
    out = np.full((n, n), np.nan, dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(c[i] - c[j]))
            out[i, j] = out[j, i] = d
    return out


def min_pairwise_distance(centroids):
    """Smallest centroid pair distance; NaN for fewer than two centroids."""
    d = pairwise_distances(centroids)
    if d.shape[0] < 2:
        return float("nan")
    return float(np.nanmin(d))


def estimate_pitch(centre, counts):
    """Median nearest-neighbour distance among hit modules (mm)."""
    mask = counts > 0
    p = np.asarray(centre, dtype=np.float64)[mask]
    if len(p) < 2:
        return float("nan")
    d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    return float(np.median(d.min(axis=1)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", default=DEFAULT_RAW,
                        help="raw MC file with G4Run0/ECALSimHit")
    parser.add_argument("--out", default=DEFAULT_POSITIONS,
                        help="output .npz position map")
    parser.add_argument("--n-hits", type=int, default=500_000,
                        help="number of MC hits to average (default 500k)")
    parser.add_argument("--n-mod", type=int, default=N_MOD)
    args = parser.parse_args()

    centre, counts = read_module_positions(args.raw, args.n_hits, args.n_mod)
    path = save_positions(args.out, centre, counts)
    hit = int((counts > 0).sum())
    pitch = estimate_pitch(centre, counts)
    print(f"read {args.n_hits} hits from {args.raw}")
    print(f"modules with hits: {hit}/{args.n_mod}")
    print(f"estimated pitch: {pitch:.1f} mm")
    print(f"saved {path}")


if __name__ == "__main__":
    main()
