"""Stage 3 sample generation: 2gamma pairs in controlled separation bins.

Two different ``gamma_iso`` events are stitched.  Both primary photons have
the fixed 53.085 MeV initial energy; **no deposited-energy selection** is
applied.  The only technical requirement is that each source event has a
non-zero deposited energy so that its energy-weighted centroid is defined.

For each event the two centroids are computed from the crystal-centre map
(``data/utilities/ecal_module_positions.npz``), the separation is measured in
crystal pitches, and events are generated until every bin has ``N_PER_BIN``
entries.  The bins are 0-0.5, 0.5-1, 1-2, 2-3, 3-5, >5 pitch.

One background draw is added per stitched event (same model as
``pytrial/bkg_est_add``); the signal-only and +background files share the
signal truth and the separation.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
PYTRIAL = HERE.parents[2]
REPO = HERE.parents[3]
BKG_DIR = PYTRIAL / "bkg_est_add"
sys.path.insert(0, str(PYTRIAL))
sys.path.insert(0, str(BKG_DIR))

import bkg_model  # noqa: E402
import geometry  # noqa: E402
import stitching  # noqa: E402
from add_background import draw_injections  # noqa: E402

DATA = HERE / "data"
META = HERE / "meta.json"
POSITIONS = REPO / "data/utilities/ecal_module_positions.npz"
GAMMA = REPO / "data/processed/edep_of_each_evt_gamma_iso.root"

POOL_EVENTS = 100_000
N_PER_BIN = 2_500
SEED = 20261106
PRIMARY_ENERGY_MEV = 53.085
N_MOD = 568
EDGES_PITCH = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 5.0, np.inf])
N_BINS = len(EDGES_PITCH) - 1


def build_truth(primary_arrays):
    """Sparse per-primary truth for a list of (n_events, n_mod) arrays."""
    n_events = primary_arrays[0].shape[0]
    truth_crystal, truth_primary, truth_energy = [], [], []
    for k in range(n_events):
        tc, tp, te = [], [], []
        for p, arr in enumerate(primary_arrays):
            row = arr[k]
            nz = np.flatnonzero(row > 0)
            tc.extend(nz.tolist())
            tp.extend([p] * len(nz))
            te.extend(row[nz].astype(float).tolist())
        truth_crystal.append(tc)
        truth_primary.append(tp)
        truth_energy.append(te)
    return truth_crystal, truth_primary, truth_energy


def sample_knn_pairs(tree, centroids, rng, lo, hi, target, k=500):
    """Sample unique pairs whose centroid distance is in [lo, hi) via KNN."""
    m = len(centroids)
    pairs = []
    seen = set()
    k = min(int(k), m)
    while len(pairs) < target:
        seeds = rng.integers(0, m, size=256)
        for s in seeds:
            d, j = tree.query(centroids[s], k=k)
            mask = (d >= lo) & (d < hi)
            if not mask.any():
                continue
            cand = j[mask]
            if len(cand) == 0:
                continue
            p = int(cand[rng.integers(0, len(cand))])
            if p == s:
                continue
            key = (min(s, p), max(s, p))
            if key in seen:
                continue
            seen.add(key)
            pairs.append(key)
            if len(pairs) >= target:
                break
    return pairs


def sample_random_pairs(centroids, rng, lo, hi, target, batch=20_000):
    """Sample unique pairs whose centroid distance is in [lo, hi)."""
    m = len(centroids)
    pairs = []
    seen = set()
    while len(pairs) < target:
        i = rng.integers(0, m, size=batch)
        j = rng.integers(0, m, size=batch)
        d = np.linalg.norm(centroids[i] - centroids[j], axis=1)
        valid = (i != j) & (d >= lo) & (d < hi)
        for a, b in zip(i[valid], j[valid]):
            key = (min(int(a), int(b)), max(int(a), int(b)))
            if key in seen:
                continue
            seen.add(key)
            pairs.append(key)
            if len(pairs) >= target:
                break
    return pairs


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    bkg_meta = json.loads((BKG_DIR / "bkg_estimation.json").read_text())
    n_bkg = float(bkg_meta["rates"]["n_bkg_per_event"])
    spectrum = bkg_model.DepositionSpectrum()

    centre, counts = geometry.load_positions(POSITIONS)
    pitch = geometry.estimate_pitch(centre, counts)
    print(f"crystal map: {int((counts > 0).sum())}/{N_MOD} modules, pitch={pitch:.2f} mm")

    print(f"loading first {POOL_EVENTS} gamma_iso events")
    ga = stitching.load_edeps(GAMMA, n_events=POOL_EVENTS)
    e_ga = stitching.event_energies(ga)
    positive = e_ga > 0
    pool = np.flatnonzero(positive)
    print(f"positive-deposit gamma events: {len(pool)}/{len(e_ga)} "
          f"(zero-deposit events have no centroid)")

    # Energy-weighted centroids of all positive-deposit events.
    centroids = (ga[pool] @ centre) / e_ga[pool, None]
    tree = cKDTree(centroids)
    m = len(centroids)

    edges_mm = EDGES_PITCH * pitch
    rng = np.random.default_rng([SEED, 3])
    pairs_by_bin = []
    for b in range(N_BINS):
        lo, hi = edges_mm[b], edges_mm[b + 1]
        if b <= 2:
            pairs = sample_knn_pairs(tree, centroids, rng, lo, hi, N_PER_BIN,
                                     k=500 if b == 2 else 200)
        else:
            pairs = sample_random_pairs(centroids, rng, lo, hi, N_PER_BIN)
        pairs_by_bin.append(pairs)
        print(f"bin {b} ({EDGES_PITCH[b]:g}-{EDGES_PITCH[b+1]:g} pitch, "
              f"{lo:.1f}-{hi:.1f} mm): {len(pairs)} pairs")

    all_pairs = []
    sep_mm, sep_bin = [], []
    for b, pairs in enumerate(pairs_by_bin):
        for a, c in pairs:
            d = float(np.linalg.norm(centroids[a] - centroids[c]))
            all_pairs.append((a, c))
            sep_mm.append(d)
            sep_bin.append(b)
    sep_mm = np.asarray(sep_mm, dtype=np.float32)
    sep_bin = np.asarray(sep_bin, dtype=np.int32)
    sep_pitch = sep_mm / pitch
    total = len(all_pairs)

    src_a = np.array([pool[a] for a, _ in all_pairs], dtype=np.int32)
    src_b = np.array([pool[c] for _, c in all_pairs], dtype=np.int32)
    e_a = e_ga[src_a].astype(np.float32)
    e_b = e_ga[src_b].astype(np.float32)
    signal_energy = (e_a + e_b).astype(np.float32)
    species = np.ones((total, 2), dtype=np.int32)
    src = np.stack([src_a, src_b], axis=1).astype(np.int32)
    src_e = np.stack([e_a, e_b], axis=1).astype(np.float32)
    primary_e = np.full((total, 2), PRIMARY_ENERGY_MEV, dtype=np.float32)

    print("building per-primary truth and summed Edeps")
    edeps_sig = (ga[src_a] + ga[src_b]).astype(np.float32)
    truth = build_truth([ga[src_a], ga[src_b]])

    bg_draw = draw_injections(total, N_MOD, n_bkg, spectrum, SEED, 3)
    bg_vecs, bg_mods, bg_energies = stitching.background_vectors(*bg_draw)
    bkg_total = bg_vecs.sum(axis=1).astype(np.float32)
    edeps_bkg = (edeps_sig + bg_vecs).astype(np.float32)

    ntrack = np.full(total, 2, dtype=np.int32)
    zero_nbkg = np.zeros(total, dtype=np.int32)
    zero_bkg = np.zeros(total, dtype=np.float32)
    empty = [[] for _ in range(total)]

    sig_path = DATA / "stage3_2gamma_sig.root"
    bkg_path = DATA / "stage3_2gamma_bkg.root"
    stitching.write_multiparticle_sample(
        sig_path, edeps_sig, ntrack, species, src, primary_e, src_e,
        signal_energy, zero_bkg, zero_nbkg, *truth, empty, empty,
        separation=sep_mm, sep_bin=sep_bin)
    stitching.write_multiparticle_sample(
        bkg_path, edeps_bkg, ntrack, species, src, primary_e, src_e,
        signal_energy, bkg_total, bg_draw[0], *truth, bg_mods, bg_energies,
        separation=sep_mm, sep_bin=sep_bin)

    max_truth = float(np.max(np.abs(
        np.array([sum(e) for e in truth[2]], dtype=np.float64) - signal_energy)))
    max_bkg = float(np.max(np.abs(
        edeps_bkg.sum(axis=1) - (signal_energy + bkg_total))))
    bin_counts = np.bincount(sep_bin, minlength=N_BINS)
    bin_medians = [float(np.median(sep_mm[sep_bin == b])) for b in range(N_BINS)]

    meta = {
        "stage": "3",
        "date": "2026-11-09",
        "seed": SEED,
        "pool_events": POOL_EVENTS,
        "positive_pool": int(m),
        "n_events": int(total),
        "n_per_bin": int(N_PER_BIN),
        "pitch_mm": float(pitch),
        "edges_pitch": EDGES_PITCH.tolist(),
        "edges_mm": edges_mm.tolist(),
        "bin_counts": bin_counts.tolist(),
        "bin_median_mm": bin_medians,
        "separation_min_mm": float(sep_mm.min()),
        "separation_max_mm": float(sep_mm.max()),
        "signal_energy_mean_mev": float(signal_energy.mean()),
        "signal_energy_min_mev": float(signal_energy.min()),
        "signal_energy_max_mev": float(signal_energy.max()),
        "source_energy_mean_by_primary_mev": [float(e_a.mean()), float(e_b.mean())],
        "background": {
            "n_bkg_per_event": n_bkg,
            "spectrum_mean_mev": spectrum.mean_deposition_mev,
            "frac_events_with_bkg": float(np.mean(bg_draw[0] > 0)),
        },
        "max_abs_truth_minus_signal_mev": max_truth,
        "max_abs_edeps_bkg_minus_signal_plus_bkg_mev": max_bkg,
        "note": ("Two different gamma_iso events; no deposited-energy selection. "
                 "Only zero-deposit events are excluded because their centroid is undefined."),
        "sig_file": str(sig_path.relative_to(REPO)),
        "bkg_file": str(bkg_path.relative_to(REPO)),
    }
    META.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    print(f"wrote {total} events; bin counts={bin_counts.tolist()}; "
          f"E_sig mean={signal_energy.mean():.2f} MeV")
    print(f"wrote {META}")


if __name__ == "__main__":
    main()
