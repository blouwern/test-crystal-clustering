"""Stage 2 sample generation: 2gamma and positron+2gamma, with gamma background.

Signal scenarios (the primary particles are monochromatic at 53.085 MeV; **no
deposited-energy selection is applied**):

1. ``2gamma``   : two different gamma_iso events, stitched.
2. ``ep2gamma`` : one ep_iso (positron) + two different gamma_iso events.

The primary energies are fixed by the samples (53.085 MeV each); the deposited
signal energy ``E_sig = sum(Edeps)`` varies only through detector leakage.  It
is used as the energy denominator in the metrics, not as an event-selection
criterion.

Background (same model as ``pytrial/bkg_est_add``): ``k ~ Poisson(n_bkg)``,
``n_bkg = 2.42466`` hits per event, uniform random crystals, hit energies from
N(5, 2^2) truncated at zero.  **One background draw per stitched event**, not
one per signal particle.

For each scenario two files are written with identical signal truth:

* ``*_sig.root``: signal only (control);
* ``*_bkg.root``: signal + background.

Sparse per-primary truth is stored in ``TruthCrystal`` / ``TruthPrimary`` /
``TruthEnergy``; background truth in ``BkgCrystal`` / ``BkgEnergy``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PYTRIAL = HERE.parents[2]
REPO = HERE.parents[3]
BKG_DIR = PYTRIAL / "bkg_est_add"
sys.path.insert(0, str(PYTRIAL))
sys.path.insert(0, str(BKG_DIR))

import bkg_model  # noqa: E402
import stitching  # noqa: E402
from add_background import draw_injections  # noqa: E402

DATA = HERE / "data"
META = HERE / "meta.json"

POOL_EVENTS = 100_000
N_EVENTS = 20_000
SEED = 20261106
PRIMARY_ENERGY_MEV = 53.085
N_MOD = 568


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


def write_pair(name, edeps_sig, edeps_bkg, species, source_idx, source_energy,
               signal_energy, truth, bkg):
    """Write the signal-only and background versions of one scenario."""
    n_events = edeps_sig.shape[0]
    ntrack = np.full(n_events, species.shape[1], dtype=np.int32)
    primary_energy = np.full(species.shape, PRIMARY_ENERGY_MEV, dtype=np.float32)
    bkg_vecs = bkg["vecs"]
    bkg_mods = bkg["mods"]
    bkg_energies = bkg["energies"]
    bkg_total = bkg_vecs.sum(axis=1).astype(np.float32)
    zero_nbkg = np.zeros(n_events, dtype=np.int32)
    zero_bkg = np.zeros(n_events, dtype=np.float32)
    empty = [[] for _ in range(n_events)]
    tc, tp, te = truth

    sig_path = DATA / f"stitched_{name}_sig.root"
    bkg_path = DATA / f"stitched_{name}_bkg.root"
    stitching.write_multiparticle_sample(
        sig_path, edeps_sig, ntrack, species, source_idx, primary_energy,
        source_energy, signal_energy, zero_bkg, zero_nbkg,
        tc, tp, te, empty, empty)
    stitching.write_multiparticle_sample(
        bkg_path, edeps_bkg, ntrack, species, source_idx, primary_energy,
        source_energy, signal_energy, bkg_total, bkg["counts"],
        tc, tp, te, bkg_mods, bkg_energies)

    # Validation is done on the in-memory arrays before writing; record it.
    max_truth = float(np.max(np.abs(
        np.array([float(np.sum(e)) for e in te], dtype=np.float64)
        - np.asarray(signal_energy, dtype=np.float64))))
    max_sig = float(np.max(np.abs(
        edeps_sig.sum(axis=1) - np.asarray(signal_energy, dtype=np.float64))))
    max_bkg = float(np.max(np.abs(
        edeps_bkg.sum(axis=1) - (np.asarray(signal_energy) + bkg_total))))
    return {
        "sig_file": str(sig_path.relative_to(REPO)),
        "bkg_file": str(bkg_path.relative_to(REPO)),
        "n_events": int(n_events),
        "ntrack": int(species.shape[1]),
        "max_abs_truth_minus_signal_mev": max_truth,
        "max_abs_edeps_sig_minus_signal_mev": max_sig,
        "max_abs_edeps_bkg_minus_signal_plus_bkg_mev": max_bkg,
        "signal_energy_mean_mev": float(np.mean(signal_energy)),
        "signal_energy_min_mev": float(np.min(signal_energy)),
        "signal_energy_max_mev": float(np.max(signal_energy)),
        "bkg_energy_mean_mev": float(np.mean(bkg_total)),
        "nbkg_mean": float(np.mean(bkg["counts"])),
        "frac_events_with_bkg": float(np.mean(bkg["counts"] > 0)),
    }


def main():
    DATA.mkdir(parents=True, exist_ok=True)

    bkg_meta = json.loads((BKG_DIR / "bkg_estimation.json").read_text())
    n_bkg = float(bkg_meta["rates"]["n_bkg_per_event"])
    spectrum = bkg_model.DepositionSpectrum()

    print(f"loading first {POOL_EVENTS} events of ep_iso and gamma_iso")
    ep = stitching.load_edeps(
        REPO / "data/processed/edep_of_each_evt_ep_iso.root", n_events=POOL_EVENTS)
    ga = stitching.load_edeps(
        REPO / "data/processed/edep_of_each_evt_gamma_iso.root", n_events=POOL_EVENTS)
    e_ep = stitching.event_energies(ep)
    e_ga = stitching.event_energies(ga)
    print(f"ep: mean E_sig={e_ep.mean():.2f} MeV; gamma: mean E_sig={e_ga.mean():.2f} MeV")
    print(f"background: n_bkg={n_bkg:.6f} hits/event, spectrum mean="
          f"{spectrum.mean_deposition_mev:.4f} MeV")

    # ------------------------------------------------------------------ #
    # Scenario 1: two different gamma events, no deposited-energy selection
    # ------------------------------------------------------------------ #
    rng = np.random.default_rng([SEED, 1])
    pool = np.arange(len(e_ga), dtype=np.int64)
    print(f"2gamma pool: all {len(pool)} gamma events (no E_sig cut)")
    i1 = rng.choice(pool, size=N_EVENTS, replace=True)
    i2 = rng.choice(pool, size=N_EVENTS, replace=True)
    same = i1 == i2
    while same.any():
        i2[same] = rng.choice(pool, size=int(same.sum()), replace=True)
        same = i1 == i2

    sig_2g = ga[i1] + ga[i2]
    species_2g = np.ones((N_EVENTS, 2), dtype=np.int32)
    src_2g = np.stack([i1, i2], axis=1).astype(np.int32)
    src_e_2g = np.stack([e_ga[i1], e_ga[i2]], axis=1).astype(np.float32)
    sig_e_2g = src_e_2g.sum(axis=1)
    truth_2g = build_truth([ga[i1], ga[i2]])
    bg_2g_draw = draw_injections(N_EVENTS, N_MOD, n_bkg, spectrum, SEED, 1)
    bg_2g_vecs, bg_2g_mods, bg_2g_energies = stitching.background_vectors(*bg_2g_draw)
    bkg_2g = {"counts": bg_2g_draw[0], "vecs": bg_2g_vecs,
              "mods": bg_2g_mods, "energies": bg_2g_energies}
    edeps_2g_sig = sig_2g.astype(np.float32)
    edeps_2g_bkg = (sig_2g + bg_2g_vecs).astype(np.float32)
    meta_2g = write_pair("2gamma", edeps_2g_sig, edeps_2g_bkg, species_2g,
                         src_2g, src_e_2g, sig_e_2g, truth_2g, bkg_2g)
    meta_2g.update({
        "selection": "two different gamma_iso events; no deposited-energy selection",
        "pool_size": int(len(pool)),
        "source_energy_mean_by_primary_mev": [
            float(e_ga[i1].mean()), float(e_ga[i2].mean())],
    })
    print(f"2gamma: wrote {meta_2g['n_events']} events; "
          f"E_sig mean={meta_2g['signal_energy_mean_mev']:.2f} MeV")

    # ------------------------------------------------------------------ #
    # Scenario 2: one positron + two different gammas, no energy selection
    # ------------------------------------------------------------------ #
    rng = np.random.default_rng([SEED, 2])
    ep_pool = np.arange(len(e_ep), dtype=np.int64)
    ga_pool = np.arange(len(e_ga), dtype=np.int64)
    ie = rng.choice(ep_pool, size=N_EVENTS, replace=True)
    j1 = rng.choice(ga_pool, size=N_EVENTS, replace=True)
    j2 = rng.choice(ga_pool, size=N_EVENTS, replace=True)
    same = j1 == j2
    while same.any():
        j2[same] = rng.choice(ga_pool, size=int(same.sum()), replace=True)
        same = j1 == j2
    sig_e_ep2g = (e_ep[ie] + e_ga[j1] + e_ga[j2]).astype(np.float32)

    sig_ep2g = ep[ie] + ga[j1] + ga[j2]
    species_ep2g = np.tile(np.array([0, 1, 1], dtype=np.int32), (N_EVENTS, 1))
    src_ep2g = np.stack([ie, j1, j2], axis=1).astype(np.int32)
    src_e_ep2g = np.stack([e_ep[ie], e_ga[j1], e_ga[j2]], axis=1).astype(np.float32)
    truth_ep2g = build_truth([ep[ie], ga[j1], ga[j2]])
    bg_ep2g_draw = draw_injections(N_EVENTS, N_MOD, n_bkg, spectrum, SEED, 2)
    bg_ep2g_vecs, bg_ep2g_mods, bg_ep2g_energies = stitching.background_vectors(*bg_ep2g_draw)
    bkg_ep2g = {"counts": bg_ep2g_draw[0], "vecs": bg_ep2g_vecs,
                "mods": bg_ep2g_mods, "energies": bg_ep2g_energies}
    edeps_ep2g_sig = sig_ep2g.astype(np.float32)
    edeps_ep2g_bkg = (sig_ep2g + bg_ep2g_vecs).astype(np.float32)
    meta_ep2g = write_pair("ep2gamma", edeps_ep2g_sig, edeps_ep2g_bkg, species_ep2g,
                           src_ep2g, src_e_ep2g, sig_e_ep2g, truth_ep2g, bkg_ep2g)
    meta_ep2g.update({
        "selection": ("one ep_iso positron + two different gamma_iso events; "
                      "no deposited-energy selection"),
        "pool_size_ep": int(len(ep_pool)),
        "pool_size_gamma": int(len(ga_pool)),
        "source_energy_mean_by_primary_mev": [
            float(e_ep[ie].mean()), float(e_ga[j1].mean()), float(e_ga[j2].mean())],
    })
    print(f"ep2gamma: wrote {meta_ep2g['n_events']} events; "
          f"E_sig mean={meta_ep2g['signal_energy_mean_mev']:.2f} MeV")

    meta = {
        "stage": "2",
        "date": "2026-11-09",
        "seed": SEED,
        "pool_events": POOL_EVENTS,
        "n_events": N_EVENTS,
        "primary_energy_mev": PRIMARY_ENERGY_MEV,
        "primary_energy_note": (
            "SimPrimaryVertex.Ek0 is 53.085 MeV for every ep_iso/gamma_iso event. "
            "No deposited-energy selection is applied; the deposited signal energy "
            "E_sig = sum(Edeps) is the detector response and is used only in the metrics."
        ),
        "background": {
            "model": "pytrial/bkg_est_add",
            "n_bkg_per_event": n_bkg,
            "spectrum_mean_mev": spectrum.mean_deposition_mev,
            "spectrum_sigma_mev": spectrum.sigma_mev,
            "draw": "one Poisson draw per stitched event",
        },
        "scenarios": {"2gamma": meta_2g, "ep2gamma": meta_ep2g},
    }
    META.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    print(f"wrote {META}")


if __name__ == "__main__":
    main()
