#!/usr/bin/env python3
"""Validate the background-injected samples produced by ``add_background.py``.

For every output file the script

* reads back the parameters stored in the ``bkg_model`` metadata object;
* regenerates the random draws deterministically from the stored seed;
* re-reads the original and the injected file event by event and checks that
  (a) events with ``NBkg == 0`` are bit-identical,
  (b) events with ``NBkg > 0`` differ exactly by the regenerated deposits,
  (c) the per-event ``NBkg`` values match the Poisson draws;
* compares the observed ``NBkg`` distribution with Poisson(n_bkg), the mean
  energy per hit with the assumed spectrum, and the crystal-ID usage with a
  uniform distribution.

Usage
-----
    python verify_background.py
    python verify_background.py --datasets mu_rest --max-events 100000
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import ROOT
from scipy import stats

import add_background
import bkg_model

ROOT.gROOT.SetBatch(True)


def read_metadata(path: Path) -> dict:
    f = ROOT.TFile.Open(str(path), "READ")
    obj = f.Get("bkg_model")
    if obj is None:
        f.Close()
        raise RuntimeError(f"{path} has no bkg_model metadata (was it written by add_background.py?)")
    meta = json.loads(str(obj.GetTitle()))
    f.Close()
    return meta


def verify_dataset(input_path: Path, output_path: Path, max_events: int = 0) -> dict:
    meta = read_metadata(output_path)
    n_bkg = float(meta["n_bkg_per_event"])
    seed = int(meta["seed"])
    stream = int(meta["stream"])

    fin = ROOT.TFile.Open(str(input_path), "READ")
    fout = ROOT.TFile.Open(str(output_path), "READ")
    tin, tout = fin.Get(add_background.TREE_NAME), fout.Get(add_background.TREE_NAME)
    n_entries = int(tout.GetEntries())
    if int(tin.GetEntries()) < n_entries:
        raise RuntimeError(f"{input_path} has fewer entries than {output_path}")
    n_events = n_entries if max_events <= 0 else min(int(max_events), n_entries)

    vin = ROOT.std.vector("float")()
    vout = ROOT.std.vector("float")()
    tin.SetBranchAddress(add_background.BRANCH_NAME, vin)
    tout.SetBranchAddress(add_background.BRANCH_NAME, vout)
    nb_buffer = np.zeros(1, dtype=np.int32)
    tout.SetBranchAddress("NBkg", nb_buffer)

    tin.GetEntry(0)
    n_crystals = int(np.asarray(vin).size)

    counts, mods, energies = add_background.draw_injections(
        n_events, n_crystals, n_bkg, bkg_model.DepositionSpectrum(
            meta["estimate"]["spectrum"]["nominal_mean_mev"],
            meta["estimate"]["spectrum"]["variance_mev2"],
        ), seed, stream,
    )

    nbkg_seen = np.zeros(8, dtype=np.int64)
    zero_events_mismatch = 0
    deposit_mismatch = 0
    nbkg_mismatch = 0
    max_abs_error = 0.0
    changed_total = 0
    pos = 0

    for i in range(n_events):
        tin.GetEntry(i)
        tout.GetEntry(i)
        k = int(counts[i])
        nbkg_seen[min(k, 7)] += 1
        if int(nb_buffer[0]) != k:
            nbkg_mismatch += 1

        a = np.asarray(vin)
        b = np.asarray(vout)
        if k == 0:
            if not np.array_equal(a, b):
                zero_events_mismatch += 1
            continue

        idx = mods[pos:pos + k]
        dep = energies[pos:pos + k]
        pos += k
        expected = a.copy()
        np.add.at(expected, idx, dep)
        error = float(np.abs(b - expected).max())
        max_abs_error = max(max_abs_error, error)
        if error > 1e-3:
            deposit_mismatch += 1
        changed_total += int(np.count_nonzero(b != a))

    fin.Close()
    fout.Close()

    n_hits = int(counts.sum())                      # total Poisson draws
    n_events_affected = int(nbkg_seen[1:].sum())    # events with >= 1 hit
    n_hits_expected = n_events * n_bkg
    frac_expected = 1.0 - np.exp(-n_bkg)
    frac_observed = n_events_affected / n_events

    # Poisson goodness of fit.  Histogram bins 0..6 plus a ">= 7" tail; merge
    # the last bins until every expected count is at least 5.
    n_bins = len(nbkg_seen)
    obs = nbkg_seen.astype(float).copy()
    probs = [float(np.exp(-n_bkg) * n_bkg**k / math.factorial(k)) for k in range(n_bins - 1)]
    probs.append(max(0.0, 1.0 - sum(probs)))
    exp = n_events * np.array(probs)
    while len(obs) > 2 and (exp[-1] < 5.0 or exp[-2] < 5.0):
        obs[-2] += obs[-1]
        exp[-2] += exp[-1]
        obs, exp = obs[:-1], exp[:-1]
    chi2 = float(((obs - exp) ** 2 / np.maximum(exp, 1e-12)).sum())
    dof = int(len(exp) - 1)
    p_value = float(stats.chi2.sf(chi2, dof)) if dof > 0 else float("nan")

    # Uniformity of the crystal IDs.
    if n_hits > 0:
        usage = np.bincount(mods, minlength=n_crystals).astype(float)
        expected_usage = n_hits / n_crystals
        uniform_chi2 = float(((usage - expected_usage) ** 2 / expected_usage).sum())
        uniform_p = float(stats.chi2.sf(uniform_chi2, n_crystals - 1))
    else:
        uniform_chi2, uniform_p = float("nan"), float("nan")

    mean_e = float(energies.astype(np.float64).mean())
    std_e = float(energies.astype(np.float64).std(ddof=1)) if n_hits > 1 else float("nan")
    sem_e = std_e / np.sqrt(n_hits) if n_hits > 1 else float("nan")
    spectrum_mean = float(bkg_model.DepositionSpectrum(
        meta["estimate"]["spectrum"]["nominal_mean_mev"],
        meta["estimate"]["spectrum"]["variance_mev2"],
    ).mean_deposition_mev)

    checks = {
        "entry_count": bool(n_entries == int(meta["n_events"])),
        "n_crystals": bool(n_crystals == int(meta["n_crystals"])),
        "zero_events_identical": bool(zero_events_mismatch == 0),
        "deposits_match": bool(deposit_mismatch == 0),
        "nbkg_match": bool(nbkg_mismatch == 0),
        "poisson_ok": bool(p_value > 0.001),
        "uniform_crystals": bool(not (uniform_p <= 0.001)),
        "mean_energy_ok": bool(abs(mean_e - spectrum_mean) < 4 * sem_e + 1e-6),
    }

    return {
        "dataset": input_path.name,
        "n_events": n_events,
        "n_crystals": n_crystals,
        "n_bkg": n_bkg,
        "n_hits": n_hits,
        "n_events_affected": n_events_affected,
        "n_hits_expected": n_hits_expected,
        "frac_events_with_bkg": frac_observed,
        "frac_expected": frac_expected,
        "nbkg_counts": nbkg_seen.tolist(),
        "poisson_expected": exp.tolist(),
        "poisson_chi2": chi2,
        "poisson_p": p_value,
        "uniform_chi2": uniform_chi2,
        "uniform_p": uniform_p,
        "mean_energy_mev": mean_e,
        "std_energy_mev": std_e,
        "sem_energy_mev": sem_e,
        "spectrum_mean_mev": spectrum_mean,
        "changed_crystals": changed_total,
        "max_abs_error_mev": max_abs_error,
        "checks": checks,
        "pass": all(bool(v) for v in checks.values()),
    }


def print_report(r: dict) -> None:
    print(f"\n[{r['dataset']}]")
    print(f"  events / crystals            {r['n_events']:,} / {r['n_crystals']}")
    print(f"  background hits (draws)      {r['n_hits']:,}  (expected {r['n_hits_expected']:.1f})")
    print(f"  events with background       {r['n_events_affected']:,}"
          f"  ({r['frac_events_with_bkg']:.5f}, expected {r['frac_expected']:.5f})")
    counts = r["nbkg_counts"]
    hist = " / ".join(f"{counts[k]:,}" for k in range(min(5, len(counts))))
    tail = sum(counts[5:]) if len(counts) > 5 else 0
    if tail:
        hist += f" / {tail:,} (>=5)"
    print(f"  NBkg histogram               {hist}")
    print(f"  Poisson chi2 / p             {r['poisson_chi2']:.3f} / {r['poisson_p']:.3f}")
    print(f"  crystal uniformity chi2 / p  {r['uniform_chi2']:.1f} / {r['uniform_p']:.3f}")
    print(f"  mean energy per hit          {r['mean_energy_mev']:.4f} +- {r['sem_energy_mev']:.4f} MeV"
          f"  (model {r['spectrum_mean_mev']:.4f})")
    print(f"  max |delta - expected|       {r['max_abs_error_mev']:.3e} MeV")
    print(f"  changed crystal entries      {r['changed_crystals']:,}")
    failed = [k for k, v in r["checks"].items() if not v]
    print("  checks                       " + ("PASS" if not failed else f"FAIL: {failed}"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets", nargs="+", default=list(add_background.DEFAULT_DATASETS))
    parser.add_argument("--data-dir", type=Path, default=add_background.default_data_dir())
    parser.add_argument("--suffix", default="_bkg")
    parser.add_argument("--max-events", type=int, default=0,
                        help="verify only the first N events (0 = all)")
    parser.add_argument("--json", type=Path, default=None,
                        help="optionally write the verification summary as JSON")
    args = parser.parse_args()

    processed = args.data_dir / "processed"
    reports = []
    for name in args.datasets:
        reports.append(verify_dataset(
            processed / f"edep_of_each_evt_{name}.root",
            processed / f"edep_of_each_evt_{name}{args.suffix}.root",
            max_events=args.max_events,
        ))

    for report in reports:
        print_report(report)

    ok = all(r["pass"] for r in reports)
    print("\n" + ("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED"))
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(reports, indent=2), encoding="utf-8")
        print(f"wrote {args.json}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
