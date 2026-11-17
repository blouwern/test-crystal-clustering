#!/usr/bin/env python3
"""Step 7: inject the conservative neutron background into the processed samples.

For every event the script

1. draws the number of captures ``k ~ Poisson(n_capture)``;
2. draws ``k`` crystal IDs uniformly in ``[0, N_crystals)`` (the flux is
   isotropic and the array is a spherical shell, so every crystal is equally
   likely to be hit from outside);
3. draws ``k`` deposited energies from the capture-depth model
   ``E = Q_gamma * g(d)`` with ``d ~ Exp(lambda_c)`` (mean 4.39 MeV, range
   3.4-6.1 MeV; see ``neutron_model.CaptureEnergy``);
4. adds those energies to the ``Edeps`` vector of the event.

All the deposited energy of one capture is put into a **single** crystal.  This
is conservative (it maximises the single-crystal energy), while the physical
cascade usually spreads over 1-2 neighbouring crystals; even so the single
crystal stays at or below ``Q_gamma = 6.86 MeV``, i.e. below the seed threshold.

The original files are never modified: the output is written next to them with
the ``_bkg_n`` suffix (``edep_of_each_evt_ep_iso_bkg_n.root`` etc.).  Each output
file keeps the ``EdepOfEachEvt`` tree with the ``Edeps`` branch, so it is a
drop-in replacement for the clustering pipeline, and adds an ``NNeutron`` branch
with the number of captures.  A ``bkg_model_neutron`` metadata object records the
parameters, the seed and the model result used for the file.

Usage
-----
    python add_neutron_background.py                       # all three samples, full 1e6 events
    python add_neutron_background.py --datasets mu_rest    # one sample only
    python add_neutron_background.py --max-events 10000    # quick test
    python add_neutron_background.py --store-bkg           # also store the background-only vector
"""

from __future__ import annotations

import argparse
import json
import os
import time
import zlib
from pathlib import Path

import numpy as np
import ROOT

import neutron_model

ROOT.gROOT.SetBatch(True)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEFAULT_DATASETS = ("ep_iso", "gamma_iso", "mu_rest")
TREE_NAME = "EdepOfEachEvt"
BRANCH_NAME = "Edeps"
COUNT_BRANCH = "NNeutron"


def default_data_dir() -> Path:
    env = os.environ.get("ECAL_CLUSTERING_DATA_DIR")
    return Path(env) if env else REPO / "data"


def draw_injections(
    n_events: int,
    n_crystals: int,
    n_capture: float,
    energy: neutron_model.CaptureEnergy,
    seed: int,
    stream: int,
):
    """Draw the neutron background for a whole file (deterministic given seed).

    Returns ``(counts, mods, energies)`` where ``counts[i]`` is the number of
    captures of event ``i`` and the flat arrays ``mods``/``energies`` hold the
    crystal IDs and deposits in event order.
    """
    rng = np.random.default_rng([int(seed), int(stream)])
    counts = rng.poisson(n_capture, size=int(n_events)).astype(np.int64)
    n_hits = int(counts.sum())
    mods = rng.integers(0, int(n_crystals), size=n_hits, dtype=np.int64).astype(np.int32)
    energies = energy.rvs(n_hits, rng)
    return counts, mods, energies


def inject_file(
    input_path: Path,
    output_path: Path,
    n_capture: float,
    energy: neutron_model.CaptureEnergy,
    seed: int,
    stream: int,
    max_events: int = 0,
    store_bkg: bool = False,
    metadata: dict | None = None,
) -> dict:
    """Add neutron background to one ROOT file and return a small summary."""
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    fin = ROOT.TFile.Open(str(input_path), "READ")
    if not fin or fin.IsZombie():
        raise RuntimeError(f"cannot open {input_path}")
    tin = fin.Get(TREE_NAME)
    if not tin:
        raise RuntimeError(f"{input_path} has no {TREE_NAME} tree")
    n_total = int(tin.GetEntries())

    vin = ROOT.std.vector("float")()
    tin.SetBranchAddress(BRANCH_NAME, vin)
    tin.GetEntry(0)
    n_crystals = int(np.asarray(vin).size)
    if n_crystals == 0:
        raise RuntimeError(f"{input_path}: empty {BRANCH_NAME} vector")

    n_events = n_total if max_events <= 0 else min(int(max_events), n_total)

    counts, mods, energies = draw_injections(
        n_events, n_crystals, n_capture, energy, seed, stream
    )
    n_hits = int(counts.sum())

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fout = ROOT.TFile.Open(str(output_path), "RECREATE")
    tout = ROOT.TTree(TREE_NAME, "Edep of each event (conservative neutron background added)")
    tout.Branch(BRANCH_NAME, vin)
    n_buffer = np.zeros(1, dtype=np.int32)
    tout.Branch(COUNT_BRANCH, n_buffer, f"{COUNT_BRANCH}/I")
    vbkg = None
    if store_bkg:
        vbkg = ROOT.std.vector("float")(n_crystals)
        tout.Branch("BkgEdeps", vbkg)

    step = max(1, n_events // 10)
    t0 = time.time()
    pos = 0
    for i in range(n_events):
        tin.GetEntry(i)
        start = pos
        pos += int(counts[i])
        n_buffer[0] = pos - start
        if store_bkg:
            np.asarray(vbkg)[:] = 0.0
        if pos > start:
            index = mods[start:pos]
            deposit = energies[start:pos]
            np.add.at(np.asarray(vin), index, deposit)
            if store_bkg:
                np.add.at(np.asarray(vbkg), index, deposit)
        tout.Fill()
        if (i + 1) % step == 0 and i + 1 < n_events:
            done = i + 1
            rate = done / (time.time() - t0)
            print(f"    {done:>9d}/{n_events} events  ({rate:,.0f} evt/s)", flush=True)

    tout.Write()

    summary = {
        "source": str(input_path),
        "output": str(output_path),
        "n_events": n_events,
        "n_events_in_source": n_total,
        "n_crystals": n_crystals,
        "n_capture_per_event": float(n_capture),
        "n_hits_injected": n_hits,
        "total_energy_injected_mev": float(energies.astype(np.float64).sum()),
        "mean_energy_per_hit_mev": float(energies.astype(np.float64).mean()) if n_hits else 0.0,
        "seed": int(seed),
        "stream": int(stream),
        "store_bkg_branch": bool(store_bkg),
    }
    if metadata is not None:
        meta = dict(metadata)
        meta.update(summary)
        obj = ROOT.TNamed("bkg_model_neutron", json.dumps(meta, indent=2, ensure_ascii=False))
        obj.Write()
    fout.Close()
    fin.Close()

    elapsed = time.time() - t0
    print(f"  wrote {output_path}  ({n_hits} captures, {elapsed:.1f} s)")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS),
                        help="dataset suffixes, e.g. ep_iso gamma_iso mu_rest")
    parser.add_argument("--data-dir", type=Path, default=default_data_dir(),
                        help="root data directory (default: $ECAL_CLUSTERING_DATA_DIR or <repo>/data)")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="output directory (default: <data-dir>/processed)")
    parser.add_argument("--suffix", default="_bkg_n",
                        help="output filename suffix (default: _bkg_n)")
    parser.add_argument("--seed", type=int, default=20261116,
                        help="random seed (default: 20261116; independent of the gamma run)")
    parser.add_argument("--max-events", type=int, default=0,
                        help="process only the first N events (0 = all)")
    parser.add_argument("--n-capture", type=float, default=None,
                        help="override the estimated captures per event (default: from neutron_model)")
    parser.add_argument("--store-bkg", action="store_true",
                        help="also store a BkgEdeps branch with the background-only vector")
    args = parser.parse_args()

    processed = args.data_dir / "processed"
    out_dir = args.out_dir if args.out_dir is not None else processed

    result = neutron_model.estimate()
    n_capture = result["rate"]["n_capture_per_event"]
    if args.n_capture is not None:
        n_capture = float(args.n_capture)
    energy = neutron_model.CaptureEnergy()

    print(f"data dir   : {processed}")
    print(f"out dir    : {out_dir}")
    print(f"n_capture  : {n_capture:.6f} captures per event (P_cap = 1, all thermal)")
    print(f"capture    : lambda_c = {energy.lam_c_cm:.3f} cm, "
          f"mean deposit {energy.mean_mev:.4f} MeV "
          f"(range {energy.summary()['energy_range_mev'][0]:.2f}-"
          f"{energy.summary()['energy_range_mev'][1]:.2f} MeV)")
    print(f"seed       : {args.seed}")

    summaries = []
    for name in args.datasets:
        input_path = processed / f"edep_of_each_evt_{name}.root"
        output_path = out_dir / f"edep_of_each_evt_{name}{args.suffix}.root"
        stream = zlib.crc32(name.encode("utf-8"))
        print(f"\n[{name}] {input_path.name} -> {output_path.name}")

        summary = inject_file(
            input_path=input_path,
            output_path=output_path,
            n_capture=n_capture,
            energy=energy,
            seed=args.seed,
            stream=stream,
            max_events=args.max_events,
            store_bkg=args.store_bkg,
            metadata={"estimate": result, "model": "pytrial/bkg_est_add/neutron/neutron_model.py"},
        )
        summaries.append(summary)

    print("\nsummary")
    for s in summaries:
        print(f"  {Path(s['output']).name:<45s} events={s['n_events']:>8d}  "
              f"captures={s['n_hits_injected']:>6d}  "
              f"mean E={s['mean_energy_per_hit_mev']:.3f} MeV")


if __name__ == "__main__":
    main()
