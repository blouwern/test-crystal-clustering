#!/usr/bin/env python3
"""Step 5: inject the estimated gamma background into the processed samples.

For every event the script

1. draws the number of background hits ``n ~ Poisson(n_bkg)``;
2. draws ``n`` crystal IDs uniformly in ``[0, N_crystals)``;
3. draws ``n`` energies from the assumed deposition spectrum
   (normal, mean 5 MeV, variance 2 MeV^2, truncated at zero);
4. adds those energies to the ``Edeps`` vector of the event.

The original files are never modified: the output is written next to them with
the ``_bkg`` suffix (``edep_of_each_evt_ep_iso_bkg.root`` etc.).  Each output
file keeps the ``EdepOfEachEvt`` tree with the ``Edeps`` branch, so it is a
drop-in replacement for the clustering pipeline, and adds an ``NBkg`` branch
with the number of injected hits.  A ``bkg_model`` metadata object records the
parameters, the seed and the model result used for the file.

Usage
-----
    python add_background.py                       # all three samples, full 1e6 events
    python add_background.py --datasets mu_rest    # one sample only
    python add_background.py --max-events 10000    # quick test
    python add_background.py --store-bkg           # also store the background-only vector
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

import bkg_model

ROOT.gROOT.SetBatch(True)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
DEFAULT_DATASETS = ("ep_iso", "gamma_iso", "mu_rest")
TREE_NAME = "EdepOfEachEvt"
BRANCH_NAME = "Edeps"


def default_data_dir() -> Path:
    env = os.environ.get("ECAL_CLUSTERING_DATA_DIR")
    return Path(env) if env else REPO / "data"


def draw_injections(
    n_events: int,
    n_crystals: int,
    n_bkg: float,
    spectrum: bkg_model.DepositionSpectrum,
    seed: int,
    stream: int,
):
    """Draw the background for a whole file (deterministic given the seed).

    Returns ``(counts, mods, energies)`` where ``counts[i]`` is the number of
    background hits of event ``i`` and the flat arrays ``mods``/``energies``
    hold the crystal IDs and deposits in event order.
    """
    rng = np.random.default_rng([int(seed), int(stream)])
    counts = rng.poisson(n_bkg, size=int(n_events)).astype(np.int64)
    n_hits = int(counts.sum())
    mods = rng.integers(0, int(n_crystals), size=n_hits, dtype=np.int64).astype(np.int32)
    energies = spectrum.rvs(n_hits, rng).astype(np.float32)
    return counts, mods, energies


def inject_file(
    input_path: Path,
    output_path: Path,
    n_bkg: float,
    spectrum: bkg_model.DepositionSpectrum,
    seed: int,
    stream: int,
    max_events: int = 0,
    store_bkg: bool = False,
    metadata: dict | None = None,
) -> dict:
    """Add background to one ROOT file and return a small summary."""
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
        n_events, n_crystals, n_bkg, spectrum, seed, stream
    )
    n_hits = int(counts.sum())

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fout = ROOT.TFile.Open(str(output_path), "RECREATE")
    tout = ROOT.TTree(TREE_NAME, "Edep of each event (estimated gamma background added)")
    tout.Branch(BRANCH_NAME, vin)
    nbkg_buffer = np.zeros(1, dtype=np.int32)
    tout.Branch("NBkg", nbkg_buffer, "NBkg/I")
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
        nbkg_buffer[0] = pos - start
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
        "n_bkg_per_event": float(n_bkg),
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
        obj = ROOT.TNamed("bkg_model", json.dumps(meta, indent=2, ensure_ascii=False))
        obj.Write()
    fout.Close()
    fin.Close()

    elapsed = time.time() - t0
    print(f"  wrote {output_path}  ({n_hits} hits, {elapsed:.1f} s)")
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
    parser.add_argument("--suffix", default="_bkg", help="output filename suffix (default: _bkg)")
    parser.add_argument("--seed", type=int, default=20261105, help="random seed (default: 20261105)")
    parser.add_argument("--max-events", type=int, default=0,
                        help="process only the first N events (0 = all)")
    parser.add_argument("--n-bkg", type=float, default=None,
                        help="override the estimated n_bkg (default: from bkg_model)")
    parser.add_argument("--store-bkg", action="store_true",
                        help="also store a BkgEdeps branch with the background-only vector")
    args = parser.parse_args()

    processed = args.data_dir / "processed"
    out_dir = args.out_dir if args.out_dir is not None else processed

    result = bkg_model.estimate()
    n_bkg = result["rates"]["n_bkg_per_event"]
    if args.n_bkg is not None:
        n_bkg = float(args.n_bkg)
    spectrum = bkg_model.DepositionSpectrum(
        result["spectrum"]["nominal_mean_mev"], result["spectrum"]["variance_mev2"]
    )

    print(f"data dir   : {processed}")
    print(f"out dir    : {out_dir}")
    print(f"n_bkg      : {n_bkg:.6f} hits per event "
          f"(V_total = {result['geometry']['n_crystals']} x "
          f"{result['geometry']['volume_cm3']:.1f} cm^3)")
    print(f"spectrum   : mean {spectrum.mean_deposition_mev:.4f} MeV, sigma {spectrum.sigma_mev:.4f} MeV")
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
            n_bkg=n_bkg,
            spectrum=spectrum,
            seed=args.seed,
            stream=stream,
            max_events=args.max_events,
            store_bkg=args.store_bkg,
            metadata={"estimate": result, "model": "pytrial/bkg_est_add/bkg_model.py"},
        )
        summaries.append(summary)

    print("\nsummary")
    for s in summaries:
        print(f"  {Path(s['output']).name:<45s} events={s['n_events']:>8d}  "
              f"hits={s['n_hits_injected']:>6d}  "
              f"mean E={s['mean_energy_per_hit_mev']:.3f} MeV")


if __name__ == "__main__":
    main()
