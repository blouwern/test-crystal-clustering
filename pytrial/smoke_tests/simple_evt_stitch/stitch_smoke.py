"""Stage 1 (Scheme A) smoke test: count-only N=2 stitching.

Reads the ep_iso single-particle sample, keeps a positive-energy pool, sums two
events per stitched event, and writes a self-contained ROOT file with
``Edeps``, ``NTrack`` and the source indices.  No per-primary truth is stored:
that is the Stage 2 (Scheme B) deliverable.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PYTRIAL = HERE.parents[1]
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(PYTRIAL))

import stitching  # noqa: E402

INPUT = REPO_ROOT / "data" / "processed" / "edep_of_each_evt_ep_iso.root"
OUTPUT = HERE / "data" / "stitched_ep_iso_n2_2k.root"
META = HERE / "meta.json"

POOL_EVENTS = 20_000
N_PRIMARY = 2
N_EVENTS = 2_000
SEED = 20261018


def main():
    rng = np.random.default_rng(SEED)

    print(f"loading first {POOL_EVENTS} events from {INPUT}")
    pool_all = stitching.load_edeps(INPUT, n_events=POOL_EVENTS)
    mask = stitching.positive_energy_mask(pool_all)
    pool = pool_all[mask]
    print(f"positive-energy pool: {len(pool)}/{len(pool_all)} events")

    stitched, ntrack, source_idx = stitching.stitch_count_only(
        pool, N_PRIMARY, N_EVENTS, rng)

    # Energy is additive by construction; verify it explicitly.
    source_sum = pool[source_idx[:, 0]].sum(axis=1) + pool[source_idx[:, 1]].sum(axis=1)
    stitched_sum = stitched.sum(axis=1)
    max_abs_diff = float(np.max(np.abs(source_sum - stitched_sum)))

    stitching.write_stitched(OUTPUT, stitched, ntrack, source_idx)
    print(f"wrote {OUTPUT} ({len(stitched)} events, N={N_PRIMARY})")

    meta = {
        "stage": "1",
        "scheme": "A",
        "input": str(INPUT.relative_to(REPO_ROOT)),
        "output": str(OUTPUT.relative_to(REPO_ROOT)),
        "pool_events_read": POOL_EVENTS,
        "pool_events_positive": int(len(pool)),
        "n_primary": N_PRIMARY,
        "n_events": N_EVENTS,
        "seed": SEED,
        "rng": "numpy.random.default_rng(PCG64)",
        "energy_check": {
            "max_abs_diff_mev": max_abs_diff,
            "tolerance_mev": 1e-3,
            "passed": max_abs_diff < 1e-3,
        },
        "note": "count-only stitching; no per-primary truth (Stage 2)",
    }
    META.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    print(f"wrote {META}")
    print(f"energy conservation: max |diff| = {max_abs_diff:.3e} MeV")


if __name__ == "__main__":
    main()
