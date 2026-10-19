"""Configuration for the two-threshold ECAL clustering sweep.

The clustering takes a *seed* threshold and a (looser) *wavefront* threshold:

  * a crystal becomes a seed only if ``E >= Eth_seed`` **and** at least
    ``Nth`` of its neighbours are above ``Eth_wave``;
  * any crystal with ``E >= Eth_wave`` joins the wavefront.

Ground truth for the single-particle samples is exactly one primary particle.

Grid rationale (2026-10-16 refinement)
-------------------------------------
``Eth_seed`` spans the same coarse range as the first sweep (the metric depends
only weakly on it), but ``Eth_wave`` is sampled far more finely between 0.1 and
2 MeV — that is where the physically interesting transition happens: the
wavefront either absorbs the whole connected shower core (1 cluster) or breaks
it into pieces (>= 2 clusters).  ``Nth`` is trimmed to ``0, 1, 2`` because the
first sweep showed ``Nth = 3`` destroys ~33 % of the acceptance for no gain.

Grid rationale (2026-10-16-refined follow-up)
--------------------------------------------
Two open questions from that campaign are addressed by overriding the grid:

  * all three ``Nth`` optima sat on the lower wave-threshold edge, so 0.02 and
    0.05 MeV are added to locate the plateau and to check whether sub-0.1 MeV
    noise seeding starts to over-split again;
  * the seed axis was flat from 7-15 MeV, so 6, 8 and 9 MeV are added between
    5 and 10 MeV (relevant for the softer ``mu_rest`` spectrum).

Both axes are overridable so the same code reproduces every campaign:

    ECAL_SWEEP_SEEDS=5,6,7,8,9,10,12,15,20
    ECAL_SWEEP_WAVES=0.02,0.05,0.1,0.2,0.3,0.5,0.75,1.0,1.5,2.0,3.0,5.0,6.0,8.0,10.0,12.0
    ECAL_SWEEP_RESULTS_DATE=26-10-16-refined   # deliverable folder
    ECAL_SWEEP_TAG=refined                     # working-area suffix (cache)
"""

import os
from pathlib import Path

SWEEP_DIR = Path(__file__).resolve().parent
REPO_ROOT = SWEEP_DIR.parents[1]

DATA_DIR = Path(os.environ.get("ECAL_CLUSTERING_DATA_DIR") or (REPO_ROOT / "data"))
TAG = os.environ.get("ECAL_SWEEP_TAG", "")
PROCESSED_DIR = DATA_DIR / "processed"
UTILITIES_DIR = DATA_DIR / "utilities"

# regenerable working area (subsets + cluster trees + per-point metrics)
WORK_DIR = SWEEP_DIR / "subsets"
SUBSET_DIR = WORK_DIR
CLUSTER_DIR = WORK_DIR / f"clusters_{TAG}" if TAG else WORK_DIR / "clusters"
TASK_DIR = WORK_DIR / f"tasks_{TAG}" if TAG else WORK_DIR / "tasks"

# deliverables for this run (date identifies the campaign; override to rebuild
# an older one, e.g. ECAL_SWEEP_RESULTS_DATE=26-09-23).  TAG only changes the
# regenerable working directories, never the deliverable folder.
RESULTS_DATE = os.environ.get("ECAL_SWEEP_RESULTS_DATE", "26-10-16")
RESULTS_DIR = SWEEP_DIR.parent / "results" / "sweep" / RESULTS_DATE
SUMMARY_CSV = RESULTS_DIR / "sweep_summary.csv"

DATASETS = ["ep_iso", "gamma_iso", "mu_rest"]
EXPECTED_LAMBDA = 1.0


def _float_list(env, default):
    raw = os.environ.get(env)
    if not raw:
        return list(default)
    return [float(x) for x in raw.split(",") if x.strip()]


# --- grid ------------------------------------------------------------------ #
# seed threshold must stay strictly above the wave threshold
ETH_SEED_VALUES = _float_list(
    "ECAL_SWEEP_SEEDS", [5.0, 7.0, 10.0, 12.0, 15.0, 20.0])
ETH_WAVE_VALUES = _float_list(
    "ECAL_SWEEP_WAVES",
    [0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0],
)
NTH_VALUES = [int(x) for x in os.environ.get("ECAL_SWEEP_NTHS", "0,1,2").split(",")]
N_NEIGHBOR_THRESHOLD = NTH_VALUES[0]   # kept for backwards compatibility


def valid_pairs():
    """(seed, wave) pairs with seed > wave."""
    return [(s, w) for s in ETH_SEED_VALUES for w in ETH_WAVE_VALUES if s > w]


N_EVENTS = int(os.environ.get("ECAL_SWEEP_N_EVENTS", "50000"))
N_WORKERS = int(os.environ.get("ECAL_SWEEP_WORKERS", "18"))
TASK_TIMEOUT = int(os.environ.get("ECAL_SWEEP_TASK_TIMEOUT", "1800"))

# topology consumed by reconstruction(); the shipped *_added.root file does not
# match the MC geometry (see data/utilities/build_neighbor_info_from_geometry.py)
NEIGHBOR_FILE = os.environ.get(
    "ECAL_SWEEP_NEIGHBOR_FILE", "ecal_neighbor_info.root")
