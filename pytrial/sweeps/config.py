"""Configuration for the two-threshold ECAL clustering sweep.

The clustering now takes a *seed* threshold and a (looser) *wavefront* threshold:

  * a crystal becomes a seed only if ``E >= Eth_seed`` **and** at least
    ``Nth`` of its neighbours are above ``Eth_wave``;
  * any crystal with ``E >= Eth_wave`` joins the wavefront.

Ground truth for the single-particle samples is exactly one primary particle.
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

# deliverables for this run
RESULTS_DIR = SWEEP_DIR.parent / "results" / "sweep" / "26-09-16" / TAG
SUMMARY_CSV = RESULTS_DIR / "sweep_summary.csv"

DATASETS = ["ep_iso", "gamma_iso", "mu_rest"]
EXPECTED_LAMBDA = 1.0

# --- grid ------------------------------------------------------------------ #
# seed threshold must stay strictly above the wave threshold
ETH_SEED_VALUES = [5.0, 7.0, 10.0, 12.0, 15.0, 20.0]
ETH_WAVE_VALUES = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]
NTH_VALUES = [int(x) for x in os.environ.get("ECAL_SWEEP_NTHS", "0,1,2,3").split(",")]
N_NEIGHBOR_THRESHOLD = NTH_VALUES[0]   # kept for backwards compatibility

def valid_pairs():
    """(seed, wave) pairs with seed > wave."""
    return [(s, w) for s in ETH_SEED_VALUES for w in ETH_WAVE_VALUES if s > w]


N_EVENTS = int(os.environ.get("ECAL_SWEEP_N_EVENTS", "50000"))
N_WORKERS = int(os.environ.get("ECAL_SWEEP_WORKERS", "18"))

# topology consumed by reconstruction(); the shipped *_added.root file does not
# match the MC geometry (see data/utilities/build_neighbor_info_from_geometry.py)
NEIGHBOR_FILE = os.environ.get(
    "ECAL_SWEEP_NEIGHBOR_FILE", "ecal_neighbor_info_geo.root")
