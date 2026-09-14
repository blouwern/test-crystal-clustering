"""Configuration for the ECAL clustering parameter sweep.

The sweep runs on the single-particle Monte-Carlo samples that are available in
``data/processed`` (e+ / gamma / mu-).  Ground truth for those samples is exactly
one primary particle, hence ``EXPECTED_LAMBDA = 1``.
"""

import os
from pathlib import Path

SWEEP_DIR = Path(__file__).resolve().parent
REPO_ROOT = SWEEP_DIR.parents[1]

# Data root: prefer the environment variable, otherwise <repo>/data.
DATA_DIR = Path(os.environ.get("ECAL_CLUSTERING_DATA_DIR") or (REPO_ROOT / "data"))
PROCESSED_DIR = DATA_DIR / "processed"
UTILITIES_DIR = DATA_DIR / "utilities"

# Working area for the sweep (subsets + cluster trees + per-point metrics).
# These are regenerable intermediates and deliberately kept out of results/.
WORK_DIR = SWEEP_DIR / "subsets"
SUBSET_DIR = WORK_DIR
CLUSTER_DIR = WORK_DIR / "clusters"
TASK_DIR = WORK_DIR / "tasks"

# Deliverables requested by the user live under pytrial/results/.
RESULTS_DIR = SWEEP_DIR.parent / "results" / "sweep"
SUMMARY_CSV = RESULTS_DIR / "sweep_summary.csv"

# Single-particle samples available after the raw-data processing step.
DATASETS = ["ep_iso", "gamma_iso", "mu_rest"]
EXPECTED_LAMBDA = 1.0

# Fine grid from requirement.md NO.6 (36 combinations).
ETH_VALUES = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 20.0]
NTH_VALUES = [1.0, 2.0, 3.0]

# Events per dataset used for the sweep (out of 1,000,000 available).
N_EVENTS = int(os.environ.get("ECAL_SWEEP_N_EVENTS", "50000"))
N_WORKERS = int(os.environ.get("ECAL_SWEEP_WORKERS", "18"))

# Nearest-neighbour topology actually consumed by pytrial/reconstruction.py.
NEIGHBOR_FILE = "ecal_neighbor_info_added.root"
