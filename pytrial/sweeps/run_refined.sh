#!/usr/bin/env bash
# Refined follow-up sweep (2026-10-16-refined).
#
# Addresses the two open questions of the 26-10-16 campaign:
#   * Eth_wave below 0.1 MeV (0.02 / 0.05) to locate the acceptance plateau and
#     to check whether sub-0.1 MeV noise seeding over-splits the shower again;
#   * Eth_seed densified between 5 and 10 MeV (6, 8, 9) for the softer mu_rest
#     spectrum.
#
# Usage:
#   bash sweeps/run_refined.sh              # sweep + analysis (cache aware)
#   bash sweeps/run_refined.sh --force      # recompute every point
#   bash sweeps/run_refined.sh --csv-only   # rebuild CSV + figures from cache
#
# Deliverables land in results/sweep/26-10-16-refined/.
set -euo pipefail
cd "$(dirname "$0")/.."

export ECAL_SWEEP_RESULTS_DATE="${ECAL_SWEEP_RESULTS_DATE:-26-10-16-refined}"
export ECAL_SWEEP_TAG="${ECAL_SWEEP_TAG:-refined}"
export ECAL_SWEEP_SEEDS="${ECAL_SWEEP_SEEDS:-5,6,7,8,9,10,12,15,20}"
export ECAL_SWEEP_WAVES="${ECAL_SWEEP_WAVES:-0.02,0.05,0.1,0.2,0.3,0.5,0.75,1.0,1.5,2.0,3.0,5.0,6.0,8.0,10.0,12.0}"
export ECAL_SWEEP_NTHS="${ECAL_SWEEP_NTHS:-0,1,2}"
export ECAL_SWEEP_N_EVENTS="${ECAL_SWEEP_N_EVENTS:-50000}"
export ECAL_SWEEP_NEIGHBOR_FILE="${ECAL_SWEEP_NEIGHBOR_FILE:-ecal_neighbor_info.root}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/.mplconfig}"

python3 sweeps/run_sweep.py "$@"
python3 sweeps/analyze.py
