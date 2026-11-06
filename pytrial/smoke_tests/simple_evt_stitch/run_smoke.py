"""Stage 1 (Scheme A) smoke test: run the current reconstruction and count.

The stitched sample is produced by ``stitch_smoke.py``.  This script only
exercises the I/O + cluster-count path; per-primary truth, efficiency, purity,
merge and split are Stage 2/3 deliverables.
"""

from __future__ import annotations

import array
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PYTRIAL = HERE.parents[1]
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(PYTRIAL))

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(REPO_ROOT / ".mplconfig"))

import ROOT  # noqa: E402
from ROOT import vector  # noqa: E402

import reconstruction  # noqa: E402
from evaluation import counts  # noqa: E402

STITCHED = HERE / "data" / "stitched_ep_iso_n2_2k.root"
META = HERE / "meta.json"
WORK = HERE / "work"
RESULTS = HERE / "results"
CLUSTERS = RESULTS / "clusters_ep_iso_n2_2k_nth0_seed7_wave0.05.root"

N_EVENTS = 2_000
N_TRUTH = 2
ETH_SEED = 7.0
ETH_WAVE = 0.05
NTH = 0


def prepare_workdir():
    """Create the ``processed/`` + ``utilities/`` layout reconstruction expects."""
    proc = WORK / "processed"
    util = WORK / "utilities"
    proc.mkdir(parents=True, exist_ok=True)
    util.mkdir(parents=True, exist_ok=True)

    link = proc / STITCHED.name
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(STITCHED.resolve())

    topo_src = REPO_ROOT / "data" / "utilities" / "ecal_neighbor_info.root"
    topo_link = util / topo_src.name
    if topo_link.is_symlink() or topo_link.exists():
        topo_link.unlink()
    topo_link.symlink_to(topo_src.resolve())
    return WORK


def inspect_stitched():
    """Check NTrack / Edeps length / SourceIdx length on the stitched sample."""
    ntrack_ok = True
    length_ok = True
    source_ok = True
    f = ROOT.TFile.Open(str(STITCHED), "READ")
    tree = f.Get("EdepOfEachEvt")
    n = tree.GetEntries()
    ntrack = array.array("i", [0])
    edeps = vector("float")()
    source = vector("int")()
    tree.SetBranchAddress("NTrack", ntrack)
    tree.SetBranchAddress("Edeps", edeps)
    tree.SetBranchAddress("SourceIdx", source)
    for i in range(n):
        tree.GetEntry(i)
        ntrack_ok &= int(ntrack[0]) == N_TRUTH
        length_ok &= edeps.size() == 568
        source_ok &= source.size() == N_TRUTH
    f.Close()
    return {
        "n_entries": int(n),
        "ntrack_all_2": bool(ntrack_ok),
        "edeps_length_568": bool(length_ok),
        "sourceidx_length_2": bool(source_ok),
    }


def write_plot(cluster_counts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    max_k = max(int(cluster_counts.max()), N_TRUTH)
    hist = np.bincount(cluster_counts, minlength=max_k + 1)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.bar(range(len(hist)), hist, color="steelblue", edgecolor="black")
    ax.axvline(N_TRUTH, color="red", linestyle="--", label=f"truth N={N_TRUTH}")
    ax.set_xlabel("clusters per event")
    ax.set_ylabel("events")
    ax.set_title(f"Stage 1 smoke test: {N_EVENTS} stitched events (N={N_TRUTH})")
    ax.set_xticks(range(len(hist)))
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    out = RESULTS / "count_distribution.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main():
    if not STITCHED.exists():
        raise SystemExit(f"missing {STITCHED}; run stitch_smoke.py first")
    if not META.exists():
        raise SystemExit(f"missing {META}; run stitch_smoke.py first")

    RESULTS.mkdir(parents=True, exist_ok=True)
    prepare_workdir()

    print(f"reconstructing {STITCHED.name} with Nth={NTH}, "
          f"Eth_seed={ETH_SEED}, Eth_wave={ETH_WAVE}")
    reconstruction.reconstruction(
        filename_edep_of_each_evt=STITCHED.name,
        energy_threshold_seed=ETH_SEED,
        energy_threshold_wave=ETH_WAVE,
        n_neighbor_threshold=NTH,
        data_dir=str(WORK),
        output_path=str(CLUSTERS),
    )

    cluster_counts = counts.read_cluster_counts(CLUSTERS, n_events=N_EVENTS)
    metrics = counts.count_metrics(cluster_counts, n_truth=N_TRUTH)

    (RESULTS / "metrics.json").write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False) + "\n")
    flat = {k: v for k, v in metrics.items() if k != "histogram"}
    with open(RESULTS / "metrics.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(flat))
        writer.writeheader()
        writer.writerow(flat)
    plot_path = write_plot(cluster_counts)

    meta = json.loads(META.read_text())
    stitched_info = inspect_stitched()
    checks = {
        "stage": "1",
        "scheme": "A",
        "stitched_info": stitched_info,
        "energy_check": meta["energy_check"],
        "cluster_output_exists": CLUSTERS.exists(),
        "n_events_read": int(len(cluster_counts)),
        "all_events_accounted": int(len(cluster_counts)) == N_EVENTS,
        "metrics_finite": bool(all(
            np.isfinite(v) for k, v in metrics.items() if k != "histogram")),
        "reconstruction": {
            "Nth": NTH,
            "Eth_seed_mev": ETH_SEED,
            "Eth_wave_mev": ETH_WAVE,
            "cluster_file": str(CLUSTERS.relative_to(REPO_ROOT)),
        },
    }
    checks["passed"] = bool(
        stitched_info["ntrack_all_2"]
        and stitched_info["edeps_length_568"]
        and stitched_info["sourceidx_length_2"]
        and meta["energy_check"]["passed"]
        and checks["cluster_output_exists"]
        and checks["all_events_accounted"]
        and checks["metrics_finite"]
    )
    (RESULTS / "checks.json").write_text(
        json.dumps(checks, indent=2, ensure_ascii=False) + "\n")

    print("\n=== Stage 1 smoke test ===")
    print(f"events: {metrics['n_events']}  truth N: {metrics['n_truth']}")
    print(f"exact-{N_TRUTH}: {metrics['exact_rate']:.4f}")
    print(f"under:    {metrics['under_rate']:.4f}")
    print(f"over:     {metrics['over_rate']:.4f}")
    print(f"MAE:      {metrics['MAE']:.4f}")
    print(f"RMSE:     {metrics['RMSE']:.4f}")
    print(f"bias:     {metrics['bias']:+.4f}")
    print(f"mean clusters: {metrics['mean_clusters']:.4f}")
    print(f"checks passed: {checks['passed']}")
    print(f"metrics: {RESULTS / 'metrics.json'}")
    print(f"plot:    {plot_path}")


if __name__ == "__main__":
    main()
