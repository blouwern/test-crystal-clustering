"""Run one sweep point: reconstruct + score a single (dataset, Eth, Nth).

Usage:  python sweeps/run_task.py <dataset> <Eth> <Nth>

Writes the cluster tree to ``sweeps/subsets/clusters/<tag>.root`` and a JSON
metric record to ``sweeps/subsets/tasks/<tag>.json``.
"""

import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import ROOT
from ROOT import vector
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import reconstruction  # noqa: E402
import scoring  # noqa: E402
from sweeps.config import (  # noqa: E402
    CLUSTER_DIR, EXPECTED_LAMBDA, SUBSET_DIR, TASK_DIR,
)


def chisquare_vs_poisson(counts, lam):
    """Chi-squared GOF against Poisson(lam) with full probability mass.

    Bins 0..kmax plus a tail bin are merged until every bin has an expected
    count of at least 5, so no probability mass is dropped.
    """
    n = counts.size
    kmax = int(counts.max())
    k = np.arange(0, kmax + 1)
    obs = np.array([np.sum(counts == i) for i in k], dtype=float)
    exp = n * stats.poisson.pmf(k, lam)
    obs = np.append(obs, float(np.sum(counts > kmax)))
    exp = np.append(exp, n * (1.0 - stats.poisson.cdf(kmax, lam)))

    merged_obs, merged_exp = [], []
    pool_o = pool_e = 0.0
    for o, e in zip(obs, exp):
        pool_o += o
        pool_e += e
        if pool_e >= 5.0:
            merged_obs.append(pool_o)
            merged_exp.append(pool_e)
            pool_o = pool_e = 0.0
    if pool_e > 0 and merged_exp:
        merged_obs[-1] += pool_o
        merged_exp[-1] += pool_e

    if len(merged_obs) < 2 or min(merged_exp) <= 0:
        return float("nan"), float("nan"), 0
    chi2, p_value = stats.chisquare(merged_obs, merged_exp, ddof=0)
    return float(chi2), float(p_value), len(merged_obs) - 1


def energy_metrics(dataset, n_events, cluster_path):
    """Energy containment of the leading cluster (single-particle truth)."""
    clusters = defaultdict(list)
    with ROOT.TFile(str(cluster_path), "READ") as f:
        tree = f.Get("ECALClusters")
        for entry in tree:
            clusters[int(entry.EvtID)].append([int(m) for m in entry.ModIDList])

    subset = SUBSET_DIR / "processed" / f"subset_{dataset}.root"
    fin = ROOT.TFile.Open(str(subset), "READ")
    tree = fin.Get("EdepOfEachEvt")
    edeps = vector("float")()
    tree.SetBranchAddress("Edeps", edeps)

    containments = []
    total_energy = 0.0
    leading_energy = 0.0
    captured_energy = 0.0
    for evt in range(n_events):
        tree.GetEntry(evt)
        arr = np.array(list(edeps), dtype=np.float64)
        e_tot = float(arr.sum())
        if e_tot <= 0.0:
            continue
        e_max = 0.0
        e_in_clusters = 0.0
        for mods in clusters.get(evt, ()):
            e_cl = float(arr[mods].sum())
            e_in_clusters += e_cl
            if e_cl > e_max:
                e_max = e_cl
        total_energy += e_tot
        leading_energy += e_max
        captured_energy += e_in_clusters
        containments.append(e_max / e_tot)
    fin.Close()

    containments = np.array(containments)
    if containments.size == 0:
        return {}
    return {
        # fraction of the event energy sitting in the *leading* cluster
        "energy_containment_mean": float(containments.mean()),
        "energy_containment_median": float(np.median(containments)),
        "energy_containment_std": float(containments.std(ddof=1)),
        "frac_events_containment_ge_90": float(np.mean(containments >= 0.90)),
        "frac_events_containment_ge_99": float(np.mean(containments >= 0.99)),
        # fraction of the event energy sitting in *some* cluster (threshold loss)
        "energy_captured_fraction": float(captured_energy / total_energy),
        # energy in clusters other than the leading one (splitting loss)
        "energy_split_fraction": float((captured_energy - leading_energy) / total_energy),
        # energy in no cluster at all (sub-threshold crystals)
        "energy_unclustered_fraction": float(1.0 - captured_energy / total_energy),
    }


def run_point(dataset, eth, nth):
    CLUSTER_DIR.mkdir(parents=True, exist_ok=True)
    TASK_DIR.mkdir(parents=True, exist_ok=True)

    tag = f"{dataset}_Eth{eth:g}_Nth{nth:g}"
    out = CLUSTER_DIR / f"{tag}.root"

    t0 = time.time()
    reconstruction.reconstruction(
        f"subset_{dataset}.root", eth, nth,
        data_dir=str(SUBSET_DIR), output_path=str(out),
    )
    recon_seconds = time.time() - t0

    counts = scoring._read_cluster_counts(out).astype(np.int64)
    n = int(counts.size)
    lam = EXPECTED_LAMBDA
    diff = counts.astype(np.float64) - lam

    chi2, chi2_p, ndf = chisquare_vs_poisson(counts, lam)
    ks_d, ks_p = stats.kstest(counts, "poisson", args=(lam,))

    row = {
        "dataset": dataset,
        "Eth": eth,
        "Nth": nth,
        "N_events": n,
        "mean_clusters": float(counts.mean()),
        "std_clusters": float(counts.std(ddof=1)),
        "median_clusters": float(np.median(counts)),
        "max_clusters": int(counts.max()),
        "frac_0": float(np.mean(counts == 0)),
        "frac_1": float(np.mean(counts == 1)),
        "frac_2": float(np.mean(counts == 2)),
        "frac_3plus": float(np.mean(counts >= 3)),
        "exact_match_rate": float(np.mean(diff == 0)),
        "under_cluster_rate": float(np.mean(diff < 0)),
        "over_cluster_rate": float(np.mean(diff > 0)),
        "bias": float(counts.mean() - lam),
        "MAE": float(np.mean(np.abs(diff))),
        "RMSE": float(np.sqrt(np.mean(diff ** 2))),
        "chi2_vs_poisson1": chi2,
        "chi2_ndf": (chi2 / ndf) if ndf > 0 else float("nan"),
        "chi2_ndf_dof": ndf,
        "chi2_p": chi2_p,
        "ks_D": float(ks_d),
        "ks_p": float(ks_p),
        "recon_seconds": recon_seconds,
    }
    row.update(energy_metrics(dataset, n, out))

    with open(TASK_DIR / f"{tag}.json", "w") as f:
        json.dump(row, f, indent=1, sort_keys=True)
    print(f"{tag:28s} mean={row['mean_clusters']:.4f} "
          f"exact={row['exact_match_rate']:.4f} MAE={row['MAE']:.4f} "
          f"cont={row.get('energy_containment_mean', float('nan')):.4f} "
          f"({recon_seconds:.1f}s)")
    return row


def main():
    dataset, eth, nth = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    run_point(dataset, eth, nth)


if __name__ == "__main__":
    main()
