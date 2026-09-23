"""Run one sweep point: reconstruct + score a (dataset, Eth_seed, Eth_wave).

Usage:  python sweeps/run_task.py <dataset> <Eth_seed> <Eth_wave>
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
    CLUSTER_DIR, EXPECTED_LAMBDA, NEIGHBOR_FILE, N_NEIGHBOR_THRESHOLD,
    SUBSET_DIR, TASK_DIR,
)


def chisquare_vs_poisson(counts, lam):
    """Chi-squared GOF against Poisson(lam), merging bins to expected >= 5."""
    n = counts.size
    kmax = int(counts.max())
    k = np.arange(0, kmax + 1)
    obs = np.array([np.sum(counts == i) for i in k], dtype=float)
    exp = n * stats.poisson.pmf(k, lam)
    obs = np.append(obs, float(np.sum(counts > kmax)))
    exp = np.append(exp, n * (1.0 - stats.poisson.cdf(kmax, lam)))

    mo, me = [], []
    po = pe = 0.0
    for o, e in zip(obs, exp):
        po += o
        pe += e
        if pe >= 5.0:
            mo.append(po)
            me.append(pe)
            po = pe = 0.0
    if pe > 0 and me:
        mo[-1] += po
        me[-1] += pe
    if len(mo) < 2 or min(me) <= 0:
        return float("nan"), float("nan"), 0
    chi2, p = stats.chisquare(mo, me, ddof=0)
    return float(chi2), float(p), len(mo) - 1


def energy_metrics(dataset, n_events, cluster_path):
    """Energy budget of the leading cluster / of all clusters."""
    clusters = defaultdict(list)
    with ROOT.TFile(str(cluster_path), "READ") as f:
        for entry in f.Get("ECALClusters"):
            clusters[int(entry.EvtID)].append([int(m) for m in entry.ModIDList])

    fin = ROOT.TFile.Open(str(SUBSET_DIR / "processed" / f"subset_{dataset}.root"), "READ")
    tree = fin.Get("EdepOfEachEvt")
    edeps = vector("float")()
    tree.SetBranchAddress("Edeps", edeps)

    containments = []
    total = leading = captured = 0.0
    for evt in range(n_events):
        tree.GetEntry(evt)
        arr = np.array(list(edeps), dtype=np.float64)
        e_tot = float(arr.sum())
        if e_tot <= 0.0:
            continue
        e_max = e_in = 0.0
        for mods in clusters.get(evt, ()):
            e_cl = float(arr[mods].sum())
            e_in += e_cl
            e_max = max(e_max, e_cl)
        total += e_tot
        leading += e_max
        captured += e_in
        containments.append(e_max / e_tot)
    fin.Close()

    if not containments:
        return {}
    c = np.array(containments)
    return {
        "energy_containment_mean": float(c.mean()),
        "energy_captured_fraction": float(captured / total),
        "energy_split_fraction": float((captured - leading) / total),
        "energy_unclustered_fraction": float(1.0 - captured / total),
        "frac_events_containment_ge_90": float(np.mean(c >= 0.90)),
    }


def subset_n_events(dataset):
    """Number of events in the subset file (needed when the tree is empty)."""
    f = ROOT.TFile.Open(str(SUBSET_DIR / "processed" / f"subset_{dataset}.root"), "READ")
    n = f.Get("EdepOfEachEvt").GetEntries()
    f.Close()
    return int(n)


def run_point(dataset, eth_seed, eth_wave, nth=0):
    CLUSTER_DIR.mkdir(parents=True, exist_ok=True)
    TASK_DIR.mkdir(parents=True, exist_ok=True)

    tag = f"{dataset}_seed{eth_seed:g}_wave{eth_wave:g}_nth{nth:g}"
    out = CLUSTER_DIR / f"{tag}.root"

    t0 = time.time()
    reconstruction.reconstruction(
        f"subset_{dataset}.root", eth_seed, eth_wave, nth,
        neighbor_file_name=NEIGHBOR_FILE,
        data_dir=str(SUBSET_DIR), output_path=str(out),
    )
    recon_seconds = time.time() - t0

    n_expected = subset_n_events(dataset)
    counts = scoring._read_cluster_counts(out, n_events=n_expected).astype(np.int64)
    n = int(counts.size)
    lam = EXPECTED_LAMBDA
    diff = counts.astype(np.float64) - lam
    chi2, chi2_p, ndf = chisquare_vs_poisson(counts, lam)
    ks_d, ks_p = stats.kstest(counts, "poisson", args=(lam,))

    row = {
        "dataset": dataset,
        "Eth_seed": eth_seed,
        "Eth_wave": eth_wave,
        "Nth": nth,
        "neighbor_file": NEIGHBOR_FILE,
        "N_events": n,
        "mean_clusters": float(counts.mean()),
        "std_clusters": float(counts.std(ddof=1)) if n > 1 else 0.0,
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
        "chi2_ndf": (chi2 / ndf) if ndf > 0 else float("nan"),
        "ks_D": float(ks_d),
        "recon_seconds": recon_seconds,
        "empty_tree": bool(counts.sum() == 0),
    }
    row.update(energy_metrics(dataset, n, out))

    with open(TASK_DIR / f"{tag}.json", "w") as f:
        json.dump(row, f, indent=1, sort_keys=True)
    print(f"{tag:32s} mean={row['mean_clusters']:.3f} "
          f"exact={row['exact_match_rate']:.3f} "
          f"cont={row.get('energy_containment_mean', float('nan')):.3f} "
          f"({recon_seconds:.1f}s)")
    return row


if __name__ == "__main__":
    run_point(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]),
              int(float(sys.argv[4])) if len(sys.argv) > 4 else N_NEIGHBOR_THRESHOLD)
