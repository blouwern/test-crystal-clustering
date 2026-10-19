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


def energy_metrics(dataset, n_events, cluster_path, eth_wave):
    """Energy budget of the leading cluster / of all clusters.

    In addition to the aggregate containment/capture fractions this records how
    much energy sits in crystals that are above ``Eth_wave`` (i.e. physically
    part of the shower) but were never absorbed into any cluster.  That
    "missed lit energy" is the quantity the wavefront threshold trades against
    over-splitting, so it makes the threshold scan physically interpretable.
    """
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
    lit_total = lit_unclustered = 0.0
    n_all_one = 0
    for evt in range(n_events):
        tree.GetEntry(evt)
        arr = np.array(list(edeps), dtype=np.float64)
        e_tot = float(arr.sum())
        if e_tot <= 0.0:
            continue
        e_max = e_in = e_lit = e_lit_out = 0.0
        claimed = np.zeros(arr.size, dtype=bool)
        for mods in clusters.get(evt, ()):
            e_cl = float(arr[mods].sum())
            e_in += e_cl
            if e_cl > e_max:
                e_max = e_cl
            claimed[mods] = True
        lit_mask = arr >= eth_wave
        if lit_mask.any():
            e_lit = float(arr[lit_mask].sum())
            e_lit_out = float(arr[lit_mask & ~claimed].sum())
        total += e_tot
        leading += e_max
        captured += e_in
        lit_total += e_lit
        lit_unclustered += e_lit_out
        containments.append(e_max / e_tot)
        if e_in >= 0.999 * e_tot:
            n_all_one += 1
    fin.Close()

    if not containments:
        return {}
    c = np.array(containments)
    return {
        "energy_containment_mean": float(c.mean()),
        "energy_containment_median": float(np.median(c)),
        "energy_captured_fraction": float(captured / total),
        "energy_split_fraction": float((captured - leading) / total),
        "energy_unclustered_fraction": float(1.0 - captured / total),
        # energy carried by crystals above Eth_wave that no cluster absorbed
        "lit_energy_total_mev": float(lit_total),
        "lit_energy_missed_fraction": (float(lit_unclustered / lit_total)
                                       if lit_total > 0 else float("nan")),
        "frac_events_containment_ge_90": float(np.mean(c >= 0.90)),
        # NOTE: the historical column name says "one_cluster", but the
        # implemented quantity is ">=99.9 % of the energy captured by any
        # cluster(s)" with the denominator restricted to events with e_tot>0.
        # For the single-particle samples the strict one-cluster version
        # differs by <0.2 pp because over-splitting is rare.
        "frac_events_all_energy_one_cluster": float(n_all_one / len(containments)),
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
    row.update(energy_metrics(dataset, n, out, eth_wave))

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
