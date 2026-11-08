"""Run one sweep point: reconstruct + score a (dataset, Eth_seed, Eth_wave).

For background-injected samples (``*_bkg``) the same event index in the
matching signal-only sample is read as well, so every reconstructed cluster can
be split crystal-by-crystal into signal and background energy.  This powers the
background-aware metrics (purity, fake-cluster rate, energy response) and the
per-energy / per-NBkg binned curves written next to the scalar JSON.

Usage:  python sweeps/run_task.py <dataset> <Eth_seed> <Eth_wave> [Nth]
"""

import array
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
    SUBSET_DIR, TASK_DIR, has_bkg, signal_dataset,
)

# signal-energy binning for the efficiency / response curves (MeV)
E_SIG_EDGES = np.array(
    [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 80], dtype=float)
# NBkg bins: 0, 1, 2, 3, >=4
NBKG_VALUES = [0, 1, 2, 3, 4]
# a cluster with < 1 MeV of signal energy is treated as background-like
FAKE_SIGNAL_THRESHOLD = 1.0


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


def _new_stats(nbins):
    """Accumulators for one set of bins (energy bins or NBkg bins)."""
    zeros = lambda: np.zeros(nbins, dtype=float)  # noqa: E731
    return {
        "n": np.zeros(nbins, dtype=np.int64),
        "n_with": np.zeros(nbins, dtype=np.int64),
        "n_exact": np.zeros(nbins, dtype=np.int64),
        "n_over": np.zeros(nbins, dtype=np.int64),
        "sum_clusters": zeros(),
        "sum_fake": zeros(),
        "sum_esig": zeros(),
        "sum_ebkg": zeros(),
        "sum_lead_sig": zeros(),
        "sum_lead_tot": zeros(),
        "sum_lead_bkg": zeros(),
        "sum_sig_in": zeros(),
        "sum_bkg_in": zeros(),
        "sum_bias_tot": zeros(),
        "sum_bias_tot2": zeros(),
        "n_bias_tot": np.zeros(nbins, dtype=np.int64),
        "sum_bias_sig": zeros(),
        "sum_bias_sig2": zeros(),
        "n_bias_sig": np.zeros(nbins, dtype=np.int64),
    }


def _ratio(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    out = np.full(a.shape, np.nan)
    np.divide(a, b, out=out, where=b > 0)
    return out


def _finalize_bins(stats):
    """Turn accumulators into per-bin rates / means (lists for JSON)."""
    n = stats["n"]
    out = {
        "n": n.astype(int).tolist(),
        "efficiency": _ratio(stats["n_with"], n).tolist(),
        "exact": _ratio(stats["n_exact"], n).tolist(),
        "over": _ratio(stats["n_over"], n).tolist(),
        "mean_clusters": _ratio(stats["sum_clusters"], n).tolist(),
        "fake_per_event": _ratio(stats["sum_fake"], n).tolist(),
        "lead_sig_containment": _ratio(
            stats["sum_lead_sig"], stats["sum_esig"]).tolist(),
        "lead_bkg_fraction": _ratio(
            stats["sum_lead_bkg"], stats["sum_lead_tot"]).tolist(),
        "signal_capture": _ratio(
            stats["sum_sig_in"], stats["sum_esig"]).tolist(),
        "bkg_capture": _ratio(
            stats["sum_bkg_in"], stats["sum_ebkg"]).tolist(),
        "signal_purity": _ratio(
            stats["sum_sig_in"],
            stats["sum_sig_in"] + stats["sum_bkg_in"]).tolist(),
    }
    bias_t = _ratio(stats["sum_bias_tot"], stats["n_bias_tot"])
    var_t = _ratio(stats["sum_bias_tot2"], stats["n_bias_tot"]) - bias_t ** 2
    out["bias_total"] = bias_t.tolist()
    out["rms_total"] = np.sqrt(np.maximum(var_t, 0.0)).tolist()
    bias_s = _ratio(stats["sum_bias_sig"], stats["n_bias_sig"])
    var_s = _ratio(stats["sum_bias_sig2"], stats["n_bias_sig"]) - bias_s ** 2
    out["bias_signal"] = bias_s.tolist()
    out["rms_signal"] = np.sqrt(np.maximum(var_s, 0.0)).tolist()
    return out


def _mean(values):
    return float(np.mean(values)) if len(values) else float("nan")


def _rms(values):
    if not len(values):
        return float("nan")
    a = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean((a - a.mean()) ** 2)))


def energy_metrics(dataset, n_events, cluster_path, eth_wave):
    """Energy / background budget of the clusters for one sweep point.

    Returns ``(scalar, bins)``.  The scalar dict keeps the historical energy
    keys and adds background-aware ones.  ``bins`` holds the signal-energy and
    NBkg-binned curves used by the background report.
    """
    is_bkg = has_bkg(dataset)
    sig_name = signal_dataset(dataset)

    clusters = defaultdict(list)
    with ROOT.TFile(str(cluster_path), "READ") as f:
        for entry in f.Get("ECALClusters"):
            clusters[int(entry.EvtID)].append([int(m) for m in entry.ModIDList])

    fin_b = ROOT.TFile.Open(
        str(SUBSET_DIR / "processed" / f"subset_{dataset}.root"), "READ")
    tree_b = fin_b.Get("EdepOfEachEvt")
    edeps_b = vector("float")()
    tree_b.SetBranchAddress("Edeps", edeps_b)
    has_nbkg = any(b.GetName() == "NBkg" for b in tree_b.GetListOfBranches())
    nbkg_vec = array.array("i", [0])
    if has_nbkg:
        tree_b.SetBranchAddress("NBkg", nbkg_vec)

    if is_bkg:
        fin_s = ROOT.TFile.Open(
            str(SUBSET_DIR / "processed" / f"subset_{sig_name}.root"), "READ")
        tree_s = fin_s.Get("EdepOfEachEvt")
        edeps_s = vector("float")()
        tree_s.SetBranchAddress("Edeps", edeps_s)
    else:
        fin_s, tree_s, edeps_s = fin_b, tree_b, edeps_b

    # historical energy accumulators (total energy, includes background)
    containments = []
    total = leading = captured = 0.0
    lit_total = lit_unclustered = 0.0
    n_all_one = 0

    # background-aware accumulators
    nbkg_list, esig_list, ebkg_list = [], [], []
    esig_all = ebkg_all = 0.0
    sig_in_all = bkg_in_all = 0.0
    fake_all = 0
    n_fake_events = 0
    lead_sig_cont, lead_bkg_frac, purity = [], [], []
    bias_tot, bias_sig = [], []
    bins_e = _new_stats(len(E_SIG_EDGES) - 1)
    bins_n = _new_stats(len(NBKG_VALUES))

    for evt in range(n_events):
        tree_b.GetEntry(evt)
        arr_b = np.array(list(edeps_b), dtype=np.float64)
        if is_bkg:
            tree_s.GetEntry(evt)
            arr_s = np.array(list(edeps_s), dtype=np.float64)
        else:
            arr_s = arr_b
        e_tot = float(arr_b.sum())
        e_sig = float(arr_s.sum())
        e_bkg = e_tot - e_sig
        nb = int(nbkg_vec[0]) if has_nbkg else 0
        cls = clusters.get(evt, [])
        n_cl = len(cls)

        # ---- historical energy metrics (denominator e_tot > 0) ------------
        if e_tot > 0.0:
            e_max = e_in = 0.0
            claimed = np.zeros(arr_b.size, dtype=bool)
            for mods in cls:
                e_cl = float(arr_b[mods].sum())
                e_in += e_cl
                if e_cl > e_max:
                    e_max = e_cl
                claimed[mods] = True
            lit_mask = arr_b >= eth_wave
            e_lit = e_lit_out = 0.0
            if lit_mask.any():
                e_lit = float(arr_b[lit_mask].sum())
                e_lit_out = float(arr_b[lit_mask & ~claimed].sum())
            total += e_tot
            leading += e_max
            captured += e_in
            lit_total += e_lit
            lit_unclustered += e_lit_out
            containments.append(e_max / e_tot)
            if e_in >= 0.999 * e_tot:
                n_all_one += 1

        # ---- background-aware per-event quantities ------------------------
        e_lead_sig = e_lead_tot = e_lead_bkg = 0.0
        sig_in = bkg_in = 0.0
        fake = 0
        for mods in cls:
            e_s = float(arr_s[mods].sum())
            e_t = float(arr_b[mods].sum())
            e_b = e_t - e_s
            sig_in += e_s
            bkg_in += e_b
            if e_s < FAKE_SIGNAL_THRESHOLD:
                fake += 1
            if e_t > e_lead_tot:
                e_lead_tot, e_lead_sig, e_lead_bkg = e_t, e_s, e_b

        esig_all += e_sig
        ebkg_all += e_bkg
        sig_in_all += sig_in
        bkg_in_all += bkg_in
        fake_all += fake
        if fake > 0:
            n_fake_events += 1
        nbkg_list.append(nb)
        esig_list.append(e_sig)
        ebkg_list.append(e_bkg)

        if e_sig > 0 and n_cl > 0:
            lead_sig_cont.append(e_lead_sig / e_sig)
            bias_sig.append((e_lead_sig - e_sig) / e_sig)
            if e_lead_tot > 0:
                lead_bkg_frac.append(e_lead_bkg / e_lead_tot)
                bias_tot.append((e_lead_tot - e_sig) / e_sig)
        if n_cl > 0:
            den = sig_in + bkg_in
            if den > 0:
                purity.append(sig_in / den)

        # ---- binned curves ------------------------------------------------
        ei = int(np.searchsorted(E_SIG_EDGES, e_sig, side="right") - 1)
        ei = min(max(ei, 0), len(E_SIG_EDGES) - 2)
        bi = min(nb, NBKG_VALUES[-1])
        for st, b in ((bins_e, ei), (bins_n, bi)):
            st["n"][b] += 1
            if n_cl > 0:
                st["n_with"][b] += 1
            if n_cl == 1:
                st["n_exact"][b] += 1
            if n_cl >= 2:
                st["n_over"][b] += 1
            st["sum_clusters"][b] += n_cl
            st["sum_fake"][b] += fake
            st["sum_esig"][b] += e_sig
            st["sum_ebkg"][b] += e_bkg
            st["sum_lead_sig"][b] += e_lead_sig
            st["sum_lead_tot"][b] += e_lead_tot
            st["sum_lead_bkg"][b] += e_lead_bkg
            st["sum_sig_in"][b] += sig_in
            st["sum_bkg_in"][b] += bkg_in
            if e_sig > 0 and n_cl > 0:
                st["sum_bias_sig"][b] += (e_lead_sig - e_sig) / e_sig
                st["sum_bias_sig2"][b] += ((e_lead_sig - e_sig) / e_sig) ** 2
                st["n_bias_sig"][b] += 1
                if e_lead_tot > 0:
                    st["sum_bias_tot"][b] += (e_lead_tot - e_sig) / e_sig
                    st["sum_bias_tot2"][b] += ((e_lead_tot - e_sig) / e_sig) ** 2
                    st["n_bias_tot"][b] += 1

    if fin_s is not fin_b:
        fin_s.Close()
    fin_b.Close()

    scalar = {
        "has_bkg": int(is_bkg),
        "n_bkg_mean": _mean(nbkg_list),
        "frac_events_nbkg_ge1": (float(np.mean(np.asarray(nbkg_list) >= 1))
                                 if nbkg_list else 0.0),
        "signal_energy_mean": _mean(esig_list),
        "bkg_energy_mean": _mean(ebkg_list),
        "signal_capture_fraction": (sig_in_all / esig_all
                                    if esig_all > 0 else float("nan")),
        "bkg_capture_fraction": (bkg_in_all / ebkg_all
                                 if ebkg_all > 0 else float("nan")),
        "lead_signal_containment_mean": _mean(lead_sig_cont),
        "lead_bkg_fraction_mean": _mean(lead_bkg_frac),
        "signal_purity_mean": _mean(purity),
        "fake_clusters_per_event": fake_all / n_events,
        "frac_events_with_fake_cluster": n_fake_events / n_events,
        "energy_bias_total": _mean(bias_tot),
        "energy_rms_total": _rms(bias_tot),
        "energy_bias_signal": _mean(bias_sig),
        "energy_rms_signal": _rms(bias_sig),
    }

    # response summary restricted to the well-measured signal region
    # (E_sig >= 20 MeV), where the ratio is not dominated by floor events
    hi = E_SIG_EDGES[:-1] >= 20.0

    def _agg(stats, sum_key, n_key):
        num = float(stats[sum_key][hi].sum())
        den = int(stats[n_key][hi].sum())
        return num / den if den else float("nan")

    bias_t_hi = _agg(bins_e, "sum_bias_tot", "n_bias_tot")
    var_t_hi = _agg(bins_e, "sum_bias_tot2", "n_bias_tot") - bias_t_hi ** 2
    bias_s_hi = _agg(bins_e, "sum_bias_sig", "n_bias_sig")
    var_s_hi = _agg(bins_e, "sum_bias_sig2", "n_bias_sig") - bias_s_hi ** 2
    scalar.update({
        "energy_bias_total_hi": bias_t_hi,
        "energy_rms_total_hi": float(np.sqrt(max(var_t_hi, 0.0))),
        "energy_bias_signal_hi": bias_s_hi,
        "energy_rms_signal_hi": float(np.sqrt(max(var_s_hi, 0.0))),
    })

    if containments:
        c = np.array(containments)
        scalar.update({
            "energy_containment_mean": float(c.mean()),
            "energy_containment_median": float(np.median(c)),
            "energy_captured_fraction": float(captured / total),
            "energy_split_fraction": float((captured - leading) / total),
            "energy_unclustered_fraction": float(1.0 - captured / total),
            "lit_energy_total_mev": float(lit_total),
            "lit_energy_missed_fraction": (float(lit_unclustered / lit_total)
                                           if lit_total > 0 else float("nan")),
            "frac_events_containment_ge_90": float(np.mean(c >= 0.90)),
            # historical name; implemented as "all energy captured by any
            # cluster(s) >= 99.9%", denominator e_tot>0
            "frac_events_all_energy_one_cluster": float(
                n_all_one / len(containments)),
        })

    bins = {
        "e_sig_edges": E_SIG_EDGES.tolist(),
        "e_sig": _finalize_bins(bins_e),
        "nbkg_values": NBKG_VALUES,
        "nbkg": _finalize_bins(bins_n),
    }
    return scalar, bins


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
        "cluster_efficiency": float(np.mean(counts >= 1)),
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
    scalar, bins = energy_metrics(dataset, n, out, eth_wave)
    row.update(scalar)
    row["has_bkg"] = int(has_bkg(dataset))

    bins_dir = TASK_DIR / "bins"
    bins_dir.mkdir(parents=True, exist_ok=True)
    with open(bins_dir / f"{tag}.json", "w") as f:
        json.dump(bins, f, indent=1, sort_keys=True)

    with open(TASK_DIR / f"{tag}.json", "w") as f:
        json.dump(row, f, indent=1, sort_keys=True)
    extra = ""
    if row["has_bkg"]:
        extra = (f" fake/ev={row['fake_clusters_per_event']:.2f}"
                 f" purity={row['signal_purity_mean']:.3f}")
    print(f"{tag:34s} mean={row['mean_clusters']:.3f} "
          f"exact={row['exact_match_rate']:.3f} "
          f"sigcap={row.get('signal_capture_fraction', float('nan')):.3f}"
          f"{extra} ({recon_seconds:.1f}s)")
    return row


if __name__ == "__main__":
    run_point(sys.argv[1], float(sys.argv[2]), float(sys.argv[3]),
              int(float(sys.argv[4])) if len(sys.argv) > 4 else N_NEIGHBOR_THRESHOLD)
