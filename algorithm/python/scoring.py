import ROOT
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def _read_cluster_counts(input_path):
    with ROOT.TFile(str(input_path), "READ") as f:
        tree = f.Get("ECALClusters")
        n_evt = int(tree.GetMaximum("EvtID") + 1)
        n_cluster_list = [0] * n_evt
        for entry in tree:
            n_cluster_list[entry.EvtID] += 1
    return np.array(n_cluster_list)


def _read_ntrack(data_path, n_expected):
    with ROOT.TFile(str(data_path), "READ") as f:
        data_tree = f.Get("EdepOfEachEvt")
        if data_tree.GetBranch("NTrack") is None:
            return None
        n_track_list = []
        for entry in data_tree:
            n_track_list.append(entry.NTrack)
            if len(n_track_list) >= n_expected:
                break
    return np.array(n_track_list)


def _run_chisquare(observed_arr, expected_arr, ddof, label):
    valid_mask = expected_arr >= 5.0
    if np.sum(valid_mask) < 2:
        print(f"\n  Chi-squared test ({label}): insufficient bins with expected >= 5")
        return

    merged_obs = []
    merged_exp = []
    pool_o, pool_e = 0.0, 0.0
    for o, e in zip(observed_arr, expected_arr):
        if e >= 5.0:
            merged_obs.append(o + pool_o)
            merged_exp.append(e + pool_e)
            pool_o, pool_e = 0.0, 0.0
        else:
            pool_o += o
            pool_e += e
    if pool_e > 0 and merged_obs:
        merged_obs[-1] += pool_o
        merged_exp[-1] += pool_e

    chi2, p_value = stats.chisquare(merged_obs, merged_exp, ddof=ddof)
    print(f"\n  Chi-squared test ({label}):")
    print(f"    chi2 = {chi2:.4f}")
    print(f"    p    = {p_value:.4f}")
    if p_value < 0.05:
        print(f"    => Significant deviation (p < 0.05)")
    else:
        print(f"    => No significant deviation (p >= 0.05)")


def _run_event_by_event(data, n_track_arr, N, input_path, output_plot):
    diff = data - n_track_arr
    abs_diff = np.abs(diff)
    mae = np.mean(abs_diff)
    rmse = np.sqrt(np.mean(diff.astype(np.float64) ** 2))
    exact_match = np.sum(diff == 0)
    over_cluster = np.sum(diff > 0)
    under_cluster = np.sum(diff < 0)

    print("=" * 60)
    print("  Event-by-event scoring")
    print("=" * 60)
    print(f"  Number of events:             {N}")
    print(f"  Mean predicted clusters:      {np.mean(data):.4f}")
    print(f"  Mean ground-truth tracks:     {np.mean(n_track_arr):.4f}")
    print(f"  Exact match rate:             {exact_match}/{N} ({exact_match / N * 100:.2f}%)")
    print(f"  Over-clustered events:        {over_cluster}/{N} ({over_cluster / N * 100:.2f}%)")
    print(f"  Under-clustered events:       {under_cluster}/{N} ({under_cluster / N * 100:.2f}%)")
    print(f"  MAE:                          {mae:.4f}")
    print(f"  RMSE:                         {rmse:.4f}")
    print(f"  Pearson correlation:          {np.corrcoef(data, n_track_arr)[0, 1]:.4f}")
    print()

    max_k = max(np.max(data), np.max(n_track_arr))
    print(f"  {'Diff':>6s}  {'Count':>6s}  {'Fraction':>9s}")
    print(f"  {'-' * 6}  {'-' * 6}  {'-' * 9}")
    for d in sorted(set(diff)):
        c = int(np.sum(diff == d))
        print(f"  {d:6d}  {c:6d}  {c / N:9.4f}")
    print("=" * 60)

    bins = np.arange(-0.5, max_k + 1.5, 1)
    plt.figure(figsize=(8, 6))
    plt.hist(
        data, bins=bins, rwidth=0.4, color="steelblue", edgecolor="black",
        alpha=0.7, label="Predicted clusters",
    )
    plt.hist(
        n_track_arr, bins=bins, rwidth=0.4, color="orange", edgecolor="black",
        alpha=0.7, label="Ground-truth tracks",
    )
    plt.xlabel("Count per event")
    plt.ylabel("Number of events")
    plt.title(f"Predicted vs Ground-truth (N={N})")
    plt.xticks(range(int(max_k) + 1))
    plt.ylim(0, None)
    plt.legend()
    plt.grid(axis="y", linestyle="--", alpha=0.4)
    plt.tight_layout()
    if output_plot is None:
        output_plot = "event_by_event_" + input_path.stem + ".png"
    plt.savefig(output_plot, dpi=150)


def _run_statistical(data, expected_lambda, n_track_arr, N, input_path, output_plot):
    mean_val = np.mean(data)
    std_val = np.std(data, ddof=1)
    median_val = np.median(data)
    min_val = np.min(data)
    max_val = np.max(data)

    unique, freq = np.unique(data, return_counts=True)

    print("=" * 60)
    print("  Cluster-count distribution analysis")
    print("=" * 60)
    print(f"  Number of events:           {N}")
    print(f"  Mean number of clusters:    {mean_val:.4f}")
    print(f"  Standard deviation:         {std_val:.4f}")
    print(f"  Median:                     {median_val}")
    print(f"  Range:                      [{min_val}, {max_val}]")
    print()
    print(f"  {'Clusters':>10s}  {'Count':>6s}  {'Fraction':>9s}")
    print(f"  {'-' * 10}  {'-' * 6}  {'-' * 9}")
    for k, c in zip(unique, freq):
        print(f"  {k:10d}  {c:6d}  {c / N:9.4f}")
    print()

    if n_track_arr is not None:
        _run_statistical_empirical(data, n_track_arr, N, max_val)
        exp_val_int = int(round(np.mean(n_track_arr)))
        ks_stat, ks_p = stats.ks_2samp(data, n_track_arr)
        print(f"\n  Two-sample KS test (clusters vs ground-truth tracks):")
        print(f"    D = {ks_stat:.4f}")
        print(f"    p = {ks_p:.4f}")
        if ks_p < 0.05:
            print(f"    => Distributions differ significantly (p < 0.05)")
        else:
            print(f"    => Distributions are consistent (p >= 0.05)")
    else:
        if expected_lambda is None:
            expected_lambda = mean_val
            ddof = 1
        else:
            ddof = 0
        print(f"  Expected Poisson lambda:    {expected_lambda:.4f}")
        _run_poisson_tests(data, expected_lambda, ddof, N, max_val)
        exp_val_int = int(round(expected_lambda))
        ks_stat, ks_p = stats.kstest(data, "poisson", args=(expected_lambda,))
        print(f"\n  Kolmogorov-Smirnov test against Poisson:")
        print(f"    D = {ks_stat:.4f}")
        print(f"    p = {ks_p:.4f}")
        if ks_p < 0.05:
            print(f"    => Significant deviation from Poisson (p < 0.05)")
        else:
            print(f"    => No significant deviation from Poisson (p >= 0.05)")

    print("=" * 60)

    bins = np.arange(-0.5, max_val + 1.5, 1)
    plt.figure(figsize=(8, 6))
    plt.hist(
        data, bins=bins, rwidth=0.8, color="steelblue",
        edgecolor="black", alpha=0.7,
    )
    plt.axvline(
        x=exp_val_int, color="red", linestyle="--", linewidth=2,
        label=f"Expected value ({exp_val_int})",
    )
    plt.xlabel("Number of clusters per event")
    plt.ylabel("Number of events")
    plt.title(f"ECAL cluster distribution (N={N})")
    plt.xticks(range(int(max_val) + 1))
    plt.ylim(0, None)
    plt.legend()
    plt.grid(axis="y", linestyle="--", alpha=0.4)
    plt.tight_layout()
    if output_plot is None:
        output_plot = "distribution_" + input_path.stem + ".png"
    plt.savefig(output_plot, dpi=150)


def _run_statistical_empirical(data, n_track_arr, N, max_val):
    unique_truth, freq_truth = np.unique(n_track_arr, return_counts=True)
    mean_truth = np.mean(n_track_arr)

    print(f"  Ground-truth mean tracks:   {mean_truth:.4f}")
    print(f"  Ground-truth distribution:")
    print(f"  {'Tracks':>10s}  {'Count':>6s}  {'Fraction':>9s}")
    print(f"  {'-' * 10}  {'-' * 6}  {'-' * 9}")
    for k, c in zip(unique_truth, freq_truth):
        print(f"  {k:10d}  {c:6d}  {c / N:9.4f}")
    print()

    k_values = np.arange(0, int(max_val) + 1)
    observed_arr = np.array([np.sum(data == k) for k in k_values], dtype=np.float64)
    expected_arr = np.array([
        N * (np.sum(n_track_arr == k) / N) if k <= np.max(n_track_arr) else 0.0
        for k in k_values
    ], dtype=np.float64)

    tail_observed = np.sum(data > max_val)
    tail_expected = max(0.0, N - np.sum(expected_arr))
    observed_arr = np.append(observed_arr, tail_observed)
    expected_arr = np.append(expected_arr, tail_expected)

    _run_chisquare(observed_arr, expected_arr, ddof=0, label="vs empirical ground-truth")


def _run_poisson_tests(data, expected_lambda, ddof, N, max_val):
    k_values = np.arange(0, int(max_val) + 1)
    observed_arr = np.array([np.sum(data == k) for k in k_values], dtype=np.float64)
    expected_arr = N * stats.poisson.pmf(k_values, expected_lambda)

    tail_expected = N * (1.0 - stats.poisson.cdf(max_val, expected_lambda))
    tail_observed = np.sum(data > max_val)
    observed_arr = np.append(observed_arr, tail_observed)
    expected_arr = np.append(expected_arr, tail_expected)

    _run_chisquare(
        observed_arr, expected_arr, ddof=ddof,
        label=f"vs Poisson(lambda={expected_lambda:.4f})",
    )


def run(input_path, expected_lambda=None, output_plot=None,
        event_by_event=False, data_path=None):
    data = _read_cluster_counts(input_path)
    N = len(data)

    n_track_arr = None
    if data_path is not None:
        n_track_arr = _read_ntrack(data_path, N)
        if n_track_arr is None:
            print("Warning: data file has no NTrack branch, "
                  "falling back to statistical mode")
            event_by_event = False

    if event_by_event and n_track_arr is not None:
        _run_event_by_event(data, n_track_arr, N, input_path, output_plot)
    else:
        _run_statistical(data, expected_lambda, n_track_arr, N, input_path, output_plot)
