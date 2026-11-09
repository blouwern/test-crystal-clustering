"""Stage 2 evaluation: reconstruct the stitched samples and measure multi-particle metrics.

Grid: seeds {7, 8, 10, 13, 14} x waves {0.02, 0.05} x Nth {0, 1}, for both
scenarios (2gamma, ep2gamma) and both versions (signal-only control, +background).
80 reconstruction points total.

The signal-only and background versions share the same signal truth; the
background version additionally has per-event background arrays from the
``BkgCrystal`` / ``BkgEnergy`` branches.
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PYTRIAL = HERE.parents[2]
REPO = HERE.parents[3]
sys.path.insert(0, str(PYTRIAL))
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(REPO / ".mplconfig"))

import reconstruction  # noqa: E402
from evaluation import background, counts, energy, io, matching  # noqa: E402

DATA = HERE / "data"
RESULTS = HERE / "results"
WORK = HERE / "work"
POINTS = RESULTS / "points"

SCENARIOS = {
    "2gamma": {"ntrack": 2, "bkg": DATA / "stitched_2gamma_bkg.root"},
    "ep2gamma": {"ntrack": 3, "bkg": DATA / "stitched_ep2gamma_bkg.root"},
}
VERSIONS = ("sig", "bkg")
SEEDS = (5, 6, 7, 8, 10, 13, 14)
WAVES = (0.02, 0.05)
NTHS = (0, 1)
N_EVENTS = 20_000
FAKE_THRESHOLD = 1.0
MIN_PRIMARY_ENERGY = 20.0


def prepare_workdir():
    proc = WORK / "processed"
    util = WORK / "utilities"
    proc.mkdir(parents=True, exist_ok=True)
    util.mkdir(parents=True, exist_ok=True)
    for scenario in SCENARIOS:
        for version in VERSIONS:
            src = DATA / f"stitched_{scenario}_{version}.root"
            dst = proc / src.name
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src.resolve())
    topo = REPO / "data" / "utilities" / "ecal_neighbor_info.root"
    dst = util / topo.name
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(topo.resolve())


def evaluate_point(truth, clusters, version, n_truth):
    """Aggregate all Stage 2 metrics for one reconstructed point."""
    n_events = truth["n_events"]
    zero_bkg = np.zeros(truth["signal_arrays"][0].shape[1], dtype=np.float32)

    counts_total = np.zeros(n_events, dtype=np.int64)
    counts_signal = np.zeros(n_events, dtype=np.int64)

    eff_sum = 0.0
    eff_n = 0
    mp_sum = 0.0
    mp_n = 0
    matched_ok = 0
    merge_prim = 0
    merge_events = 0
    split_prim = 0
    fake_total = 0
    fake_events = 0

    sig_in_all = bkg_in_all = 0.0
    e_sig_all = e_bkg_all = 0.0
    purity_sum = purity_n = 0.0
    lead_cont_sum = lead_cont_n = 0.0
    lead_bkg_sum = lead_bkg_n = 0.0
    energy_biases = []

    for i in range(n_events):
        sig = truth["signal_arrays"][i]
        bkg = truth["bkg_arrays"][i] if version == "bkg" else zero_bkg
        e_sig = float(sig.sum())
        e_bkg = float(bkg.sum()) if version == "bkg" else 0.0
        cluster_mods = clusters[i]

        if cluster_mods:
            c_sig = np.array([float(sig[:, mods].sum()) for mods in cluster_mods])
            c_bkg = np.array([float(bkg[mods].sum()) for mods in cluster_mods])
        else:
            c_sig = np.zeros(0)
            c_bkg = np.zeros(0)

        n_total = len(cluster_mods)
        n_signal = int(np.sum(c_sig >= FAKE_THRESHOLD))
        counts_total[i] = n_total
        counts_signal[i] = n_signal

        overlap = matching.build_overlap_matrix(cluster_mods, sig)
        assignment = matching.assign(overlap)
        eff, match_purity = matching.efficiency_purity(
            overlap, assignment, truth["source_energy"][i])
        ms = matching.merge_split(
            overlap, truth["source_energy"][i])
        bg = background.event_background(
            c_sig, c_bkg, e_sig, e_bkg, FAKE_THRESHOLD)
        en = energy.event_energy(
            c_sig, c_bkg, e_sig, overlap, assignment,
            truth["source_energy"][i], MIN_PRIMARY_ENERGY)

        eff_sum += float(eff.sum())
        eff_n += int(eff.size)

        signal_like = c_sig >= FAKE_THRESHOLD
        mp_sum += float(match_purity[signal_like].sum()) if signal_like.any() else 0.0
        mp_n += int(signal_like.sum())

        for k, j in enumerate(assignment):
            if j >= 0 and eff[j] > 0.5 and match_purity[k] > 0.5:
                matched_ok += 1

        merge_prim += ms["merged_primaries"]
        merge_events += int(ms["merge_event"])
        split_prim += ms["split_primaries"]
        fake_total += bg["n_fake"]
        fake_events += int(bg["has_fake"])

        sig_in = float(c_sig.sum())
        bkg_in = float(c_bkg.sum())
        sig_in_all += sig_in
        bkg_in_all += bkg_in
        e_sig_all += e_sig
        e_bkg_all += e_bkg

        if n_total > 0:
            purity_sum += bg["signal_purity"]
            purity_n += 1
            if e_sig > 0:
                lead_cont_sum += en["lead_signal_containment"]
                lead_cont_n += 1
            lead_bkg_sum += bg["lead_bkg_fraction"]
            lead_bkg_n += 1
        energy_biases.extend(en["energy_biases"])

    cm_total = counts.count_metrics(counts_total, n_truth)
    cm_signal = counts.count_metrics(counts_signal, n_truth)
    total_primaries = float(eff_n)  # N per event

    bias_mean = float(np.mean(energy_biases)) if energy_biases else float("nan")
    bias_rms = (float(np.sqrt(np.mean(np.square(energy_biases))))
                if energy_biases else float("nan"))

    row = {
        "n_events": int(n_events),
        "ntrack": int(n_truth),
        "exact_total": cm_total["exact_rate"],
        "under_total": cm_total["under_rate"],
        "over_total": cm_total["over_rate"],
        "mae_total": cm_total["MAE"],
        "mean_clusters_total": cm_total["mean_clusters"],
        "exact_signal": cm_signal["exact_rate"],
        "under_signal": cm_signal["under_rate"],
        "over_signal": cm_signal["over_rate"],
        "mae_signal": cm_signal["MAE"],
        "mean_clusters_signal": cm_signal["mean_clusters"],
        "fake_clusters_per_event": fake_total / n_events,
        "frac_events_with_fake": fake_events / n_events,
        "signal_purity_energy": (sig_in_all / (sig_in_all + bkg_in_all)
                                 if sig_in_all + bkg_in_all > 0 else float("nan")),
        "signal_purity_mean": purity_sum / purity_n if purity_n else float("nan"),
        "bkg_capture": (bkg_in_all / e_bkg_all
                        if e_bkg_all > 0 else float("nan")),
        "signal_capture": sig_in_all / e_sig_all if e_sig_all > 0 else float("nan"),
        "lead_signal_containment": (lead_cont_sum / lead_cont_n
                                    if lead_cont_n else float("nan")),
        "lead_bkg_fraction": (lead_bkg_sum / lead_bkg_n
                              if lead_bkg_n else float("nan")),
        "efficiency_mean": eff_sum / total_primaries if total_primaries else float("nan"),
        "match_purity_mean": mp_sum / mp_n if mp_n else float("nan"),
        "matched_primary_fraction": (matched_ok / total_primaries
                                     if total_primaries else float("nan")),
        "merge_rate_primaries": merge_prim / total_primaries if total_primaries else float("nan"),
        "merge_rate_events": merge_events / n_events,
        "split_rate_primaries": split_prim / total_primaries if total_primaries else float("nan"),
        "energy_bias_mean": bias_mean,
        "energy_rms": bias_rms,
        "n_energy": len(energy_biases),
    }
    detail = {
        "row": row,
        "hist_total": cm_total["histogram"],
        "hist_signal": cm_signal["histogram"],
    }
    return row, detail


def pick_best(rows):
    return max(rows, key=lambda r: (r["exact_signal"],
                                    r["signal_purity_energy"],
                                    -r["fake_clusters_per_event"]))


def plot_scenario(rows, scenario):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {("sig", 0): "#1f77b4", ("sig", 1): "#17becf",
              ("bkg", 0): "#d62728", ("bkg", 1): "#ff9896"}
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    for version in VERSIONS:
        for nth in NTHS:
            sel = [r for r in rows if r["scenario"] == scenario
                   and r["version"] == version and r["wave"] == 0.02
                   and r["nth"] == nth]
            sel.sort(key=lambda r: r["seed"])
            if not sel:
                continue
            ax.plot([r["seed"] for r in sel], [r["exact_signal"] for r in sel],
                    "o-", color=colors[(version, nth)],
                    label=f"{version}, Nth={nth}")
    ax.set_xlabel("Eth_seed (MeV)")
    ax.set_ylabel("exact-N (signal clusters)")
    ax.set_title(f"{scenario}: exact-N vs seed (wave=0.02)")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = RESULTS / f"exact_signal_vs_seed_{scenario}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def plot_fake_purity(rows, scenario):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for nth in NTHS:
        sel = [r for r in rows if r["scenario"] == scenario
               and r["version"] == "bkg" and r["wave"] == 0.02
               and r["nth"] == nth]
        sel.sort(key=lambda r: r["seed"])
        if not sel:
            continue
        xs = [r["seed"] for r in sel]
        axes[0].plot(xs, [r["fake_clusters_per_event"] for r in sel],
                     "o-", label=f"Nth={nth}")
        axes[1].plot(xs, [r["signal_purity_energy"] for r in sel],
                     "o-", label=f"Nth={nth}")
    axes[0].set_xlabel("Eth_seed (MeV)")
    axes[0].set_ylabel("fake clusters / event")
    axes[0].set_title("background-only clusters")
    axes[0].grid(alpha=0.3)
    axes[0].legend()
    axes[1].set_xlabel("Eth_seed (MeV)")
    axes[1].set_ylabel("signal purity (energy-weighted)")
    axes[1].set_title("signal purity")
    axes[1].grid(alpha=0.3)
    axes[1].legend()
    fig.suptitle(f"{scenario}: background-induced failures (wave=0.02)")
    fig.tight_layout()
    out = RESULTS / f"fake_purity_vs_seed_{scenario}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def plot_best_bars(best):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [f"{b['scenario']}\n{b['version']}" for b in best]
    x = np.arange(len(best))
    width = 0.38
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ax.bar(x - width / 2, [b["merge_rate_events"] for b in best], width,
           label="merge rate (events)")
    ax.bar(x + width / 2, [b["split_rate_primaries"] for b in best], width,
           label="split rate (primaries)")
    ax.set_xticks(x, labels, fontsize=8)
    ax.set_ylabel("rate")
    ax.set_title("Best Stage 2 points: merge and split")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = RESULTS / "merge_split_best.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def plot_count_hist(best):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, item in zip(axes, best):
        detail = json.loads((POINTS / f"{item['tag']}.json").read_text())
        hist_total = {int(k): int(v) for k, v in detail["hist_total"].items()}
        hist_signal = {int(k): int(v) for k, v in detail["hist_signal"].items()}
        max_k = max(max(hist_total, default=0), max(hist_signal, default=0),
                    item["ntrack"])
        for hist, offset, label in ((hist_total, -0.2, "all clusters"),
                                    (hist_signal, 0.2, "signal-like clusters")):
            ks = sorted(hist)
            ax.bar([k + offset for k in ks], [hist[k] for k in ks],
                   width=0.38, label=label)
        ax.axvline(item["ntrack"], color="red", linestyle="--", label="truth N")
        ax.set_xlabel("clusters per event")
        ax.set_ylabel("events")
        ax.set_title(f"{item['scenario']} ({item['version']}: "
                     f"seed={item['seed']:g}, w={item['wave']:g}, Nth={item['nth']})",
                     fontsize=9)
        ax.set_xticks(range(0, max_k + 2))
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    out = RESULTS / "count_hist_best.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    POINTS.mkdir(parents=True, exist_ok=True)
    prepare_workdir()

    truth = {sc: io.load_multiparticle_sample(cfg["bkg"])
             for sc, cfg in SCENARIOS.items()}
    for sc, t in truth.items():
        print(f"{sc}: {t['n_events']} events, N={SCENARIOS[sc]['ntrack']}, "
              f"mean E_sig={t['signal_energy'].mean():.2f} MeV")

    rows = []
    total_points = len(SCENARIOS) * len(VERSIONS) * len(SEEDS) * len(WAVES) * len(NTHS)
    done = 0
    for scenario, cfg in SCENARIOS.items():
        n_truth = cfg["ntrack"]
        for version in VERSIONS:
            sample = DATA / f"stitched_{scenario}_{version}.root"
            for seed in SEEDS:
                for wave in WAVES:
                    for nth in NTHS:
                        tag = (f"{scenario}_{version}_seed{seed:g}_"
                               f"wave{wave:g}_nth{nth}")
                        cluster_path = RESULTS / f"clusters_{tag}.root"
                        reconstruction.reconstruction(
                            filename_edep_of_each_evt=sample.name,
                            energy_threshold_seed=float(seed),
                            energy_threshold_wave=float(wave),
                            n_neighbor_threshold=int(nth),
                            data_dir=str(WORK),
                            output_path=str(cluster_path),
                        )
                        clusters = io.load_clusters(cluster_path, N_EVENTS)
                        row, detail = evaluate_point(
                            truth[scenario], clusters, version, n_truth)
                        row.update({
                            "scenario": scenario,
                            "version": version,
                            "seed": float(seed),
                            "wave": float(wave),
                            "nth": int(nth),
                            "tag": tag,
                        })
                        (POINTS / f"{tag}.json").write_text(
                            json.dumps(detail, indent=2) + "\n")
                        rows.append(row)
                        done += 1
                        print(f"[{done:2d}/{total_points}] {tag}: "
                              f"exact_signal={row['exact_signal']:.4f}, "
                              f"fake={row['fake_clusters_per_event']:.4f}, "
                              f"merge_ev={row['merge_rate_events']:.4f}",
                              flush=True)

    fieldnames = list(rows[0].keys())
    with open(RESULTS / "stage2_metrics.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    best = []
    for scenario in SCENARIOS:
        for version in VERSIONS:
            sub = [r for r in rows if r["scenario"] == scenario
                   and r["version"] == version]
            pick = pick_best(sub)
            best.append(pick)
    best_payload = []
    for b in best:
        payload = {k: b[k] for k in fieldnames}
        best_payload.append(payload)
    (RESULTS / "best_points.json").write_text(
        json.dumps(best_payload, indent=2) + "\n")

    # background impact at the best background point (paired signal-only row)
    impact = []
    for scenario in SCENARIOS:
        b = pick_best([r for r in rows if r["scenario"] == scenario
                       and r["version"] == "bkg"])
        s = [r for r in rows if r["scenario"] == scenario
             and r["version"] == "sig" and r["seed"] == b["seed"]
             and r["wave"] == b["wave"] and r["nth"] == b["nth"]][0]
        impact.append({
            "scenario": scenario,
            "seed": b["seed"], "wave": b["wave"], "nth": b["nth"],
            "exact_signal_bkg": b["exact_signal"],
            "exact_signal_sig": s["exact_signal"],
            "delta_exact_signal": b["exact_signal"] - s["exact_signal"],
            "fake_clusters_per_event": b["fake_clusters_per_event"],
            "signal_purity_energy": b["signal_purity_energy"],
            "delta_merge_rate_events": b["merge_rate_events"] - s["merge_rate_events"],
            "delta_split_rate_primaries": (b["split_rate_primaries"]
                                           - s["split_rate_primaries"]),
            "delta_efficiency_mean": b["efficiency_mean"] - s["efficiency_mean"],
        })
    (RESULTS / "background_impact.json").write_text(
        json.dumps(impact, indent=2) + "\n")

    plot_best_bars([b for b in best if b["version"] == "bkg"])
    plot_count_hist([b for b in best if b["version"] == "bkg"])
    for scenario in SCENARIOS:
        plot_scenario(rows, scenario)
        plot_fake_purity(rows, scenario)

    print(f"wrote {RESULTS / 'stage2_metrics.csv'} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
