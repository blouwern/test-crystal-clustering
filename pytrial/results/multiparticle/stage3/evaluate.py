"""Stage 3 evaluation: merge rate versus two-shower separation.

The sample contains 2 500 events in each separation bin (0-0.5, 0.5-1, 1-2,
2-3, 3-5, >5 pitch).  For every grid point we reconstruct the signal-only and
+background versions, match the two truth primaries to clusters, and measure:

* merge rate: both primaries have the same best cluster;
* split rate: a primary's best cluster holds <90% of its signal energy;
* exact-2 (signal clusters);
* efficiency and match purity;
* background fake clusters and signal purity (for the +background version).

The 10% and 50% merge radii are obtained by linear interpolation between the
bin medians.
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
from evaluation import background, counts, io, matching, separation  # noqa: E402

DATA = HERE / "data"
RESULTS = HERE / "results"
WORK = HERE / "work"
POINTS = RESULTS / "points"
STAGE2_CSV = HERE.parent / "stage2/results/stage2_metrics.csv"

VERSIONS = ("sig", "bkg")
SEEDS = (5, 6, 7, 8, 10, 13, 14)
NTHS = (0, 1)
WAVE = 0.05
N_EVENTS = 15_000
FAKE_THRESHOLD = 1.0
BINS = separation.DEFAULT_EDGES_PITCH
N_BINS = len(BINS) - 1


def prepare_workdir():
    proc = WORK / "processed"
    util = WORK / "utilities"
    proc.mkdir(parents=True, exist_ok=True)
    util.mkdir(parents=True, exist_ok=True)
    for version in VERSIONS:
        src = DATA / f"stage3_2gamma_{version}.root"
        dst = proc / src.name
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(src.resolve())
    topo = REPO / "data" / "utilities" / "ecal_neighbor_info.root"
    dst = util / topo.name
    if dst.is_symlink() or dst.exists():
        dst.unlink()
    dst.symlink_to(topo.resolve())


def evaluate_point(truth, clusters, version, pitch):
    """Return one row per separation bin for one reconstructed point."""
    n_events = truth["n_events"]
    zero_bkg = np.zeros(truth["signal_arrays"][0].shape[1], dtype=np.float32)
    stats = {
        b: {"n": 0, "merge": 0, "split": 0, "exact2": 0,
            "eff_sum": 0.0, "eff_n": 0, "mp_sum": 0.0, "mp_n": 0,
            "fake": 0, "sig_in": 0.0, "bkg_in": 0.0, "bkg_total": 0.0}
        for b in range(N_BINS)
    }

    for i in range(n_events):
        b = int(truth["sep_bin"][i])
        if b < 0 or b >= N_BINS:
            continue
        st = stats[b]
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

        overlap = matching.build_overlap_matrix(cluster_mods, sig)
        assignment = matching.assign(overlap)
        eff, match_purity = matching.efficiency_purity(
            overlap, assignment, truth["source_energy"][i])
        ms = matching.merge_split(overlap, truth["source_energy"][i])
        bg = background.event_background(
            c_sig, c_bkg, e_sig, e_bkg, FAKE_THRESHOLD)

        n_signal = int(np.sum(c_sig >= FAKE_THRESHOLD))
        st["n"] += 1
        st["merge"] += int(ms["merge_event"])
        st["split"] += int(ms["split_primaries"])
        st["exact2"] += int(n_signal == 2)
        st["eff_sum"] += float(eff.sum())
        st["eff_n"] += int(eff.size)
        signal_like = c_sig >= FAKE_THRESHOLD
        st["mp_sum"] += float(match_purity[signal_like].sum()) if signal_like.any() else 0.0
        st["mp_n"] += int(signal_like.sum())
        st["fake"] += int(bg["n_fake"])
        st["sig_in"] += float(c_sig.sum())
        st["bkg_in"] += float(c_bkg.sum())
        st["bkg_total"] += e_bkg

    rows = []
    for b in range(N_BINS):
        st = stats[b]
        n = st["n"]
        if n == 0:
            continue
        rows.append({
            "bin": b,
            "sep_lo_pitch": float(BINS[b]),
            "sep_hi_pitch": float(BINS[b + 1]),
            "median_sep_mm": float(np.median(
                truth["separation"][truth["sep_bin"] == b])),
            "median_sep_pitch": float(np.median(
                truth["separation"][truth["sep_bin"] == b]) / pitch),
            "n_events": n,
            "merge_rate": st["merge"] / n,
            "split_rate": st["split"] / max(st["eff_n"], 1),
            "exact2": st["exact2"] / n,
            "efficiency_mean": st["eff_sum"] / max(st["eff_n"], 1),
            "match_purity_mean": (st["mp_sum"] / st["mp_n"]
                                  if st["mp_n"] else float("nan")),
            "fake_clusters_per_event": st["fake"] / n,
            "signal_purity": (st["sig_in"] / (st["sig_in"] + st["bkg_in"])
                              if st["sig_in"] + st["bkg_in"] > 0 else float("nan")),
            "bkg_capture": (st["bkg_in"] / st["bkg_total"]
                            if st["bkg_total"] > 0 else float("nan")),
        })
    return rows


def pick_stage2_best():
    """Best corrected Stage 2 2gamma +background point, if available."""
    if not STAGE2_CSV.exists():
        return 14, 0
    rows = list(csv.DictReader(open(STAGE2_CSV)))
    sub = [r for r in rows if r["scenario"] == "2gamma" and r["version"] == "bkg"]
    best = max(sub, key=lambda r: float(r["exact_signal"]))
    return int(float(best["seed"])), int(float(best["nth"]))


def plot_merge_curves(rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_seeds = (5, 7, 10, 14)
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(plot_seeds)))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, version in zip(axes, VERSIONS):
        for color, seed in zip(colors, plot_seeds):
            for nth, ls in ((0, "-"), (1, "--")):
                sel = [r for r in rows if r["version"] == version
                       and r["seed"] == seed and r["nth"] == nth]
                sel.sort(key=lambda r: r["bin"])
                if not sel:
                    continue
                ax.plot([r["median_sep_pitch"] for r in sel],
                        [r["merge_rate"] for r in sel],
                        marker="o", ms=4, ls=ls, color=color,
                        label=f"seed={seed}, Nth={nth}")
        ax.set_xlabel("median separation (pitch)")
        ax.set_ylabel("merge rate (events)")
        ax.set_title("signal only" if version == "sig" else "+ background")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=6, ncol=2)
    fig.suptitle("Stage 3: merge rate vs two-shower separation (wave=0.05)")
    fig.tight_layout()
    out = RESULTS / "merge_rate_vs_separation.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def plot_best(rows, fits, seed, nth):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for version, color in (("sig", "#1f77b4"), ("bkg", "#d62728")):
        sel = [r for r in rows if r["version"] == version
               and r["seed"] == seed and r["nth"] == nth]
        sel.sort(key=lambda r: r["bin"])
        x = np.array([r["median_sep_pitch"] for r in sel])
        y = np.array([r["merge_rate"] for r in sel])
        n = np.array([r["n_events"] for r in sel])
        err = np.sqrt(y * (1 - y) / n)
        ax.errorbar(x, y, yerr=err, marker="o", color=color,
                    label="signal only" if version == "sig" else "+ background")
        fit = [f for f in fits if f["version"] == version
               and f["seed"] == seed and f["nth"] == nth][0]
        if np.isfinite(fit["r50_pitch"]):
            ax.axvline(fit["r50_pitch"], color=color, ls=":",
                       label=f"50% radius = {fit['r50_pitch']:.2f} pitch")
    ax.set_xlabel("median separation (pitch)")
    ax.set_ylabel("merge rate (events)")
    ax.set_ylim(0, 1)
    ax.set_title(f"Stage 3 best Stage 2 point: seed={seed}, Nth={nth}, wave={WAVE:g}")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = RESULTS / "merge_rate_best.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    POINTS.mkdir(parents=True, exist_ok=True)
    prepare_workdir()

    truth = io.load_multiparticle_sample(DATA / "stage3_2gamma_bkg.root")
    pitch = float(json.loads((HERE / "meta.json").read_text())["pitch_mm"])
    print(f"{truth['n_events']} events, pitch={pitch:.2f} mm, "
          f"bins={np.bincount(truth['sep_bin'], minlength=N_BINS).tolist()}")

    rows = []
    total = len(VERSIONS) * len(SEEDS) * len(NTHS)
    done = 0
    for version in VERSIONS:
        sample = DATA / f"stage3_2gamma_{version}.root"
        for seed in SEEDS:
            for nth in NTHS:
                tag = f"stage3_{version}_seed{seed:g}_nth{nth}"
                cluster_path = RESULTS / f"clusters_{tag}.root"
                reconstruction.reconstruction(
                    filename_edep_of_each_evt=sample.name,
                    energy_threshold_seed=float(seed),
                    energy_threshold_wave=float(WAVE),
                    n_neighbor_threshold=int(nth),
                    data_dir=str(WORK),
                    output_path=str(cluster_path),
                )
                clusters = io.load_clusters(cluster_path, N_EVENTS)
                point_rows = evaluate_point(truth, clusters, version, pitch)
                for r in point_rows:
                    r.update({"version": version, "seed": float(seed),
                              "nth": int(nth), "wave": float(WAVE)})
                rows.extend(point_rows)
                done += 1
                best_bin = max(point_rows, key=lambda r: r["merge_rate"])
                print(f"[{done:2d}/{total}] {tag}: "
                      f"merge(bin0)={point_rows[0]['merge_rate']:.3f}, "
                      f"merge(max)={best_bin['merge_rate']:.3f}", flush=True)

    fieldnames = list(rows[0].keys())
    with open(RESULTS / "stage3_metrics.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    fits = []
    for version in VERSIONS:
        for seed in SEEDS:
            for nth in NTHS:
                sel = [r for r in rows if r["version"] == version
                       and r["seed"] == seed and r["nth"] == nth]
                sel.sort(key=lambda r: r["bin"])
                x = [r["median_sep_pitch"] for r in sel]
                y = [r["merge_rate"] for r in sel]
                fits.append({
                    "version": version, "seed": seed, "nth": nth,
                    "r50_pitch": separation.crossing_radius(x, y, 0.5),
                    "r10_pitch": separation.crossing_radius(x, y, 0.1),
                })
    for f in fits:
        for key in ("r50_pitch", "r10_pitch"):
            f[key + "_mm"] = f[key] * pitch if np.isfinite(f[key]) else float("nan")
    (RESULTS / "fits.json").write_text(json.dumps(fits, indent=2) + "\n")

    seed, nth = pick_stage2_best()
    plot_merge_curves(rows)
    plot_best(rows, fits, seed, nth)
    print(f"wrote {RESULTS / 'stage3_metrics.csv'} ({len(rows)} rows); "
          f"best Stage 2 point seed={seed}, Nth={nth}")


if __name__ == "__main__":
    main()
