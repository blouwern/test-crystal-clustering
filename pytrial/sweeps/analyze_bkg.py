"""Analyze the gamma-background parameter scan and produce report figures.

Run with the same grid environment used for the scan, e.g.

    ECAL_SWEEP_RESULTS_DATE=26-11-06 \
    ECAL_SWEEP_TAG=bkg \
    ECAL_SWEEP_SEEDS=7,8,9,10,11,12,13,14,15,18 \
    ECAL_SWEEP_WAVES=0.02,0.05,0.1,0.2,0.3,0.5,0.75,1.0,2.0 \
    ECAL_SWEEP_NTHS=0,1,2 \
    ECAL_SWEEP_DATASETS=ep_iso_bkg,gamma_iso_bkg,mu_rest_bkg \
    python3 sweeps/analyze_bkg.py

Outputs go to ``results/sweep/<RESULTS_DATE>/``.  The no-background comparison
uses the 26-10-16-refined CSV by default (``--nobkg-csv``), and any same-code
control rows present in the current CSV are preferred automatically.
"""

import argparse
import array
import csv
import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import analyze as base  # noqa: E402
from sweeps.config import (  # noqa: E402
    DATASETS, ETH_SEED_VALUES, ETH_WAVE_VALUES, NTH_VALUES, RESULTS_DIR,
    SUMMARY_CSV, SUBSET_DIR, TASK_DIR, has_bkg, signal_dataset,
)

BKG_DATASETS = [d for d in DATASETS if has_bkg(d)]
NOBKG_DATASETS = [d for d in DATASETS if not has_bkg(d)]

# figure style (kept close to the 2026-10 deck)
INK = "#20242A"
BODY = "#3A4048"
MUTED = "#6B7280"
GRID = "#E7EAEE"
ACCENT = "#2F5D8C"
WARM = "#B08A62"
GREEN = "#5B8C6E"
RED = "#B5544A"
GREY = "#8A9199"
COLORS = {0: "#2F5D8C", 1: "#5B8C6E", 2: "#B08A62"}
plt.rcParams.update({
    "font.size": 9.5,
    "axes.edgecolor": "#C9CED6",
    "axes.linewidth": 0.8,
    "axes.labelcolor": BODY,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.unicode_minus": False,
    "legend.frameon": False,
})


# --------------------------------------------------------------------------- #
# data loading
# --------------------------------------------------------------------------- #
def load_rows(path=SUMMARY_CSV):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in list(r.items()):
            if k in ("dataset", "neighbor_file"):
                continue
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return rows


def mean_matrix(rows, datasets, key):
    mats = [base.metric_matrix(rows, ds, key) for ds in datasets]
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(np.stack(mats), axis=0)


def diff_matrix(rows_bkg, rows_nobkg, ds, key, nth):
    b = base.metric_matrix(
        [r for r in rows_bkg if int(r["Nth"]) == nth], ds, key)
    # the no-bkg rows use the un-suffixed dataset name
    n = base.metric_matrix(
        [r for r in rows_nobkg if int(r["Nth"]) == nth],
        signal_dataset(ds), key)
    return b - n


# --------------------------------------------------------------------------- #
# heatmaps
# --------------------------------------------------------------------------- #
def heatmap_per_nth(rows, datasets, key, title, filename, *,
                    cmap="RdYlGn", vmin=0.0, vmax=1.0, fmt="{:.3f}",
                    wave_idx=None, annotate=True, cell_w=0.62, cell_h=1.6,
                    val_fs=6.0, label_stride=1):
    nrow, ncol = len(datasets) + 1, len(NTH_VALUES)
    n_wave = len(ETH_WAVE_VALUES) if wave_idx is None else len(list(wave_idx))
    fig, gs = base._heat_grid(nrow, ncol, n_wave, cell_w, cell_h)
    for r, rlab in enumerate(list(datasets) + ["mean over datasets"]):
        for c, nth in enumerate(NTH_VALUES):
            sub = [x for x in rows if int(x["Nth"]) == nth]
            mat = (mean_matrix(sub, datasets, key) if rlab.startswith("mean")
                   else base.metric_matrix(sub, rlab, key))
            ax = fig.add_subplot(gs[r, 3 * c])
            im = base._draw_matrix(
                ax, base.slice_matrix(mat, wave_idx), fmt, cmap, vmin, vmax,
                show_ylabel=(c == 0), wave=wave_idx, annotate=annotate,
                val_fs=val_fs, tick_fs=6.0, label_stride=label_stride)
            if annotate and im is not None:
                # white separators so adjacent cell values stay readable
                for x in np.arange(-0.5, n_wave, 1):
                    ax.axvline(x, color="white", lw=0.6, alpha=0.85, zorder=3)
                for y in np.arange(-0.5, len(ETH_SEED_VALUES), 1):
                    ax.axhline(y, color="white", lw=0.6, alpha=0.85, zorder=3)
            if c == 0:
                ax.set_title(f"{rlab.replace('_bkg', '')}  |  "
                             f"$N_{{th}}$ = {nth}", fontsize=9.5)
            elif r == 0:
                ax.set_title(f"$N_{{th}}$ = {nth}", fontsize=11)
            if r != nrow - 1:
                ax.set_xticks([])
                ax.set_xlabel("")
            cax = fig.add_subplot(gs[r, 3 * c + 1])
            fig.colorbar(im, cax=cax)
    fig.suptitle(title, fontsize=13)
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def comparison_heatmap(rows_bkg, rows_nobkg, key, title, filename, *,
                       cmap="RdBu_r", vmin=-0.05, vmax=0.05,
                       fmt="{:+.3f}", wave_idx=None, annotate=True):
    nths = [n for n in NTH_VALUES
            if any(int(r["Nth"]) == n for r in rows_nobkg)]
    if not nths:
        nths = list(NTH_VALUES)
    nrow, ncol = len(BKG_DATASETS), len(nths)
    n_wave = len(ETH_WAVE_VALUES) if wave_idx is None else len(list(wave_idx))
    fig, gs = base._heat_grid(nrow, ncol, n_wave, 0.72, 1.6)
    for r, ds in enumerate(BKG_DATASETS):
        for c, nth in enumerate(nths):
            mat = diff_matrix(rows_bkg, rows_nobkg, ds, key, nth)
            ax = fig.add_subplot(gs[r, 3 * c])
            im = base._draw_matrix(
                ax, base.slice_matrix(mat, wave_idx), fmt, cmap, vmin, vmax,
                show_ylabel=(c == 0), wave=wave_idx, annotate=annotate,
                val_fs=6.0, tick_fs=6.0)
            if annotate and im is not None:
                # white separators so adjacent cell values stay readable
                for x in np.arange(-0.5, n_wave, 1):
                    ax.axvline(x, color="white", lw=0.6, alpha=0.85, zorder=3)
                for y in np.arange(-0.5, len(ETH_SEED_VALUES), 1):
                    ax.axhline(y, color="white", lw=0.6, alpha=0.85, zorder=3)
            if c == 0:
                ax.set_title(f"{ds.replace('_bkg', '')}  |  "
                             f"$N_{{th}}$ = {nth}", fontsize=9.5)
            elif r == 0:
                ax.set_title(f"$N_{{th}}$ = {nth}", fontsize=11)
            if r != nrow - 1:
                ax.set_xticks([])
                ax.set_xlabel("")
            cax = fig.add_subplot(gs[r, 3 * c + 1])
            fig.colorbar(im, cax=cax)
    fig.suptitle(title, fontsize=13)
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# --------------------------------------------------------------------------- #
# point selection
# --------------------------------------------------------------------------- #
def best_for_nth(rows, nth):
    grid = {}
    for r in rows:
        if int(r["Nth"]) != nth or r["dataset"] not in BKG_DATASETS:
            continue
        grid.setdefault((r["Eth_seed"], r["Eth_wave"]), []).append(
            r["exact_match_rate"])
    if not grid:
        return None
    (s, w), v = max(((k, float(np.mean(vs))) for k, vs in grid.items()),
                    key=lambda kv: kv[1])
    at = [r for r in rows
          if int(r["Nth"]) == nth and r["Eth_seed"] == s
          and r["Eth_wave"] == w and r["dataset"] in BKG_DATASETS]
    return s, w, v, at


def best_for_dataset(rows, ds):
    sub = [r for r in rows if r["dataset"] == ds and "exact_match_rate" in r]
    return max(sub, key=lambda r: r["exact_match_rate"]) if sub else None


def load_bins(ds, seed, wave, nth):
    tag = f"{ds}_seed{seed:g}_wave{wave:g}_nth{nth:g}"
    path = TASK_DIR / "bins" / f"{tag}.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def _centers(edges):
    e = np.asarray(edges, dtype=float)
    return 0.5 * (e[1:] + e[:-1])


def _mask(n, minimum=50):
    return np.asarray(n) >= minimum


# --------------------------------------------------------------------------- #
# curves
# --------------------------------------------------------------------------- #
def energy_efficiency(rows, seed, wave, nth, filename):
    fig, axes = plt.subplots(1, len(BKG_DATASETS),
                             figsize=(4.4 * len(BKG_DATASETS), 3.9),
                             sharey=True, squeeze=False)
    for ax, ds in zip(axes[0], BKG_DATASETS):
        b = load_bins(ds, seed, wave, nth)
        if b is None:
            ax.text(0.5, 0.5, "no bins", transform=ax.transAxes, ha="center")
            continue
        cen = _centers(b["e_sig_edges"])
        m = _mask(b["e_sig"]["n"])
        eff = np.asarray(b["e_sig"]["efficiency"])
        exact = np.asarray(b["e_sig"]["exact"])
        ax.plot(cen[m], eff[m], "o-", color=ACCENT, ms=4, lw=1.6,
                label="with $\\gamma$ bkg: efficiency")
        ax.plot(cen[m], exact[m], "s--", color=ACCENT, ms=3.5, lw=1.2,
                alpha=0.75, label="with $\\gamma$ bkg: exact-1")
        nb = load_bins(signal_dataset(ds), seed, wave, nth)
        if nb is not None:
            m2 = _mask(nb["e_sig"]["n"])
            ax.plot(cen[m2], np.asarray(nb["e_sig"]["efficiency"])[m2],
                    "-", color=WARM, lw=1.8, label="no bkg: efficiency")
            ax.plot(cen[m2], np.asarray(nb["e_sig"]["exact"])[m2],
                    ":", color=WARM, lw=1.6, label="no bkg: exact-1")
        ax.set_xlabel("signal energy (MeV)")
        ax.set_title(ds.replace("_bkg", ""))
        ax.set_ylim(0, 1.04)
        ax.grid(color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    axes[0][0].set_ylabel("fraction of events")
    axes[0][0].legend(fontsize=8, loc="lower right")
    fig.suptitle(f"Event-energy vs reconstruction efficiency "
                 f"($N_{{th}}$={nth}, seed={seed:g}, wave={wave:g} MeV)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def nbkg_performance(rows, seed, wave, nth, seed_low, filename, nth_low=0):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    metrics = [("exact", "exact-1"), ("fake_per_event", "fake clusters / event"),
               ("signal_purity", "signal purity")]
    ds = BKG_DATASETS[0] if BKG_DATASETS else None
    ref = load_bins(ds, seed, wave, nth)
    series = [(seed, nth, "-", "best point"),
              (seed_low, nth_low, "--", f"seed={seed_low:g}, Nth={nth_low}")]
    for ax, (key, label) in zip(axes, metrics):
        for s, n, ls, lab in series:
            b = load_bins(ds, s, wave, n)
            if b is None:
                continue
            x = np.arange(len(b["nbkg_values"]))
            y = np.asarray(b["nbkg"][key], dtype=float)
            m = np.isfinite(y)
            ax.plot(x[m], y[m], "o" + ls, color=ACCENT if s == seed else WARM,
                    ms=4, lw=1.6, label=f"{lab} ({s:g} MeV, Nth={n})")
        ax.set_xticks(range(len(ref["nbkg_values"])) if ref else range(5))
        ax.set_xticklabels([str(v) if v < 4 else "≥4"
                            for v in (ref["nbkg_values"] if ref else [0, 1, 2, 3, 4])])
        ax.set_xlabel("injected background hits $N_{bkg}$")
        ax.set_ylabel(label)
        ax.grid(color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_ylim(0, 1.05)
        ax.legend(fontsize=8)
    fig.suptitle(f"Background multiplicity dependence "
                 f"({ds.replace('_bkg', '') if ds else ''}, "
                 f"wave={wave:g} MeV): best point vs low-seed $N_{{th}}$=0",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def purity_fake_vs_energy(rows, seed, wave, nth, filename):
    fig, axes = plt.subplots(2, len(BKG_DATASETS),
                             figsize=(4.4 * len(BKG_DATASETS), 6.4),
                             sharey="row", squeeze=False)
    for c, ds in enumerate(BKG_DATASETS):
        b = load_bins(ds, seed, wave, nth)
        if b is None:
            continue
        cen = _centers(b["e_sig_edges"])
        m = _mask(b["e_sig"]["n"])
        ax = axes[0][c]
        ax.plot(cen[m], np.asarray(b["e_sig"]["signal_purity"])[m], "o-",
                color=GREEN, ms=4, lw=1.6, label="signal purity")
        ax.plot(cen[m], np.asarray(b["e_sig"]["lead_sig_containment"])[m], "s--",
                color=ACCENT, ms=3.5, lw=1.2, label="signal containment")
        ax.set_ylim(0.7, 1.01)
        ax.set_title(ds.replace("_bkg", ""))
        ax.grid(color=GRID, lw=0.8)
        ax.legend(fontsize=8)
        ax = axes[1][c]
        ax.plot(cen[m], np.asarray(b["e_sig"]["fake_per_event"])[m], "o-",
                color=RED, ms=4, lw=1.6, label="fake clusters / event")
        ax.set_xlabel("signal energy (MeV)")
        ax.set_ylim(0, max(0.05, 1.1 * np.nanmax(
            np.asarray(b["e_sig"]["fake_per_event"])[m])))
        ax.grid(color=GRID, lw=0.8)
        ax.legend(fontsize=8)
    axes[0][0].set_ylabel("fraction")
    axes[1][0].set_ylabel("fake clusters / event")
    fig.suptitle(f"Cluster purity and fake-cluster rate vs signal energy "
                 f"($N_{{th}}$={nth}, seed={seed:g}, wave={wave:g} MeV)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def response_curves(rows, seed, wave, nth, filename):
    fig, axes = plt.subplots(2, len(BKG_DATASETS),
                             figsize=(4.4 * len(BKG_DATASETS), 6.4),
                             sharey="row", squeeze=False)
    for c, ds in enumerate(BKG_DATASETS):
        b = load_bins(ds, seed, wave, nth)
        if b is None:
            continue
        cen = _centers(b["e_sig_edges"])
        m = _mask(b["e_sig"]["n"], 50) & (cen >= 20.0)
        nb = load_bins(signal_dataset(ds), seed, wave, nth)
        ax = axes[0][c]
        ax.plot(cen[m], 100 * np.asarray(b["e_sig"]["bias_total"])[m],
                "o-", color=ACCENT, ms=4, lw=1.6, label="bkg: leading cluster")
        ax.plot(cen[m], 100 * np.asarray(b["e_sig"]["bias_signal"])[m],
                "s--", color=GREEN, ms=3.5, lw=1.3, label="bkg: signal part")
        if nb is not None:
            m2 = _mask(nb["e_sig"]["n"], 50) & (cen >= 20.0)
            ax.plot(cen[m2], 100 * np.asarray(nb["e_sig"]["bias_signal"])[m2],
                    ":", color=WARM, lw=1.5, label="no bkg")
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_title(ds.replace("_bkg", ""))
        ax.set_ylabel("relative bias (%)")
        ax.set_ylim(-15, 5)
        ax.grid(color=GRID, lw=0.8)
        ax.legend(fontsize=8)
        ax = axes[1][c]
        ax.plot(cen[m], 100 * np.asarray(b["e_sig"]["rms_total"])[m],
                "o-", color=ACCENT, ms=4, lw=1.6, label="bkg: leading cluster")
        ax.plot(cen[m], 100 * np.asarray(b["e_sig"]["rms_signal"])[m],
                "s--", color=GREEN, ms=3.5, lw=1.3, label="bkg: signal part")
        if nb is not None:
            m2 = _mask(nb["e_sig"]["n"], 50) & (cen >= 20.0)
            ax.plot(cen[m2], 100 * np.asarray(nb["e_sig"]["rms_signal"])[m2],
                    ":", color=WARM, lw=1.5, label="no bkg")
        ax.set_xlabel("signal energy (MeV)")
        ax.set_ylabel("relative resolution (%)")
        ax.set_ylim(0, 20)
        ax.grid(color=GRID, lw=0.8)
        ax.legend(fontsize=8)
    fig.suptitle(f"Reconstructed-energy response vs signal energy "
                 f"($N_{{th}}$={nth}, seed={seed:g}, wave={wave:g} MeV)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# --------------------------------------------------------------------------- #
# tables
# --------------------------------------------------------------------------- #
def _fmt(v, nd=4):
    try:
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return "—"


def best_tables(rows, nobkg_rows):
    lines = ["# Background scan — best points", ""]
    lines += ["## Best point per neighbour threshold (mean over the three "
              "signal samples)", "",
              "| Nth | seed (MeV) | wave (MeV) | exact-1 | efficiency | "
              "over-split | fake/event | purity | signal capture | "
              "bkg capture | bias (E>20 MeV) |",
              "|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|"]
    for nth in NTH_VALUES:
        best = best_for_nth(rows, nth)
        if not best:
            continue
        s, w, v, at = best
        def avg(k):
            return float(np.nanmean([r.get(k, np.nan) for r in at]))
        lines.append(
            f"| {nth} | {s:g} | {w:g} | {avg('exact_match_rate'):.4f} "
            f"| {avg('cluster_efficiency'):.4f} "
            f"| {avg('over_cluster_rate'):.4f} "
            f"| {avg('fake_clusters_per_event'):.3f} "
            f"| {avg('signal_purity_mean'):.4f} "
            f"| {avg('signal_capture_fraction'):.4f} "
            f"| {avg('bkg_capture_fraction'):.4f} "
            f"| {avg('energy_bias_total_hi'):+.4f} |")

    lines += ["", "## Best point per dataset (any Nth)", "",
              "| dataset | Nth | seed | wave | exact-1 | efficiency | "
              "purity | fake/event | signal capture |",
              "|:--|--:|--:|--:|--:|--:|--:|--:|--:|"]
    for ds in BKG_DATASETS:
        b = best_for_dataset(rows, ds)
        if b is None:
            continue
        lines.append(
            f"| {ds.replace('_bkg', '')} | {int(b['Nth'])} | {b['Eth_seed']:g} "
            f"| {b['Eth_wave']:g} | {b['exact_match_rate']:.4f} "
            f"| {b.get('cluster_efficiency', np.nan):.4f} "
            f"| {b.get('signal_purity_mean', np.nan):.4f} "
            f"| {b.get('fake_clusters_per_event', np.nan):.3f} "
            f"| {b.get('signal_capture_fraction', np.nan):.4f} |")

    if nobkg_rows:
        lines += ["", "## Same-code no-background comparison at the best "
                  "bkg points", "",
                  "| dataset | bkg exact | no-bkg exact | Δ exact | "
                  "bkg purity | no-bkg purity |", "|:--|--:|--:|--:|--:|--:|"]
        for ds in BKG_DATASETS:
            b = best_for_dataset(rows, ds)
            if b is None:
                continue
            # find the closest no-bkg point (same seed/wave/nth)
            nb = [r for r in nobkg_rows
                  if r["dataset"] == signal_dataset(ds)
                  and r["Eth_seed"] == b["Eth_seed"]
                  and r["Eth_wave"] == b["Eth_wave"]
                  and int(r["Nth"]) == int(b["Nth"])]
            if not nb:
                nb = [r for r in nobkg_rows
                      if r["dataset"] == signal_dataset(ds)
                      and int(r["Nth"]) == int(b["Nth"])]
            if nb:
                n = max(nb, key=lambda r: r["exact_match_rate"])
                lines.append(
                    f"| {ds.replace('_bkg', '')} | {b['exact_match_rate']:.4f} "
                    f"| {n['exact_match_rate']:.4f} "
                    f"| {b['exact_match_rate'] - n['exact_match_rate']:+.4f} "
                    f"| {b.get('signal_purity_mean', np.nan):.4f} "
                    f"| {n.get('signal_purity_mean', np.nan):.4f} |")
    return "\n".join(lines) + "\n"


def sample_overview(filename="sample_overview_bkg.png"):
    """Signal/total energy and NBkg distributions of the bkg samples."""
    import ROOT
    from ROOT import vector
    from scipy import stats

    def read(ds):
        f = ROOT.TFile.Open(
            str(SUBSET_DIR / "processed" / f"subset_{ds}.root"), "READ")
        t = f.Get("EdepOfEachEvt")
        e = vector("float")()
        t.SetBranchAddress("Edeps", e)
        nb = array.array("i", [0])
        has = any(b.GetName() == "NBkg" for b in t.GetListOfBranches())
        if has:
            t.SetBranchAddress("NBkg", nb)
        energy, nhit = [], []
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            energy.append(float(np.asarray(list(e), dtype=float).sum()))
            if has:
                nhit.append(int(nb[0]))
        f.Close()
        return np.asarray(energy), np.asarray(nhit, dtype=int)

    ds = BKG_DATASETS[0]
    sig, _ = read(signal_dataset(ds))
    tot, nb = read(ds)
    bkg = tot - sig

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.7))
    ax = axes[0]
    bins = np.arange(0, 82, 2)
    ax.hist(sig, bins=bins, histtype="step", lw=1.6, color=WARM,
            label="signal only")
    ax.hist(tot, bins=bins, histtype="step", lw=1.6, color=ACCENT,
            label="signal + $\\gamma$ bkg")
    ax.set_yscale("log")
    ax.set_xlabel("event total energy (MeV)")
    ax.set_ylabel("events / 2 MeV")
    ax.legend(fontsize=8)
    ax.set_title(ds.replace("_bkg", ""))

    ax = axes[1]
    n = len(nb)
    ax.hist(nb, bins=np.arange(-0.5, 13.5, 1), color=ACCENT, alpha=0.65,
            density=False, label="injected hits")
    k = np.arange(0, 13)
    lam = float(nb.mean())
    ax.plot(k, n * stats.poisson.pmf(k, lam), "o-", color=RED, ms=3,
            lw=1.2, label=f"Poisson ($\\lambda$={lam:.2f})")
    ax.set_xlabel("injected background hits $N_{bkg}$")
    ax.set_ylabel("events")
    ax.legend(fontsize=8)

    ax = axes[2]
    ax.hist(bkg, bins=np.arange(0, 62, 1), color=GREY, alpha=0.8,
            label="background energy")
    ax.axvline(bkg.mean(), color=RED, ls="--", lw=1.3,
               label=f"mean {bkg.mean():.1f} MeV")
    ax.set_xlabel("background energy per event (MeV)")
    ax.set_ylabel("events / MeV")
    ax.legend(fontsize=8)

    for ax in axes:
        ax.grid(color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
    fig.suptitle("Background-injected sample overview", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def parameter_curves(rows_bkg, ref_rows, seed_b, wave_b, nth_b,
                     filename="parameter_curves_bkg.png"):
    """1D projections at the best point: Nth, wave and seed, bkg vs no-bkg."""
    nobkg = [signal_dataset(ds) for ds in BKG_DATASETS]

    def avg(rows, datasets, nth, seed, wave, key):
        vals = [r[key] for r in rows
                if r["dataset"] in datasets and int(r["Nth"]) == nth
                and r["Eth_seed"] == seed and r["Eth_wave"] == wave
                and np.isfinite(r.get(key, np.nan))]
        return float(np.mean(vals)) if vals else np.nan

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9))

    ax = axes[0]
    xs = list(NTH_VALUES)
    yb = [avg(rows_bkg, BKG_DATASETS, n, seed_b, wave_b, "exact_match_rate")
          for n in xs]
    yn = [avg(ref_rows, nobkg, n, seed_b, wave_b, "exact_match_rate")
          for n in xs]
    ax.plot(xs, yb, "o-", color=ACCENT, lw=1.8, ms=5, label="with $\\gamma$ bkg")
    ax.plot(xs, yn, "s--", color=WARM, lw=1.5, ms=4, label="no bkg (archived)")
    ax.set_xticks(xs)
    ax.set_xlabel("$N_{th}$")
    ax.set_ylabel("exact-1")
    ax.set_title(f"neighbour threshold (seed={seed_b:g}, wave={wave_b:g})")
    ax.set_ylim(0.70, 0.90)

    ax = axes[1]
    ws = [w for w in ETH_WAVE_VALUES if w < seed_b]
    yb = [avg(rows_bkg, BKG_DATASETS, nth_b, seed_b, w, "exact_match_rate")
          for w in ws]
    yn = [avg(ref_rows, nobkg, nth_b, seed_b, w, "exact_match_rate")
          for w in ws]
    ax.plot(ws, yb, "o-", color=ACCENT, lw=1.8, ms=4, label="with $\\gamma$ bkg")
    ax.plot(ws, yn, "s--", color=WARM, lw=1.5, ms=4, label="no bkg (archived)")
    ax.set_xscale("log")
    ax.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.tick_params(axis="x", labelsize=8, labelrotation=45)
    for t in ax.get_xticklabels():
        t.set_ha("right")
    ax.set_xlabel("$E_{th}^{\\rm wave}$ (MeV)")
    ax.set_ylabel("exact-1")
    ax.set_title(f"wave threshold ($N_{{th}}$={nth_b}, seed={seed_b:g})")
    ax.set_ylim(0.6, 0.90)

    ax = axes[2]
    ss = list(ETH_SEED_VALUES)
    yb = [avg(rows_bkg, BKG_DATASETS, nth_b, s, wave_b, "exact_match_rate")
          for s in ss]
    yn = [avg(ref_rows, nobkg, nth_b, s, wave_b, "exact_match_rate")
          for s in ss]
    ax.plot(ss, yb, "o-", color=ACCENT, lw=1.8, ms=4, label="with $\\gamma$ bkg")
    ax.plot(ss, yn, "s--", color=WARM, lw=1.5, ms=4, label="no bkg (archived)")
    ax.set_xlabel("$E_{th}^{\\rm seed}$ (MeV)")
    ax.set_ylabel("exact-1")
    ax.set_title(f"seed threshold ($N_{{th}}$={nth_b}, wave={wave_b:g})")
    ax.set_ylim(0.70, 0.90)
    ax.axvspan(10, 12, color=ACCENT, alpha=0.08, lw=0)
    ax2 = ax.twinx()
    fake = [avg(rows_bkg, BKG_DATASETS, nth_b, s, wave_b,
                "fake_clusters_per_event") for s in ss]
    fake0 = ([avg(rows_bkg, BKG_DATASETS, 0, s, wave_b,
                  "fake_clusters_per_event") for s in ss]
             if nth_b != 0 else None)
    ax2.plot(ss, fake, "^:", color=RED, lw=1.3, ms=4,
             label=f"fake, $N_{{th}}$={nth_b}")
    if fake0 is not None:
        ax2.plot(ss, fake0, "v:", color=RED, lw=1.1, ms=3, alpha=0.55,
                 label="fake, $N_{th}$=0")
    ax2.set_ylabel("fake clusters / event", color=RED)
    ax2.tick_params(axis="y", colors=RED)
    top = max(np.nanmax(fake), np.nanmax(fake0) if fake0 is not None else 0)
    ax2.set_ylim(0, max(0.02, 1.15 * top))

    for i, ax in enumerate(axes):
        ax.grid(color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if i == 2:
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
        else:
            ax.legend(fontsize=8, loc="best")
    fig.suptitle("Parameter projections at the best background working point",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    global RESULTS_DIR, SUMMARY_CSV

    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=SUMMARY_CSV)
    ap.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    ap.add_argument("--nobkg-csv", type=Path,
                    default=Path(__file__).resolve().parent.parent
                    / "results/sweep/26-10-16-refined/sweep_summary.csv")
    args = ap.parse_args()

    RESULTS_DIR = args.results_dir
    SUMMARY_CSV = args.csv
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if not BKG_DATASETS:
        print("no *_bkg datasets in ECAL_SWEEP_DATASETS; nothing to do")
        return

    rows = load_rows(SUMMARY_CSV)
    rows_bkg = [r for r in rows if int(r.get("has_bkg", 0)) == 1]
    rows_nobkg = [r for r in rows if int(r.get("has_bkg", 0)) == 0]
    print(f"loaded {len(rows)} rows: {len(rows_bkg)} bkg, "
          f"{len(rows_nobkg)} same-code no-bkg")

    # Archived no-background campaign for the full comparison heatmap/table;
    # the same-code control (if run) is used for the energy-efficiency overlay
    # because its bins live in the current TASK_DIR.
    control_rows = rows_nobkg
    ref_rows = []
    ref_label = "none"
    if args.nobkg_csv.exists():
        ref_rows = load_rows(args.nobkg_csv)
        ref_label = "26-10-16-refined (archived)"
    elif control_rows:
        ref_rows = control_rows
        ref_label = "same-code control"
    print(f"no-bkg reference: {ref_label} ({len(ref_rows)} rows); "
          f"same-code control rows: {len(control_rows)}")

    best = best_for_nth(rows, 0)
    if best is None:
        print("no bkg points with Nth=0; run the scan first")
        return
    seed0, wave0, _, _ = best
    # best Nth/seed/wave overall
    best_all = None
    for nth in NTH_VALUES:
        b = best_for_nth(rows, nth)
        if b and (best_all is None or b[2] > best_all[2]):
            best_all = b
    seed_b, wave_b, _, _ = best_all
    nth_b = [n for n in NTH_VALUES
             if best_for_nth(rows, n) and
             abs(best_for_nth(rows, n)[2] - best_all[2]) < 1e-12]
    nth_b = nth_b[0] if nth_b else 0
    print(f"best bkg point: Nth={nth_b}, seed={seed_b:g}, wave={wave_b:g} "
          f"(mean exact={best_all[2]:.4f})")

    sample_overview()

    # main heatmaps
    zoom = [i for i, w in enumerate(ETH_WAVE_VALUES) if w <= 1.0]
    heatmap_per_nth(rows_bkg, BKG_DATASETS, "exact_match_rate",
                    "Exact-match rate with $\\gamma$ background "
                    "(rows: dataset; columns: $N_{th}$)",
                    "heatmap_exact_match_rate_bkg.png",
                    wave_idx=zoom, cell_w=0.75, val_fs=6.0)
    heatmap_per_nth(rows_bkg, BKG_DATASETS, "exact_match_rate",
                    "Exact-match rate with $\\gamma$ background "
                    "[full wave range]",
                    "heatmap_exact_match_rate_bkg_full.png",
                    wave_idx=None, annotate=False, cell_w=0.30)
    for key, title, fn, cmap, vmin, vmax, fmt in [
        ("fake_clusters_per_event", "Fake (background-like) clusters per event",
         "heatmap_fake_clusters_per_event_bkg.png", "RdYlGn_r", 0, 0.5, "{:.3f}"),
        ("signal_purity_mean", "Signal purity of the clusters",
         "heatmap_signal_purity_bkg.png", "RdYlGn", 0.9, 1.0, "{:.3f}"),
        ("signal_capture_fraction", "Signal energy capture",
         "heatmap_signal_capture_bkg.png", "RdYlGn", 0.8, 1.0, "{:.3f}"),
        ("over_cluster_rate", "Over-clustering rate",
         "heatmap_over_cluster_rate_bkg.png", "RdYlGn_r", 0, 0.2, "{:.3f}"),
    ]:
        heatmap_per_nth(rows_bkg, BKG_DATASETS, key, title, fn,
                        cmap=cmap, vmin=vmin, vmax=vmax, fmt=fmt,
                        wave_idx=zoom, cell_w=0.75, val_fs=6.0)

    if ref_rows:
        comparison_heatmap(
            rows_bkg, ref_rows, "exact_match_rate",
            "Δ exact-match (with $\\gamma$ bkg − no bkg), same parameters",
            "comparison_exact_bkg_minus_nobkg.png",
            vmin=-0.15, vmax=0.15, wave_idx=zoom, annotate=True)
        if ref_rows and "signal_capture_fraction" in ref_rows[0]:
            comparison_heatmap(
                rows_bkg, ref_rows, "signal_capture_fraction",
                "Δ signal capture (with $\\gamma$ bkg − no bkg)",
                "comparison_signal_capture_bkg_minus_nobkg.png",
                cmap="RdBu_r", vmin=-0.05, vmax=0.05,
                fmt="{:+.3f}", wave_idx=zoom, annotate=True)

    parameter_curves(rows_bkg, ref_rows, seed_b, wave_b, nth_b,
                     "parameter_curves_bkg.png")
    energy_efficiency(rows_bkg, seed_b, wave_b, nth_b,
                      "energy_efficiency_bkg.png")
    nbkg_performance(rows_bkg, seed_b, wave_b, nth_b, 7,
                     "nbkg_performance_bkg.png")
    purity_fake_vs_energy(rows_bkg, seed_b, wave_b, nth_b,
                          "purity_fake_vs_energy_bkg.png")
    response_curves(rows_bkg, seed_b, wave_b, nth_b,
                    "response_bkg.png")

    text = best_tables(rows_bkg, ref_rows)
    (RESULTS_DIR / "best_params_bkg.md").write_text(text)
    print(f"wrote {RESULTS_DIR / 'best_params_bkg.md'}")
    print(text)


if __name__ == "__main__":
    main()
