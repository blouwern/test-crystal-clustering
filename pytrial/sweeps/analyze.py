"""Aggregate the sweep CSV into heatmaps, trend plots and best-parameter tables.

    python sweeps/analyze.py

All figures are written next to the CSV under ``pytrial/results/sweep/`` so the
summary document can reference them with relative links.
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sweeps.config import (  # noqa: E402
    DATASETS, ETH_VALUES, EXPECTED_LAMBDA, NTH_VALUES, RESULTS_DIR, SUMMARY_CSV,
)

NTH_COLORS = {1.0: "#1f77b4", 2.0: "#d62728", 3.0: "#2ca02c"}
ETH_LABELS = [f"{e:g}" for e in ETH_VALUES]
NTH_LABELS = [f"N={n:g}" for n in NTH_VALUES]


def load_rows():
    with open(SUMMARY_CSV, newline="") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        for key, value in list(row.items()):
            if key == "dataset":
                continue
            try:
                row[key] = float(value)
            except (TypeError, ValueError):
                pass
    return rows


def matrix(rows, dataset, key):
    mat = np.full((len(ETH_VALUES), len(NTH_VALUES)), np.nan)
    for row in rows:
        if row["dataset"] != dataset:
            continue
        i = ETH_VALUES.index(row["Eth"])
        j = NTH_VALUES.index(row["Nth"])
        mat[i, j] = row.get(key, np.nan)
    return mat


def mean_matrix(rows, key):
    stack = np.stack([matrix(rows, ds, key) for ds in DATASETS])
    return np.nanmean(stack, axis=0)


def heatmap_grid(rows, key, title, filename, cmap="RdYlGn", fmt="{:.3f}",
                 vmin=None, vmax=None):
    panels = [(ds, matrix(rows, ds, key)) for ds in DATASETS]
    panels.append(("mean over datasets", mean_matrix(rows, key)))

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5))
    for ax, (label, mat) in zip(axes.ravel(), panels):
        im = ax.imshow(mat, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(NTH_VALUES)), NTH_LABELS)
        ax.set_yticks(range(len(ETH_VALUES)), ETH_LABELS)
        ax.set_xlabel(r"$N_{th}$ (lit neighbours required)")
        ax.set_ylabel(r"$E_{th}$ (MeV)")
        ax.set_title(label)
        for i in range(len(ETH_VALUES)):
            for j in range(len(NTH_VALUES)):
                val = mat[i, j]
                if np.isnan(val):
                    continue
                ax.text(j, i, fmt.format(val), ha="center", va="center", fontsize=8.5)
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.suptitle(title, fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def trend_plot(rows, key, ylabel, filename, ylim=None, hline=None):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), sharey=True)
    for ax, ds in zip(axes, DATASETS):
        for nth in NTH_VALUES:
            ys = [next((r[key] for r in rows
                        if r["dataset"] == ds and r["Eth"] == eth and r["Nth"] == nth),
                       np.nan)
                  for eth in ETH_VALUES]
            ax.plot(ETH_VALUES, ys, marker="o", ms=4, color=NTH_COLORS[nth],
                    label=f"Nth={nth:g}")
        if hline is not None:
            ax.axhline(hline, color="grey", ls="--", lw=1)
        ax.set_title(ds)
        ax.set_xlabel(r"$E_{th}$ (MeV)")
        ax.set_xticks(ETH_VALUES, ETH_LABELS, rotation=45)
        ax.grid(alpha=0.3)
        if ylim:
            ax.set_ylim(*ylim)
    axes[0].set_ylabel(ylabel)
    axes[0].legend()
    fig.tight_layout()
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def composition_plot(rows, nth, filename):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), sharey=True)
    x = np.arange(len(ETH_VALUES))
    for ax, ds in zip(axes, DATASETS):
        sel = [next((r for r in rows if r["dataset"] == ds and r["Eth"] == eth
                     and r["Nth"] == nth), None) for eth in ETH_VALUES]
        f0 = [s["frac_0"] if s else np.nan for s in sel]
        f1 = [s["frac_1"] if s else np.nan for s in sel]
        f2 = [s["frac_2"] if s else np.nan for s in sel]
        f3 = [s["frac_3plus"] if s else np.nan for s in sel]
        ax.bar(x, f0, color="#4c72b0", label="0 clusters (missed)")
        ax.bar(x, f1, bottom=f0, color="#55a868", label="1 cluster (correct)")
        ax.bar(x, f2, bottom=np.array(f0) + np.array(f1), color="#c44e52",
               label="2 clusters")
        ax.bar(x, f3, bottom=np.array(f0) + np.array(f1) + np.array(f2),
               color="#8172b2", label=r"$\geq$3 clusters")
        ax.set_xticks(x, ETH_LABELS, rotation=45)
        ax.set_xlabel(r"$E_{th}$ (MeV)")
        ax.set_ylim(0, 1)
        ax.set_title(f"{ds}  (Nth={nth:g})")
        ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("event fraction")
    axes[0].legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def dataset_comparison(rows, key, ylabel, filename):
    fig, ax = plt.subplots(figsize=(8.5, 5))
    width = 0.25
    x = np.arange(len(ETH_VALUES))
    for k, ds in enumerate(DATASETS):
        ys = [next((r[key] for r in rows if r["dataset"] == ds and r["Eth"] == eth
                    and r["Nth"] == 1.0), np.nan) for eth in ETH_VALUES]
        ax.bar(x + (k - 1) * width, ys, width, label=ds)
    ax.set_xticks(x, ETH_LABELS)
    ax.set_xlabel(r"$E_{th}$ (MeV)")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def _mean_over_datasets(rows, key, eth, nth):
    vals = [r[key] for r in rows if r["Eth"] == eth and r["Nth"] == nth]
    return float(np.mean(vals)) if vals else np.nan


def tradeoff_plot(rows, nth, filename):
    """The central result: multiplicity quality and energy quality pull apart."""
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    exact = [_mean_over_datasets(rows, "exact_match_rate", e, nth) for e in ETH_VALUES]
    cont = [_mean_over_datasets(rows, "energy_containment_mean", e, nth) for e in ETH_VALUES]
    capt = [_mean_over_datasets(rows, "energy_captured_fraction", e, nth) for e in ETH_VALUES]
    miss = [_mean_over_datasets(rows, "frac_0", e, nth) for e in ETH_VALUES]

    ax.plot(ETH_VALUES, exact, "o-", color="#2ca02c", lw=2,
            label="exact-match rate (multiplicity correct)")
    ax.plot(ETH_VALUES, capt, "^-", color="#1f77b4", lw=2,
            label="energy capture (all clusters)")
    ax.plot(ETH_VALUES, cont, "s-", color="#d62728", lw=2,
            label="energy containment (leading cluster)")
    ax.plot(ETH_VALUES, miss, "v--", color="#7f7f7f", lw=1.5,
            label="events with 0 clusters (missed)")
    ax.set_xlabel(r"$E_{th}$ (MeV)")
    ax.set_ylabel("fraction")
    ax.set_xticks(ETH_VALUES, ETH_LABELS, rotation=45)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(loc="center right", fontsize=9)
    ax.set_title(f"Threshold trade-off (mean over e+/gamma/mu-, Nth={nth:g})")
    fig.tight_layout()
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def energy_decomposition(rows, nth, filename):
    """Where the event energy goes, as a function of the threshold.

    Uses a ratio-of-sums decomposition so the three shares add up to exactly 1:
    leading = captured - split, and captured + unclustered = 1.
    """
    capt = np.array([_mean_over_datasets(rows, "energy_captured_fraction", e, nth)
                     for e in ETH_VALUES])
    split = np.array([_mean_over_datasets(rows, "energy_split_fraction", e, nth)
                      for e in ETH_VALUES])
    uncl = np.array([_mean_over_datasets(rows, "energy_unclustered_fraction", e, nth)
                     for e in ETH_VALUES])
    lead = capt - split

    fig, ax = plt.subplots(figsize=(10, 5.4))
    x = np.arange(len(ETH_VALUES))
    ax.bar(x, lead, color="#55a868", label="leading cluster (wanted)")
    ax.bar(x, split, bottom=lead, color="#c44e52",
           label="other clusters (over-splitting)")
    ax.bar(x, uncl, bottom=lead + split, color="#8172b2",
           label="no cluster: below threshold (energy loss)")
    ax.set_xticks(x, ETH_LABELS)
    ax.set_xlabel(r"$E_{th}$ (MeV)")
    ax.set_ylabel("fraction of event energy")
    ax.set_ylim(0, 1)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(loc="lower left", fontsize=9)
    ax.set_title(f"Event-energy budget vs threshold (mean over datasets, Nth={nth:g})")
    fig.tight_layout()
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def nth_effect_plot(rows, filename):
    """How much does the neighbour threshold actually matter?"""
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for ax, key, label in ((axes[0], "exact_match_rate", "exact-match rate"),
                           (axes[1], "energy_containment_mean", "energy containment")):
        for ds in DATASETS:
            for nth in NTH_VALUES:
                ys = [next((r[key] for r in rows if r["dataset"] == ds
                            and r["Eth"] == eth and r["Nth"] == nth), np.nan)
                      for eth in ETH_VALUES]
                ax.plot(ETH_VALUES, ys, marker="o", ms=3, color=NTH_COLORS[nth],
                        ls="-" if ds == "ep_iso" else "--",
                        alpha=0.9 if ds == "ep_iso" else 0.6,
                        label=f"{ds} Nth={nth:g}")
        ax.set_xlabel(r"$E_{th}$ (MeV)")
        ax.set_ylabel(label)
        ax.set_xticks(ETH_VALUES, ETH_LABELS, rotation=45)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    fig.suptitle(r"Effect of $N_{th}$: curves for the same dataset nearly coincide")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def best_tables(rows):
    lines = []
    keys = ["exact_match_rate", "energy_containment_mean", "MAE", "mean_clusters"]
    lines.append("| metric | dataset | best Eth | best Nth | value |")
    lines.append("|:--|:--|--:|--:|--:|")
    for key in keys:
        for ds in DATASETS + ["__mean__"]:
            subset = rows if ds == "__mean__" else [r for r in rows if r["dataset"] == ds]
            if key == "exact_match_rate":
                best = max(subset, key=lambda r: r.get(key, -np.inf))
            elif key == "energy_containment_mean":
                best = max(subset, key=lambda r: r.get(key, -np.inf))
            elif key == "MAE":
                best = min(subset, key=lambda r: r.get(key, np.inf))
            else:
                best = min(subset, key=lambda r: abs(r.get(key, np.inf) - EXPECTED_LAMBDA))
            ds_label = "**mean**" if ds == "__mean__" else ds
            lines.append(f"| {key} | {ds_label} | {best['Eth']:g} | {best['Nth']:g} "
                         f"| {best[key]:.4f} |")
    return "\n".join(lines)


def ranking_table(rows):
    """Average rank across the headline metrics (lower = better), per grid point."""
    metrics = [("exact_match_rate", False), ("energy_containment_mean", False),
               ("MAE", True)]
    points = {}
    for ds in DATASETS:
        for eth in ETH_VALUES:
            for nth in NTH_VALUES:
                rec = next(r for r in rows if r["dataset"] == ds
                           and r["Eth"] == eth and r["Nth"] == nth)
                points.setdefault((eth, nth), {})[ds] = rec

    agg = {}
    for (eth, nth), per_ds in points.items():
        score = 0.0
        for key, lower_better in metrics:
            for ds, rec in per_ds.items():
                values = sorted(r[key] for r in rows)
                rank = values.index(rec[key]) / (len(values) - 1)
                score += rank if not lower_better else (1.0 - rank)
        agg[(eth, nth)] = score / (len(metrics) * len(DATASETS))

    lines = ["| rank | Eth | Nth | mean normalised score | exact match (mean) | containment (mean) | MAE (mean) |",
             "|--:|--:|--:|--:|--:|--:|--:|"]
    ordered = sorted(agg.items(), key=lambda kv: -kv[1])
    for i, ((eth, nth), score) in enumerate(ordered, 1):
        per_ds = points[(eth, nth)]
        em = np.mean([per_ds[d]["exact_match_rate"] for d in DATASETS])
        ec = np.mean([per_ds[d]["energy_containment_mean"] for d in DATASETS])
        mae = np.mean([per_ds[d]["MAE"] for d in DATASETS])
        lines.append(f"| {i} | {eth:g} | {nth:g} | {score:.3f} | {em:.4f} | {ec:.4f} | {mae:.4f} |")
    return "\n".join(lines)


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    if not rows:
        print("no sweep rows found - run run_sweep.py first")
        return

    heatmap_grid(rows, "exact_match_rate",
                 "Exact-match rate: fraction of events with exactly 1 cluster",
                 "heatmap_exact_match.png", cmap="RdYlGn", vmin=0.2, vmax=0.7)
    heatmap_grid(rows, "energy_containment_mean",
                 "Energy containment: mean $E_{\\mathrm{leading\\ cluster}}/E_{\\mathrm{event}}$",
                 "heatmap_energy_containment.png", cmap="RdYlGn", vmin=0.6, vmax=0.95)
    heatmap_grid(rows, "mean_clusters",
                 "Mean number of reconstructed clusters (truth = 1)",
                 "heatmap_mean_clusters.png", cmap="RdYlBu_r", vmin=0.9, vmax=2.0)
    heatmap_grid(rows, "MAE",
                 "MAE of cluster multiplicity vs truth = 1",
                 "heatmap_mae.png", cmap="RdYlGn_r", vmin=0.3, vmax=1.5)
    heatmap_grid(rows, "energy_captured_fraction",
                 "Energy capture: fraction of event energy inside any cluster",
                 "heatmap_energy_captured.png", cmap="RdYlGn", vmin=0.55, vmax=0.92)
    heatmap_grid(rows, "energy_unclustered_fraction",
                 "Energy loss: fraction of event energy below threshold (unclustered)",
                 "heatmap_energy_unclustered.png", cmap="RdYlGn_r", vmin=0.08, vmax=0.45)

    trend_plot(rows, "exact_match_rate", "exact-match rate",
               "trend_exact_match.png", hline=0.5)
    trend_plot(rows, "energy_containment_mean", "energy containment",
               "trend_energy_containment.png", hline=0.9)
    trend_plot(rows, "mean_clusters", "mean clusters",
               "trend_mean_clusters.png", hline=EXPECTED_LAMBDA)
    trend_plot(rows, "frac_0", "fraction of events with 0 clusters",
               "trend_missed_events.png")

    composition_plot(rows, 1.0, "composition_nth1.png")
    composition_plot(rows, 2.0, "composition_nth2.png")
    dataset_comparison(rows, "exact_match_rate", "exact-match rate",
                       "compare_datasets_exact_match.png")
    tradeoff_plot(rows, 1.0, "tradeoff_nth1.png")
    tradeoff_plot(rows, 2.0, "tradeoff_nth2.png")
    energy_decomposition(rows, 1.0, "energy_budget_nth1.png")
    energy_decomposition(rows, 2.0, "energy_budget_nth2.png")
    nth_effect_plot(rows, "nth_effect.png")

    print("\n### best per metric\n")
    print(best_tables(rows))
    print("\n### grid ranking\n")
    print(ranking_table(rows))

    (RESULTS_DIR / "best_params.md").write_text(
        "# Best parameters per metric\n\n" + best_tables(rows)
        + "\n\n# Grid ranking (mean normalised score over datasets and metrics)\n\n"
        + ranking_table(rows) + "\n"
    )
    print(f"\nwrote {RESULTS_DIR / 'best_params.md'}")


if __name__ == "__main__":
    main()
