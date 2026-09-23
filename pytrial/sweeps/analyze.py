"""Aggregate the two-threshold sweep CSV into heatmaps, a trade-off scatter and
best-parameter tables.

    python sweeps/analyze.py
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
    DATASETS, ETH_SEED_VALUES, ETH_WAVE_VALUES, EXPECTED_LAMBDA, RESULTS_DIR,
    SUMMARY_CSV, valid_pairs,
)

SEED_LABELS = [f"{s:g}" for s in ETH_SEED_VALUES]
WAVE_LABELS = [f"{w:g}" for w in ETH_WAVE_VALUES]


def load_rows():
    with open(SUMMARY_CSV, newline="") as f:
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


def matrix(rows, dataset, key):
    """seed x wave matrix; invalid (wave >= seed) pairs stay NaN."""
    m = np.full((len(ETH_SEED_VALUES), len(ETH_WAVE_VALUES)), np.nan)
    for r in rows:
        if dataset is not None and r["dataset"] != dataset:
            continue
        if (r["Eth_seed"], r["Eth_wave"]) not in valid_pairs():
            continue
        i = ETH_SEED_VALUES.index(r["Eth_seed"])
        j = ETH_WAVE_VALUES.index(r["Eth_wave"])
        if dataset is None:  # average over datasets
            prev = m[i, j]
            m[i, j] = r[key] if np.isnan(prev) else np.nanmean([prev, r[key]])
        else:
            m[i, j] = r[key]
    return m


def mean_matrix(rows, key):
    mats = [matrix(rows, ds, key) for ds in DATASETS]
    return np.nanmean(np.stack(mats), axis=0)


def heatmap(rows, key, title, filename, cmap="RdYlGn", fmt="{:.3f}",
            vmin=None, vmax=None, panels="datasets"):
    if panels == "datasets":
        tiles = [(ds, matrix(rows, ds, key)) for ds in DATASETS]
        tiles.append(("mean over datasets", mean_matrix(rows, key)))
    else:
        tiles = [("mean over datasets", mean_matrix(rows, key))]
    n = len(tiles)
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2 * ncol, 5.4 * nrow),
                             squeeze=False)
    for ax, (label, mat) in zip(axes.ravel(), tiles):
        im = ax.imshow(mat, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(ETH_WAVE_VALUES)), WAVE_LABELS)
        ax.set_yticks(range(len(ETH_SEED_VALUES)), SEED_LABELS)
        ax.set_xlabel(r"$E_{th}^{\rm wave}$ (MeV)   —   looser threshold")
        ax.set_ylabel(r"$E_{th}^{\rm seed}$ (MeV)   —   tighter threshold")
        ax.set_title(label)
        for i in range(len(ETH_SEED_VALUES)):
            for j in range(len(ETH_WAVE_VALUES)):
                v = mat[i, j]
                if not np.isnan(v):
                    ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def tradeoff_scatter(rows, filename):
    fig, axes = plt.subplots(1, 4, figsize=(21, 5), sharex=True, sharey=True)
    for ax, ds in zip(axes, DATASETS + [None]):
        sel = [r for r in rows if ds is None or r["dataset"] == ds]
        if ds is not None:
            sel = [r for r in rows if r["dataset"] == ds]
        x = [r.get("energy_containment_mean", np.nan) for r in sel]
        y = [r.get("exact_match_rate", np.nan) for r in sel]
        c = [r["Eth_seed"] for r in sel]
        sc = ax.scatter(x, y, c=c, cmap="viridis", s=45, edgecolors="k", linewidths=0.4)
        ax.set_title(ds or "mean over datasets")
        ax.set_xlabel("energy containment (leading cluster)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("exact-match rate")
    fig.colorbar(sc, ax=axes[-1], label=r"$E_{th}^{\rm seed}$ (MeV)")
    fig.suptitle("Multiplicity quality vs energy containment, one point per "
                 "(seed, wave, dataset)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def best_tables(rows):
    lines = ["| metric | dataset | best seed | best wave | value |", "|:--|:--|--:|--:|--:|"]
    for key, mode in [("exact_match_rate", "max"), ("energy_containment_mean", "max"),
                      ("energy_captured_fraction", "max"), ("MAE", "min")]:
        for ds in DATASETS + [None]:
            sub = [r for r in rows if ds is None or r["dataset"] == ds]
            if not sub:
                continue
            pick = (max if mode == "max" else min)(sub, key=lambda r: r.get(key, np.nan))
            label = ds or "**mean**"
            lines.append(f"| {key} | {label} | {pick['Eth_seed']:g} | {pick['Eth_wave']:g} "
                         f"| {pick[key]:.4f} |")
    return "\n".join(lines)


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    if not rows:
        print("no rows found - run run_sweep.py first")
        return

    heatmap(rows, "exact_match_rate",
            "Exact-match rate (fraction of events with exactly 1 cluster)",
            "heatmap_exact_match.png", vmin=0, vmax=1)
    heatmap(rows, "energy_containment_mean",
            "Energy containment (leading cluster / event energy)",
            "heatmap_energy_containment.png", vmin=0.5, vmax=1.0)
    heatmap(rows, "energy_captured_fraction",
            "Energy capture (all clusters / event energy)",
            "heatmap_energy_captured.png", vmin=0.5, vmax=1.0)
    heatmap(rows, "mean_clusters",
            "Mean number of clusters (truth = 1)",
            "heatmap_mean_clusters.png", cmap="RdYlBu_r", vmin=0.5, vmax=2.5)
    heatmap(rows, "frac_0",
            "Fraction of events with 0 clusters (seed never satisfied)",
            "heatmap_frac0.png", cmap="RdYlGn_r", vmin=0, vmax=0.6)
    tradeoff_scatter(rows, "tradeoff_scatter.png")

    text = ("# Best parameters per metric (two-threshold sweep)\n\n"
            + best_tables(rows) + "\n")
    (RESULTS_DIR / "best_params.md").write_text(text)
    print("\n" + best_tables(rows))
    print(f"\nwrote {RESULTS_DIR / 'best_params.md'}")


if __name__ == "__main__":
    main()
