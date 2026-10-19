"""Aggregate the two-threshold sweep CSV into per-``Nth`` heatmaps, trade-off
scatter plots, wave-threshold profiles and best-parameter tables.

    python sweeps/analyze.py

Layout of the main figures
--------------------------
The first sweep averaged every metric over ``Nth``, which hid the fact that the
neighbour requirement is the dominant parameter and that ``Eth_wave`` only
matters once ``Nth >= 1``.  Every heatmap produced here therefore keeps ``Nth``
as an explicit panel column:

    rows    = dataset (ep_iso / gamma_iso / mu_rest) + "mean over datasets"
    columns = Nth value (0, 1, 2, ...)

so a single figure shows, for one metric, all datasets x all neighbour
thresholds on a shared colour scale.  Filled contours are overlaid on top of
the cell annotations so the metric landscape can be read without squinting at
numbers.
"""

import argparse
import csv
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sweeps.config import (  # noqa: E402
    DATASETS, ETH_SEED_VALUES, ETH_WAVE_VALUES, NTH_VALUES, RESULTS_DIR,
    SUMMARY_CSV,
)

SEED_LABELS = [f"{s:g}" for s in ETH_SEED_VALUES]
WAVE_LABELS = [f"{w:g}" for w in ETH_WAVE_VALUES]

# metrics that the sweep always records; used for graceful degradation when an
# older CSV (from a previous campaign) is re-analysed
REQUIRED_METRICS = ("exact_match_rate", "energy_captured_fraction")

# (csv key, figure title, colour map, vmin, vmax, value format)
HEATMAP_METRICS = [
    ("exact_match_rate",
     "Exact-match rate — fraction of events with exactly 1 cluster (truth = 1)",
     "RdYlGn", 0.0, 1.0, "{:.3f}"),
    ("energy_captured_fraction",
     "Energy capture — clustered energy / event energy (fraction of total)",
     "RdYlGn", 0.0, 1.0, "{:.3f}"),
    ("energy_containment_mean",
     "Energy containment — leading cluster / event energy (fraction of total)",
     "RdYlGn", 0.0, 1.0, "{:.3f}"),
    ("lit_energy_missed_fraction",
     "Missed lit energy — energy above Eth_wave left out of every cluster",
     "RdYlGn_r", 0.0, 1.0, "{:.3f}"),
    ("mean_clusters",
     "Mean number of clusters per event (truth = 1)",
     "RdYlBu_r", 0.0, 3.0, "{:.2f}"),
    ("frac_0",
     "Fraction of events with 0 clusters (no seed ever satisfied)",
     "RdYlGn_r", 0.0, 1.0, "{:.3f}"),
    ("frac_3plus",
     "Fraction of events with >= 3 clusters (shower over-split)",
     "RdYlGn_r", 0.0, 0.5, "{:.3f}"),
]


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


def metric_matrix(rows, dataset, key):
    """seed x wave matrix for one dataset; missing pairs stay NaN."""
    m = np.full((len(ETH_SEED_VALUES), len(ETH_WAVE_VALUES)), np.nan)
    for r in rows:
        if dataset is not None and r["dataset"] != dataset:
            continue
        if r["Eth_seed"] not in ETH_SEED_VALUES or r["Eth_wave"] not in ETH_WAVE_VALUES:
            continue
        if r["Eth_wave"] >= r["Eth_seed"]:
            continue
        i = ETH_SEED_VALUES.index(r["Eth_seed"])
        j = ETH_WAVE_VALUES.index(r["Eth_wave"])
        m[i, j] = r.get(key, np.nan)
    return m


def mean_matrix(rows, key):
    mats = [metric_matrix(rows, ds, key) for ds in DATASETS]
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(np.stack(mats), axis=0)


def available_nths(rows):
    present = {int(r["Nth"]) for r in rows if "Nth" in r}
    ordered = [n for n in NTH_VALUES if n in present]
    return ordered + sorted(present - set(ordered))


# --------------------------------------------------------------------------- #
# per-Nth heatmaps
# --------------------------------------------------------------------------- #
def slice_matrix(mat, wave_idx):
    """Keep only the selected wave columns (``None`` = full grid)."""
    return mat if wave_idx is None else mat[:, list(wave_idx)]


def _draw_matrix(ax, mat, fmt, cmap, vmin, vmax, show_ylabel, wave=None,
                 annotate=True, tick_fs=5.5, val_fs=5.2, label_stride=1):
    wave = list(range(len(ETH_WAVE_VALUES))) if wave is None else list(wave)
    im = ax.imshow(mat, cmap=cmap, aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(wave)),
                  [WAVE_LABELS[j] if k % label_stride == 0 else ""
                   for k, j in enumerate(wave)], fontsize=tick_fs)
    if len(wave) >= 8:
        # dense categorical axis: rotate so 0.02 / 0.05 do not collide
        for t in ax.get_xticklabels():
            t.set_rotation(45)
            t.set_ha("right")
        ax.tick_params(axis="x", pad=1)
    ax.set_yticks(range(len(ETH_SEED_VALUES)), SEED_LABELS, fontsize=6)
    ax.set_xlabel(r"$E_{th}^{\rm wave}$ (MeV)", fontsize=8)
    if show_ylabel:
        ax.set_ylabel(r"$E_{th}^{\rm seed}$ (MeV)", fontsize=8)

    finite = np.isfinite(mat)
    if not annotate and finite.any() and finite.sum() > 1 and len(wave) > 1:
        # overview view: gap-free colour field + iso-levels carry the detail,
        # since individual cells are too narrow for readable annotations
        yi = np.arange(len(ETH_SEED_VALUES))
        xi = np.arange(len(wave))
        filled = np.where(finite, mat, np.nan)
        try:
            ax.contourf(xi, yi, filled, levels=10, cmap=cmap,
                        vmin=vmin, vmax=vmax, extend="both")
            cs = ax.contour(xi, yi, filled, levels=5, colors="k",
                            linewidths=0.4, alpha=0.45)
            ax.clabel(cs, inline=True, fontsize=4.5, fmt=fmt)
        except (ValueError, TypeError):
            pass

    if annotate:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isfinite(v):
                    ax.text(j, i, fmt.format(v), ha="center", va="center",
                            fontsize=val_fs, color="black")
    return im


def _heat_grid(nrow, ncol, n_wave, cell_w, cell_h, cbar_w=0.09, spacer_w=1.0):
    """GridSpec with (matrix, colorbar, spacer) column triplets."""
    total_cols = ncol * 3
    width_ratios = []
    for _ in range(ncol):
        width_ratios += [cell_w, cbar_w, spacer_w]
    fig_w = n_wave * ncol * cell_w + (cbar_w + spacer_w) * ncol + 0.8
    fig = plt.figure(figsize=(fig_w, cell_h * nrow + 1.4))
    gs = GridSpec(nrow, total_cols, figure=fig, width_ratios=width_ratios,
                  wspace=0.06, hspace=0.40)
    return fig, gs


def heatmap_per_nth(rows, key, title, filename, cmap="RdYlGn", fmt="{:.3f}",
                    vmin=None, vmax=None, wave_idx=None, annotate=True,
                    cell_w=0.30, cell_h=1.5, val_fs=5.2, tick_fs=5.5,
                    label_stride=1):
    """One figure per metric: rows = dataset (+mean), columns = Nth.

    Each cell is a full ``Eth_seed`` x ``Eth_wave`` heatmap on a shared colour
    scale, so the neighbour threshold can be read off horizontally and the
    datasets compared vertically.  ``wave_idx`` restricts the wave axis (used
    to zoom into the low-threshold region where the digits stay readable).
    """
    nths = available_nths(rows)
    row_labels = DATASETS + ["mean over datasets"]
    nrow, ncol = len(row_labels), len(nths)
    n_wave = len(ETH_WAVE_VALUES) if wave_idx is None else len(list(wave_idx))
    fig, gs = _heat_grid(nrow, ncol, n_wave, cell_w=cell_w, cell_h=cell_h)

    for r, rlab in enumerate(row_labels):
        for c, nth in enumerate(nths):
            sub = [x for x in rows if int(x["Nth"]) == nth]
            mat = mean_matrix(sub, key) if rlab == "mean over datasets" \
                else metric_matrix(sub, rlab, key)
            ax = fig.add_subplot(gs[r, 3 * c])
            im = _draw_matrix(ax, slice_matrix(mat, wave_idx), fmt, cmap,
                              vmin, vmax, show_ylabel=(c == 0), wave=wave_idx,
                              annotate=annotate, val_fs=val_fs, tick_fs=tick_fs,
                              label_stride=label_stride)
            if c == 0:
                # dataset name lives in the first panel's title, so it can never
                # be clipped by the colourbar / figure edge
                ax.set_title(f"{rlab}  |  $N_{{th}}$ = {nth}", fontsize=9.5)
            elif r == 0:
                ax.set_title(f"$N_{{th}}$ = {nth}", fontsize=11)
            if r != nrow - 1:
                # the wave axis is identical in every row: drop the ticks above
                # the bottom row so titles get the vertical space they need
                ax.set_xticks([])
                ax.set_xlabel("")
            cax = fig.add_subplot(gs[r, 3 * c + 1])
            fig.colorbar(im, cax=cax)

    fig.suptitle(title + "\n"
                 r"(each panel: $E_{th}^{\rm seed}$ (y) $\times$ "
                 r"$E_{th}^{\rm wave}$ (x), shared colour scale)", fontsize=13)
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def heatmap_single_nth(rows, key, nth, title, filename, cmap="RdYlGn",
                       fmt="{:.3f}", vmin=None, vmax=None, wave_idx=None):
    """Compact seed x wave heatmap for one Nth (one column per dataset + mean)."""
    sub = [x for x in rows if int(x["Nth"]) == nth]
    cols = DATASETS + ["mean over datasets"]
    n_wave = len(ETH_WAVE_VALUES) if wave_idx is None else len(list(wave_idx))
    fig, gs = _heat_grid(1, len(cols), n_wave, cell_w=0.55, cell_h=1.5)
    for c, clab in enumerate(cols):
        mat = mean_matrix(sub, key) if clab == "mean over datasets" \
            else metric_matrix(sub, clab, key)
        ax = fig.add_subplot(gs[0, 3 * c])
        im = _draw_matrix(ax, slice_matrix(mat, wave_idx), fmt, cmap, vmin, vmax,
                          show_ylabel=(c == 0), wave=wave_idx)
        ax.set_title(clab, fontsize=11)
        cax = fig.add_subplot(gs[0, 3 * c + 1])
        fig.colorbar(im, cax=cax)
    fig.suptitle(f"{title}  ($N_{{th}}$ = {nth})", fontsize=12)
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")


def wave_profiles(rows, key, title, filename, ylabel, vmin=None, vmax=None):
    """Metric vs Eth_wave, one line per Nth, at the best seed threshold.

    Three panels: one per dataset plus the dataset average.  This is the plot
    that exposes the *physical* onset of over-splitting as Eth_wave rises.
    """
    nths = available_nths(rows)
    panels = DATASETS + ["mean over datasets"]

    # pick, per Nth, the seed with the best mean-over-datasets value at the
    # loosest wave threshold so every curve is anchored at its own optimum
    wave0 = min(ETH_WAVE_VALUES)
    best_seed = {}
    for nth in nths:
        best, best_val = None, -np.inf
        for s in ETH_SEED_VALUES:
            vals = []
            for ds in DATASETS:
                for r in rows:
                    if (int(r["Nth"]) == nth and r["dataset"] == ds
                            and r["Eth_seed"] == s and r["Eth_wave"] == wave0):
                        v = r.get(key, np.nan)
                        if np.isfinite(v):
                            vals.append(v)
            if vals and np.mean(vals) > best_val:
                best, best_val = s, float(np.mean(vals))
        best_seed[nth] = best

    fig, axes = plt.subplots(1, len(panels), figsize=(3.4 * len(panels), 3.4),
                             sharey=True, squeeze=False)
    colors = plt.cm.viridis(np.linspace(0, 0.85, len(nths)))
    for ax, ds in zip(axes[0], panels):
        for color, nth in zip(colors, nths):
            s = best_seed[nth]
            if s is None:
                continue
            xs, ys = [], []
            for w in ETH_WAVE_VALUES:
                if w >= s:
                    continue
                if ds == "mean over datasets":
                    vals = [r.get(key, np.nan) for r in rows
                            if int(r["Nth"]) == nth and r["Eth_seed"] == s
                            and r["Eth_wave"] == w]
                    vals = [v for v in vals if np.isfinite(v)]
                    if vals:
                        xs.append(w)
                        ys.append(float(np.mean(vals)))
                else:
                    for r in rows:
                        if (int(r["Nth"]) == nth and r["dataset"] == ds
                                and r["Eth_seed"] == s and r["Eth_wave"] == w):
                            xs.append(w)
                            ys.append(r.get(key, np.nan))
            if xs:
                ax.plot(xs, ys, "o-", color=color, lw=1.8, ms=4,
                        label=f"$N_{{th}}$={nth} (seed={s:g})")
        if vmin is not None or vmax is not None:
            ax.set_ylim(vmin, vmax)
        ax.set_xscale("log")
        ax.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.tick_params(axis="x", labelsize=8, labelrotation=45)
        for t in ax.get_xticklabels():
            t.set_ha("right")
        ax.set_xlabel(r"$E_{th}^{\rm wave}$ (MeV, log scale)")
        ax.set_title(ds)
        ax.grid(alpha=0.3)
    axes[0][0].set_ylabel(ylabel)
    axes[0][0].legend(fontsize=9)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def ntrack_histograms(rows, metric="exact_match_rate", top=6):
    """Bar chart of the per-Nth best points for the headline metric."""
    nths = available_nths(rows)
    fig, ax = plt.subplots(figsize=(1.5 * len(DATASETS) * len(nths) + 2.0, 3.8))
    width = 0.8 / len(nths)
    x = np.arange(len(DATASETS))
    for k, nth in enumerate(nths):
        vals, labels = [], []
        for ds in DATASETS:
            sub = [r for r in rows if int(r["Nth"]) == nth and r["dataset"] == ds]
            if not sub:
                vals.append(np.nan)
                labels.append("")
                continue
            pick = max(sub, key=lambda r: r.get(metric, -np.inf))
            vals.append(pick[metric])
            labels.append(f"{pick['Eth_seed']:g}/{pick['Eth_wave']:g}")
        bars = ax.bar(x + k * width - 0.4 + width / 2, vals, width,
                      label=f"$N_{{th}}$={nth}")
        for b, lab in zip(bars, labels):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.005,
                    lab, ha="center", fontsize=7, rotation=90)
    ax.set_xticks(x, DATASETS)
    ax.set_ylabel(metric)
    ax.set_title(f"Best {metric} per dataset / per $N_{{th}}$\n"
                 "(labels: best $E_{th}^{\\rm seed}$ / $E_{th}^{\\rm wave}$)")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = RESULTS_DIR / f"best_{metric}_bars.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


def tradeoff_scatter(rows, filename="tradeoff_scatter.png"):
    """Multiplicity quality vs energy capture, one point per sweep point."""
    panels = DATASETS + ["mean over datasets"]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.6 * len(panels), 3.6),
                             sharex=True, sharey=True, squeeze=False)
    sc = None
    for ax, ds in zip(axes[0], panels):
        sel = rows if ds == "mean over datasets" else \
            [r for r in rows if r["dataset"] == ds]
        x = [r.get("energy_captured_fraction", np.nan) for r in sel]
        y = [r.get("exact_match_rate", np.nan) for r in sel]
        c = [r["Eth_seed"] for r in sel]
        m = [50 * (0.5 + int(r["Nth"])) for r in sel]
        sc = ax.scatter(x, y, c=c, s=m, cmap="viridis", edgecolors="k",
                        linewidths=0.3, alpha=0.8)
        ax.set_title(ds)
        ax.set_xlabel("energy capture (all clusters / event)")
        ax.grid(alpha=0.3)
    axes[0][0].set_ylabel("exact-match rate")
    fig.colorbar(sc, ax=axes[0][-1], label=r"$E_{th}^{\rm seed}$ (MeV)")
    fig.suptitle("Multiplicity quality vs energy capture "
                 "(marker size ~ $N_{th}$; one point per sweep point)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS_DIR / filename
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"saved {out}")


# --------------------------------------------------------------------------- #
# tables
# --------------------------------------------------------------------------- #
BEST_METRICS = [
    ("exact_match_rate", "max"),
    ("energy_captured_fraction", "max"),
    # absolute energy that no cluster absorbed - larger Eth_wave can only make
    # this worse, unlike the threshold-relative lit_energy_missed_fraction
    ("energy_unclustered_fraction", "min"),
    ("energy_containment_mean", "max"),
    ("frac_events_containment_ge_90", "max"),
    ("mean_clusters", "min"),
    ("MAE", "min"),
    ("RMSE", "min"),
]


def best_tables(rows):
    nths = available_nths(rows)
    lines = []

    lines.append("## Best point per metric and per neighbour threshold\n")
    lines.append("| metric | goal | Nth | dataset | seed | wave | value |")
    lines.append("|:--|:--|--:|:--|--:|--:|--:|")
    for key, mode in BEST_METRICS:
        for nth in nths:
            for ds in DATASETS + [None]:
                if ds is not None:
                    sub = [r for r in rows
                           if int(r["Nth"]) == nth and r["dataset"] == ds
                           and np.isfinite(r.get(key, np.nan))]
                    if not sub:
                        continue
                    pick = (max if mode == "max" else min)(sub, key=lambda r: r[key])
                    s_val, w_val, v_val = (pick["Eth_seed"], pick["Eth_wave"],
                                           pick[key])
                else:
                    # genuine dataset average at each (seed, wave)
                    groups = {}
                    for r in rows:
                        if int(r["Nth"]) != nth or not np.isfinite(r.get(key, np.nan)):
                            continue
                        groups.setdefault((r["Eth_seed"], r["Eth_wave"]), []).append(
                            r[key])
                    if not groups:
                        continue
                    avg = {k: float(np.mean(v)) for k, v in groups.items()}
                    (s_val, w_val), v_val = (max if mode == "max" else min)(
                        avg.items(), key=lambda kv: kv[1])
                lines.append(
                    f"| `{key}` | {mode} | {nth} | {ds or '**mean**'} "
                    f"| {s_val:g} | {w_val:g} | {v_val:.4f} |")

    lines.append("\n## Effect of the neighbour threshold on the optimum\n")
    lines.append("| Nth | best mean exact-match | seed | wave | energy capture "
                 "at that point | unclustered energy (abs.) | missed lit energy "
                 "(threshold-relative) |")
    lines.append("|--:|--:|--:|--:|--:|--:|--:|")
    for nth in nths:
        sub = [r for r in rows if int(r["Nth"]) == nth]
        if not sub:
            continue
        # rank by the dataset-averaged exact-match at each (seed, wave)
        grid = {}
        for r in sub:
            grid.setdefault((r["Eth_seed"], r["Eth_wave"]), []).append(
                r.get("exact_match_rate", np.nan))
        grid = {k: np.nanmean(v) for k, v in grid.items() if np.isfinite(v).any()}
        if not grid:
            continue
        (s, w), v = max(grid.items(), key=lambda kv: kv[1])
        at = [r for r in sub if r["Eth_seed"] == s and r["Eth_wave"] == w]

        def avg(k):
            return np.nanmean([r.get(k, np.nan) for r in at])

        lines.append(f"| {nth} | {v:.4f} | {s:g} | {w:g} "
                     f"| {avg('energy_captured_fraction'):.4f} "
                     f"| {avg('energy_unclustered_fraction'):.4f} "
                     f"| {avg('lit_energy_missed_fraction'):.4f} |")

    lines.append("\n## Top-10 sweep points by mean exact-match rate\n")
    lines.append("| rank | Nth | seed | wave | mean exact-match "
                 "| mean energy capture | mean unclustered energy |")
    lines.append("|--:|--:|--:|--:|--:|--:|--:|")
    grid = {}
    for r in rows:
        grid.setdefault((int(r["Nth"]), r["Eth_seed"], r["Eth_wave"]), []).append(r)
    ranked = sorted(grid.items(),
                    key=lambda kv: np.nanmean([r.get("exact_match_rate", np.nan)
                                               for r in kv[1]]),
                    reverse=True)[:10]
    for rank, ((nth, s, w), rs) in enumerate(ranked, 1):
        lines.append(
            f"| {rank} | {nth} | {s:g} | {w:g} "
            f"| {np.nanmean([r.get('exact_match_rate', np.nan) for r in rs]):.4f} "
            f"| {np.nanmean([r.get('energy_captured_fraction', np.nan) for r in rs]):.4f} "
            f"| {np.nanmean([r.get('energy_unclustered_fraction', np.nan) for r in rs]):.4f} |")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
def main():
    global RESULTS_DIR, SUMMARY_CSV

    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=SUMMARY_CSV)
    ap.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    args = ap.parse_args()

    RESULTS_DIR = args.results_dir
    SUMMARY_CSV = args.csv
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    rows = load_rows(SUMMARY_CSV)
    if not rows:
        print(f"no rows found in {SUMMARY_CSV} - run run_sweep.py first")
        return

    missing = [k for k in REQUIRED_METRICS if k not in rows[0]]
    if missing:
        print(f"CSV lacks required columns {missing}; nothing to plot")
        return

    print(f"analysing {len(rows)} sweep points from {SUMMARY_CSV}")
    print(f"Nth values present: {available_nths(rows)}")

    # zoom window over the sampled wave thresholds where the behaviour turns
    # over; keeps the per-cell digits readable with this many grid columns
    zoom_idx = [i for i, w in enumerate(ETH_WAVE_VALUES) if w <= 1.5]

    for key, title, cmap, vmin, vmax, fmt in HEATMAP_METRICS:
        if key not in rows[0]:
            print(f"  skip {key} (not in CSV)")
            continue
        heatmap_per_nth(rows, key, title, f"heatmap_{key}_per_nth.png",
                        cmap=cmap, fmt=fmt, vmin=vmin, vmax=vmax,
                        wave_idx=zoom_idx, cell_w=0.55, cell_h=1.6, val_fs=6.2,
                        tick_fs=6.5)
        heatmap_per_nth(rows, key, title + "  [full wave range]",
                        f"heatmap_{key}_per_nth_full.png",
                        cmap=cmap, fmt="{:.2f}", vmin=vmin, vmax=vmax,
                        wave_idx=None, annotate=False, cell_w=0.24, cell_h=1.5,
                        label_stride=2)

    # dedicated low-threshold zoom: the refined campaign adds wave points below
    # 0.1 MeV, so this view answers "where does the plateau actually start?"
    low_idx = [i for i, w in enumerate(ETH_WAVE_VALUES) if w <= 0.5]
    if len(low_idx) >= 3:
        for key, title, cmap, vmin, vmax, fmt in HEATMAP_METRICS:
            if key not in ("exact_match_rate", "energy_captured_fraction"):
                continue
            heatmap_per_nth(rows, key, title + "  [low-threshold zoom]",
                            f"heatmap_{key}_per_nth_low.png",
                            cmap=cmap, fmt=fmt, vmin=vmin, vmax=vmax,
                            wave_idx=low_idx, cell_w=0.80, cell_h=1.6,
                            val_fs=6.5)

    # compact single-Nth headline figures for the recommended working points
    for nth in available_nths(rows):
        heatmap_single_nth(rows, "exact_match_rate", nth,
                           "Exact-match rate",
                           f"heatmap_exact_match_nth{nth}.png",
                           cmap="RdYlGn", vmin=0.0, vmax=1.0, wave_idx=zoom_idx)

    wave_profiles(rows, "exact_match_rate",
                  "Exact-match rate vs wave threshold (best seed per Nth)",
                  "profile_exact_match.png", "exact-match rate", vmin=0.0, vmax=1.0)
    if "energy_captured_fraction" in rows[0]:
        wave_profiles(rows, "energy_captured_fraction",
                      "Energy capture vs wave threshold (best seed per Nth)",
                      "profile_energy_capture.png", "energy capture", vmin=0.0, vmax=1.0)
    if "lit_energy_missed_fraction" in rows[0]:
        wave_profiles(rows, "lit_energy_missed_fraction",
                      "Missed lit energy vs wave threshold (best seed per Nth)",
                      "profile_missed_energy.png", "missed lit energy fraction",
                      vmin=0.0, vmax=0.5)

    ntrack_histograms(rows, "exact_match_rate")
    tradeoff_scatter(rows)

    text = ("# Best parameters - refined two-threshold sweep\n\n"
            + best_tables(rows) + "\n")
    (RESULTS_DIR / "best_params.md").write_text(text)
    print(f"\nwrote {RESULTS_DIR / 'best_params.md'}")
    print("\n" + best_tables(rows))


if __name__ == "__main__":
    main()
