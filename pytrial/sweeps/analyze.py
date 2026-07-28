import csv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from sweeps.config import SUMMARY_CSV, ANALYSIS_OUTPUT_DIR, ETH_VALUES, NTH_VALUES


def load_summary():
    rows = []
    with open(SUMMARY_CSV, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for k in row:
                try:
                    row[k] = float(row[k])
                except (ValueError, TypeError):
                    pass
            rows.append(row)
    return rows


def build_matrix(rows, key):
    n_eth = len(ETH_VALUES)
    n_nth = len(NTH_VALUES)
    mat = np.full((n_eth, n_nth), np.nan)
    for i, eth in enumerate(ETH_VALUES):
        for j, nth in enumerate(NTH_VALUES):
            for r in rows:
                if r["Eth"] == eth and r["Nth"] == nth:
                    mat[i, j] = r.get(key, np.nan)
                    break
    return mat


def print_table(matrix, row_labels, col_labels, title, fmt=".4f"):
    print(f"\n{title}")
    header = f"{'':>8s}  " + "  ".join(f"{c:>10s}" for c in col_labels)
    print(header)
    print("-" * len(header))
    for i, label in enumerate(row_labels):
        vals = "  ".join(
            f"{matrix[i, j]:10{fmt}}" if not np.isnan(matrix[i, j]) else f"{'N/A':>10s}"
            for j in range(len(col_labels))
        )
        print(f"  {label:<6s}  {vals}")


def analyze():
    ANALYSIS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = load_summary()

    if not rows:
        print("No sweep data found. Run run_sweep.py first.")
        return

    eth_labels = [f"E={int(e)}" for e in ETH_VALUES]
    nth_labels = [f"N={int(n)}" for n in NTH_VALUES]

    metrics = [
        ("exact_match_rate", "Exact match rate (predicted == truth)"),
        ("MAE", "Mean Absolute Error"),
        ("RMSE", "Root Mean Squared Error"),
        ("pearson_r", "Pearson correlation"),
        ("mean_clusters", "Mean predicted clusters"),
    ]

    for key, title in metrics:
        if key not in rows[0]:
            continue
        mat = build_matrix(rows, key)
        fmt = ".4f"
        if key == "mean_clusters":
            fmt = ".2f"
        print_table(mat, eth_labels, nth_labels, title, fmt=fmt)

    if "MAE" in rows[0] and "exact_match_rate" in rows[0]:
        best_mae = min(rows, key=lambda r: r["MAE"])
        best_match = max(rows, key=lambda r: r["exact_match_rate"])
        print(f"\nBest MAE:     Eth={best_mae['Eth']}, Nth={best_mae['Nth']}, "
              f"MAE={best_mae['MAE']:.4f}")
        print(f"Best match:   Eth={best_match['Eth']}, Nth={best_match['Nth']}, "
              f"match_rate={best_match['exact_match_rate']:.4f}")

    if "exact_match_rate" in rows[0]:
        mat = build_matrix(rows, "exact_match_rate")
        fig, ax = plt.subplots(figsize=(7, 5))
        im = ax.imshow(mat, cmap="RdYlGn", aspect="auto", vmin=0, vmax=1)
        ax.set_xticks(range(len(NTH_VALUES)))
        ax.set_xticklabels(nth_labels)
        ax.set_yticks(range(len(ETH_VALUES)))
        ax.set_yticklabels(eth_labels)
        ax.set_xlabel("N_th")
        ax.set_ylabel("E_th")
        ax.set_title("Exact match rate (predicted == ground-truth)")
        for i in range(len(ETH_VALUES)):
            for j in range(len(NTH_VALUES)):
                val = mat[i, j]
                if not np.isnan(val):
                    ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                            fontsize=9, color="black" if 0.3 < val < 0.7 else "white")
        plt.colorbar(im, ax=ax)
        plt.tight_layout()
        out_path = ANALYSIS_OUTPUT_DIR / "heatmap_exact_match.png"
        plt.savefig(out_path, dpi=150)
        print(f"\nHeatmap saved to {out_path}")

    if "MAE" in rows[0]:
        mat = build_matrix(rows, "MAE")
        fig, ax = plt.subplots(figsize=(7, 5))
        im = ax.imshow(mat, cmap="RdYlGn_r", aspect="auto")
        ax.set_xticks(range(len(NTH_VALUES)))
        ax.set_xticklabels(nth_labels)
        ax.set_yticks(range(len(ETH_VALUES)))
        ax.set_yticklabels(eth_labels)
        ax.set_xlabel("N_th")
        ax.set_ylabel("E_th")
        ax.set_title("Mean Absolute Error (predicted vs ground-truth)")
        for i in range(len(ETH_VALUES)):
            for j in range(len(NTH_VALUES)):
                val = mat[i, j]
                if not np.isnan(val):
                    ax.text(j, i, f"{val:.3f}", ha="center", va="center", fontsize=9)
        plt.colorbar(im, ax=ax)
        plt.tight_layout()
        out_path = ANALYSIS_OUTPUT_DIR / "heatmap_mae.png"
        plt.savefig(out_path, dpi=150)
        print(f"Heatmap saved to {out_path}")


if __name__ == "__main__":
    analyze()
