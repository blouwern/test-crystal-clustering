import sys
import csv
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import reconstruction
import scoring
from sweeps.config import (
    DATA_PATH, ETH_VALUES, NTH_VALUES,
    RESULTS_DIR, ANALYSIS_OUTPUT_DIR, SUMMARY_CSV,
)


def run_sweep(filename, event_by_event=False):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    ANALYSIS_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []

    for eth in ETH_VALUES:
        for nth in NTH_VALUES:
            out_name = f"ecal_clusters_Eth-{eth}_Nth-{nth}.root"
            out_path = RESULTS_DIR / out_name

            label = f"Eth={eth} Nth={nth}"
            print(f"\n{'=' * 60}")
            print(f"  {label}")
            print(f"{'=' * 60}")

            if out_path.exists():
                print(f"  Output exists, skipping reconstruction ({out_name})")
            else:
                reconstruction.reconstruction(
                    filename, eth, nth,
                    output_path=str(out_path),
                )

            rec_file = out_path

            n_clusters = scoring._read_cluster_counts(rec_file)
            N = len(n_clusters)

            mean_clusters = float(np.mean(n_clusters))
            std_clusters = float(np.std(n_clusters, ddof=1))

            n_track_arr = scoring._read_ntrack(str(DATA_PATH), N) if DATA_PATH.exists() else None

            row = {
                "Eth": eth,
                "Nth": nth,
                "N_events": N,
                "mean_clusters": mean_clusters,
                "std_clusters": std_clusters,
            }

            if n_track_arr is not None:
                diff = n_clusters - n_track_arr
                abs_diff = np.abs(diff)
                row["exact_match_rate"] = float(np.mean(diff == 0))
                row["over_cluster_rate"] = float(np.mean(diff > 0))
                row["under_cluster_rate"] = float(np.mean(diff < 0))
                row["MAE"] = float(np.mean(abs_diff))
                row["RMSE"] = float(np.sqrt(np.mean(diff.astype(np.float64) ** 2)))
                row["pearson_r"] = float(np.corrcoef(n_clusters, n_track_arr)[0, 1])
                row["mean_truth"] = float(np.mean(n_track_arr))

            rows.append(row)

            print(f"  Mean clusters: {mean_clusters:.4f}  "
                  f"Exact match: {row.get('exact_match_rate', 'N/A')}")

    save_summary(rows)
    print(f"\nSummary saved to {SUMMARY_CSV}")
    return rows


def save_summary(rows):
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(SUMMARY_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename", help="Input edeps filename (under processed/)")
    parser.add_argument("--event-by-event", action="store_true")
    args = parser.parse_args()
    run_sweep(args.filename, event_by_event=args.event_by_event)
