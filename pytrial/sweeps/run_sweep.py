"""Drive the (Eth_seed, Eth_wave) x dataset sweep and write the CSV summary.

    python sweeps/run_sweep.py                 # cache-aware
    python sweeps/run_sweep.py --force         # recompute everything
    python sweeps/run_sweep.py --csv-only      # rebuild CSV from cached JSONs
"""

import argparse
import csv
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sweeps.config import (  # noqa: E402
    DATASETS, N_EVENTS, NTH_VALUES, N_WORKERS, NEIGHBOR_FILE,
    RESULTS_DIR, SUMMARY_CSV, TASK_DIR, valid_pairs,
)

TASK_SCRIPT = Path(__file__).resolve().parent / "run_task.py"


def run_one(task, force):
    dataset, seed, wave, nth = task
    tag = f"{dataset}_seed{seed:g}_wave{wave:g}_nth{nth:g}"
    if (TASK_DIR / f"{tag}.json").exists() and not force:
        return tag, "cached"
    proc = subprocess.run(
        [sys.executable, str(TASK_SCRIPT), dataset, str(seed), str(wave), str(nth)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        return tag, f"FAILED rc={proc.returncode}: {proc.stderr.strip()[-300:]}"
    return tag, "ok"


def collect_rows():
    rows = []
    for path in sorted(TASK_DIR.glob("*.json")):
        with open(path) as f:
            rows.append(json.load(f))
    rows.sort(key=lambda r: (r["dataset"], r["Eth_seed"], r["Eth_wave"], r["Nth"]))
    return rows


def write_csv(rows):
    if not rows:
        print("no task results found")
        return
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())
    with open(SUMMARY_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows -> {SUMMARY_CSV}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--workers", type=int, default=N_WORKERS)
    ap.add_argument("--csv-only", action="store_true")
    args = ap.parse_args()

    TASK_DIR.mkdir(parents=True, exist_ok=True)
    pairs = valid_pairs()
    tasks = [(ds, s, w, n) for ds in DATASETS for (s, w) in pairs for n in NTH_VALUES]

    if not args.csv_only:
        print(f"sweep: {len(pairs)} (seed,wave) pairs x {len(DATASETS)} datasets "
              f"x {len(NTH_VALUES)} Nth = {len(tasks)} points | Nth={NTH_VALUES} "
              f"| topology={NEIGHBOR_FILE} | {N_EVENTS} events/dataset")
        t0 = time.time()
        failures = []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_one, t, args.force): t for t in tasks}
            for i, fut in enumerate(as_completed(futures), 1):
                tag, status = fut.result()
                if status.startswith("FAILED"):
                    failures.append((tag, status))
                print(f"[{i:3d}/{len(tasks)}] {tag:34s} {status}", flush=True)
        print(f"finished in {time.time() - t0:.1f}s, {len(failures)} failures")
        for tag, status in failures:
            print(f"  {tag}: {status}")

    write_csv(collect_rows())


if __name__ == "__main__":
    main()
