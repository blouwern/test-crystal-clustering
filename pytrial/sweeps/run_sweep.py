"""Drive the Eth x Nth x dataset sweep and write the CSV summary.

Each grid point is executed in its own subprocess (see ``run_task.py``) so that
ROOT state never leaks between points and the sweep can use all CPU cores.

    python sweeps/run_sweep.py                # run everything (cache-aware)
    python sweeps/run_sweep.py --force        # recompute every point
    python sweeps/run_sweep.py --workers 8
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
    DATASETS, ETH_VALUES, NTH_VALUES, N_EVENTS, N_WORKERS,
    RESULTS_DIR, SUMMARY_CSV, TASK_DIR,
)

TASK_SCRIPT = Path(__file__).resolve().parent / "run_task.py"


def run_one(task, force):
    dataset, eth, nth = task
    tag = f"{dataset}_Eth{eth:g}_Nth{nth:g}"
    json_path = TASK_DIR / f"{tag}.json"
    if json_path.exists() and not force:
        return tag, "cached"
    proc = subprocess.run(
        [sys.executable, str(TASK_SCRIPT), dataset, str(eth), str(nth)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        return tag, f"FAILED rc={proc.returncode}: {proc.stderr.strip()[-400:]}"
    return tag, "ok"


def collect_rows():
    rows = []
    for path in sorted(TASK_DIR.glob("*.json")):
        with open(path) as f:
            rows.append(json.load(f))
    rows.sort(key=lambda r: (r["dataset"], r["Eth"], r["Nth"]))
    return rows


def write_csv(rows):
    if not rows:
        print("no task results found")
        return
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(SUMMARY_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {len(rows)} rows -> {SUMMARY_CSV}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true",
                        help="recompute points even if a cached JSON exists")
    parser.add_argument("--workers", type=int, default=N_WORKERS)
    parser.add_argument("--csv-only", action="store_true",
                        help="only rebuild the CSV from existing JSON records")
    args = parser.parse_args()

    TASK_DIR.mkdir(parents=True, exist_ok=True)

    tasks = [(ds, eth, nth)
             for ds in DATASETS
             for eth in ETH_VALUES
             for nth in NTH_VALUES]

    if not args.csv_only:
        print(f"sweep: {len(tasks)} points "
              f"({len(DATASETS)} datasets x {len(ETH_VALUES)} Eth x {len(NTH_VALUES)} Nth), "
              f"{N_EVENTS} events/dataset, {args.workers} workers")
        t0 = time.time()
        failures = []
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(run_one, t, args.force): t for t in tasks}
            for i, future in enumerate(as_completed(futures), 1):
                tag, status = future.result()
                if status.startswith("FAILED"):
                    failures.append((tag, status))
                print(f"[{i:3d}/{len(tasks)}] {tag:30s} {status}", flush=True)
        print(f"sweep finished in {time.time() - t0:.1f}s, {len(failures)} failures")
        for tag, status in failures:
            print(f"  {tag}: {status}")

    rows = collect_rows()
    write_csv(rows)
    if rows:
        print(f"datasets={sorted({r['dataset'] for r in rows})} "
              f"points/dataset={len(rows) // len({r['dataset'] for r in rows})}")


if __name__ == "__main__":
    main()
