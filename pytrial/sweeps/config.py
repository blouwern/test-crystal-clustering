import os
from pathlib import Path

DATA_DIR = os.environ.get("ECAL_CLUSTERING_DATA_DIR", ".")
SWEEP_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SWEEP_DIR / "results"
ANALYSIS_OUTPUT_DIR = SWEEP_DIR / "analysis_output"
DATA_PATH = Path(DATA_DIR) / "processed" / "child_stitched_edeps.root"

ETH_VALUES = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 20.0]
NTH_VALUES = [1.0, 2.0, 3.0]

SUMMARY_CSV = ANALYSIS_OUTPUT_DIR / "sweep_summary.csv"
