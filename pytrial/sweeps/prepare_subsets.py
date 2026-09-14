"""Extract the first ``N_EVENTS`` events of each processed dataset.

The reconstruction reads ``<data_dir>/processed/<file>`` and
``<data_dir>/utilities/ecal_neighbor_info_added.root``, so the sweep builds a
small self-contained data directory under ``sweeps/subsets`` that holds the
event subsets plus symlinks to the topology file.
"""

import sys
from pathlib import Path

import ROOT
from ROOT import vector

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sweeps.config import (  # noqa: E402
    DATASETS, NEIGHBOR_FILE, N_EVENTS, PROCESSED_DIR, SUBSET_DIR, UTILITIES_DIR,
)


def _link_utilities():
    (SUBSET_DIR / "utilities").mkdir(parents=True, exist_ok=True)
    for name in (NEIGHBOR_FILE, "ecal_neighbor_info.root"):
        src = UTILITIES_DIR / name
        dst = SUBSET_DIR / "utilities" / name
        if src.exists() and not dst.exists():
            dst.symlink_to(src.resolve())


def make_subset(dataset, n_events):
    src = PROCESSED_DIR / f"edep_of_each_evt_{dataset}.root"
    if not src.exists():
        raise FileNotFoundError(f"missing processed dataset: {src}")

    dst = SUBSET_DIR / "processed" / f"subset_{dataset}.root"
    dst.parent.mkdir(parents=True, exist_ok=True)

    fin = ROOT.TFile.Open(str(src), "READ")
    tree = fin.Get("EdepOfEachEvt")
    total = tree.GetEntries()
    n = min(n_events, total)

    edeps = vector("float")()
    tree.SetBranchAddress("Edeps", edeps)

    fout = ROOT.TFile.Open(str(dst), "RECREATE")
    out_tree = ROOT.TTree("EdepOfEachEvt", f"{dataset} subset ({n} events)")
    out_edeps = vector("float")()
    out_tree.Branch("Edeps", out_edeps)

    for i in range(n):
        tree.GetEntry(i)
        out_edeps.clear()
        for value in edeps:
            out_edeps.push_back(float(value))
        out_tree.Fill()

    out_tree.Write()
    fout.Close()
    fin.Close()
    print(f"  {dataset:10s} {n:>7d}/{total} events -> {dst}")
    return n


def main():
    print(f"Building subsets with N_EVENTS={N_EVENTS}")
    _link_utilities()
    for dataset in DATASETS:
        make_subset(dataset, N_EVENTS)
    print("done")


if __name__ == "__main__":
    main()
