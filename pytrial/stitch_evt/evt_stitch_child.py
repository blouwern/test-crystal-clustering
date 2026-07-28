import numpy as np
import ROOT
from ROOT import vector
from pathlib import Path
import array


data_dir = Path("$ECAL_CLUSTERING_DATA_DIR")
filename_edep_of_each_evt = Path("edep_of_each_evt_20260401.root")
expected_n_track = 3
n_samples = 2000

n_track = np.random.poisson(lam=expected_n_track, size=n_samples)

input_path = data_dir / "processed" / filename_edep_of_each_evt
input_file = ROOT.TFile.Open(str(input_path), "READ")
edep_tree = input_file.Get("EdepOfEachEvt")

edeps_list = []
edeps_vec = vector("float")()
edep_tree.SetBranchAddress("Edeps", edeps_vec)

n_entry = edep_tree.GetEntries()
for i in range(n_entry):
    edep_tree.GetEntry(i)
    arr = np.array(edeps_vec, dtype=np.float32)
    edeps_list.append(arr)

input_file.Close()

lengths = [len(arr) for arr in edeps_list]
if len(set(lengths)) > 1:
    raise ValueError("Input Edeps lengths are inconsistent. Align them first.")
vec_len = lengths[0] if lengths else 0


stitched_arrays = []
for i in range(n_samples):
    k = n_track[i]
    if k == 0:
        stitched_arrays.append(np.zeros(vec_len, dtype=np.float32))
        continue
    if k > n_entry:
        raise ValueError(f"n_track={k} exceeds given {n_entry} entries.")
    indices = np.random.choice(n_entry, size=k, replace=False)
    sum_arr = np.sum([edeps_list[idx] for idx in indices], axis=0)
    stitched_arrays.append(sum_arr.astype(np.float32))


output_path = data_dir / "processed" / "child_stitched_edeps.root"
output_file = ROOT.TFile.Open(str(output_path), "RECREATE")
output_tree = ROOT.TTree("EdepOfEachEvt", "Stiched Edep events (childish)")

out_edeps = vector("float")()
n_track_val = array.array("i", [0])
output_tree.Branch("Edeps", out_edeps)
output_tree.Branch("NTrack", n_track_val, "NTrack/I")

for i, arr in enumerate(stitched_arrays):
    out_edeps.clear()
    for val in arr:
        out_edeps.push_back(float(val))
    n_track_val[0] = n_track[i]
    output_tree.Fill()

output_tree.Write()
output_file.Close()
print(f"{n_samples} childish events Edeps stiched，saved to {output_path}")
