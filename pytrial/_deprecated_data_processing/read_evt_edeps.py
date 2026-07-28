import ROOT
import get_N_mod
import filter_timeout_hit

ROOT.EnableImplicitMT()
file_read = ROOT.TFile.Open("$ECAL_CLUSTERING_DATA_DIR/raw/SimMACEPhaseI_20260401.root")
if not file_read or file_read.IsZombie():
    print("Error: Cannot open file_read!")
    exit(1)

hit_tree = file_read.Get("G4Run0/ECALSimHit")
if not hit_tree:
    print("Error: TTree not found!")
    exit(1)

max_evtid = int(hit_tree.GetMaximum("EvtID"))
N_evt = max_evtid + 1

info_dict_list = [{"ModID": [], "Edep": [], "t": []} for _ in range(N_evt)]


def fill_entry_to_dict(entry, dict):
    dict["ModID"].append(entry.ModID)
    dict["Edep"].append(entry.Edep)
    dict["t"].append(entry.t)


for entry in hit_tree:
    fill_entry_to_dict(entry, info_dict_list[entry.EvtID])

for idx, dict in enumerate(info_dict_list):
    info_dict_list[idx] = filter_timeout_hit.filter_timeout_hit(dict)

file_output = ROOT.TFile.Open("edeps_of_each_evt_20260401.root", "RECREATE")
if not file_output or file_output.IsZombie():
    print("Error: Cannot open file_output!")
    exit(1)

output_tree = ROOT.TTree("edep_of_each_evt", "Edep vector of each event")
if not output_tree:
    print("Error: TTree not created!")
    exit(1)

vec_init = ROOT.std.vector("float")(get_N_mod.get_N_mod(), 0.0)
vec_edeps_of_evt = vec_init
output_tree.Branch("Edeps", vec_edeps_of_evt)


def fill_dict_to_vec(dict, vec):
    for modid, edep, t in zip(dict["ModID"], dict["Edep"], dict["t"]):
        vec[modid] = edep


for dict in info_dict_list:
    fill_dict_to_vec(dict, vec_edeps_of_evt)
    output_tree.Fill()

file_output.Write()
file_output.Close()
file_read.Close()
