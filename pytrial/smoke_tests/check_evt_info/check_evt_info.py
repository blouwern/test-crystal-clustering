import ROOT
import random

ROOT.EnableImplicitMT()
file = ROOT.TFile.Open("$ECAL_CLUSTERING_DATA_DIR/raw/SimMACEPhaseI_20260401.root")
if not file or file.IsZombie():
    print("Error: Cannot open file!")
    exit(1)

hit_tree = file.Get("G4Run0/ECALSimHit")
vertex_tree = file.Get("G4Run0/SimDecayVertex")

if not hit_tree or not vertex_tree:
    print("Error: TTree not found!")
    exit(1)

max_evtid = int(hit_tree.GetMaximum("EvtID"))
N_evt = max_evtid + 1

vertex_dict = {}
for entry in vertex_tree:
    evtid = int(entry.EvtID)
    secpdgid = list(entry.SecPDGID)
    if evtid not in vertex_dict:
        vertex_dict[evtid] = []
        vertex_dict[evtid].append(secpdgid)

evt_modID_list = []


def print_entry(entry):
    if evt_modID_list.count(entry.ModID) > 0:
        print(f"[ABNORMAL]!!!--NO.{entry.ModID} existence greater than 1")
    print(f"[NO.{entry.ModID}:{entry.Edep}MeV]")
    evt_modID_list.append(entry.ModID)


for evtID in range(0, 5):
    secpdgid_evt = vertex_dict.get(evtID, [])

    print(f"---Selecting event NO.{evtID}...")

    sec_str = ", ".join(str(pdg) for pdg in secpdgid_evt)
    print(f"Secondary PDGID list: {sec_str}")
    for entry in hit_tree:
        if int(entry.EvtID) == evtID:
            print_entry(entry)
    evt_modID_list = []
    print("---End of this event\n")

file.Close()
