import ROOT
import random

file = ROOT.TFile.Open("data/raw/SimMACEPhaseI_20260401.root")
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

def print_entry(entry):
        print(f"[NO.{entry.ModID}:{entry.Edep}MeV]")

for times in range(10):
    rand_evtid = random.randint(0, N_evt - 1)
    print(f"---Selecting event NO.{rand_evtid}...")

    secpdgids = vertex_dict.get(rand_evtid, [])
    sec_str = ", ".join(str(pdg) for pdg in secpdgids)
    print(f"Secondary PDGID list: {sec_str}")
    for entry in hit_tree:
        if int(entry.EvtID) == rand_evtid:
            print_entry(entry)

    print("\n---End of this event\n")

file.Close()
