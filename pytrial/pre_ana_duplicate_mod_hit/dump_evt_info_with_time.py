import ROOT
from collections import Counter

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


def print_evt_info_with_time(dict):
    print(f"{'ModID':^5} {'Edep':^20} {'t':^20}")
    for row in range(0, len(dict["ModID"])):
        print(f"{dict['ModID'][row]:^5} {dict['Edep'][row]:^20} {dict['t'][row]:^20}")


def check_duplicates_and_print(info_dict, evtid):
    modid_list = info_dict["ModID"]
    time_list = info_dict["t"]

    counter = Counter(modid_list)
    duplicated_mids = [mid for mid, cnt in counter.items() if cnt > 1]

    if not duplicated_mids:
        return

    print(f"==== Event NO.{evtid} ====")
    for mid in duplicated_mids:

        times = [time_list[i] for i, m in enumerate(modid_list) if m == mid]

        times.sort()

        print(f"  >> ModID {mid} appears {len(times)} times.")
        print(f"     Hit times (ns): {times}")

        if len(times) > 1:

            consecutive_diffs = [times[i + 1] - times[i] for i in range(len(times) - 1)]
            print(f"     Consecutive diffs (ns): {consecutive_diffs}")

            first_t = times[0]
            diffs_from_first = [t - first_t for t in times[1:]]
            print(f"     Diffs from 1st hit (ns): {diffs_from_first}")


for evtid, info_dict in enumerate(info_dict_list):
    check_duplicates_and_print(info_dict, evtid)

file_read.Close()
