import ROOT
from ROOT import std, vector
import sys
from pathlib import Path
import array


def reconstruction(
    filename_edep_of_each_evt,
    energy_threshold,
    n_neighbor_threshold,
    data_dir="$ECAL_CLUSTERING_DATA_DIR",
    output_path=None,
):
    """
    Args:
        filename_edep_of_each_evt (str): input file name under data_dir/processed/
        energy_threshold (float): minimum energy for a seed / cluster membership
        n_neighbor_threshold (int): minimum number of lit neighbors for propagation
        data_dir (str): root directory containing processed/, utilities/, results/
        output_path (str or None): output ROOT file path, defaults to auto-naming
    """
    # Step 1: Load neighbor info
    neighbor_file = ROOT.TFile.Open(
        str(Path(data_dir) / "utilities" / "ecal_neighbor_info.root"), "READ"
    )
    neighbor_tree = neighbor_file.Get("ECALCrystalNeighbors")
    nMod = neighbor_tree.GetEntries()

    neighbors = []
    # vector to receive branch data
    neighbor_vec = vector("int")()
    neighbor_tree.SetBranchAddress("neighbors", neighbor_vec)
    for i in range(nMod):
        neighbor_tree.GetEntry(i)
        neighbors.append([neighbor_vec[j] for j in range(neighbor_vec.size())])

    neighbor_file.Close()
    print(f"Number of crystals: {nMod}")

    # ===============================================================================

    # Step 2: Prepare output
    if output_path is None:
        output_path = Path(
            "ecal_clusters"
            + "_Eth-"
            + str(energy_threshold)
            + "_Nth-"
            + str(n_neighbor_threshold)
            + ".root"
        )
    else:
        output_path = Path(output_path)
    output_file = ROOT.TFile.Open(str(output_path), "RECREATE")
    cluster_tree = ROOT.TTree("ECALClusters", "ECALClusters")
    evt_id = array.array("i", [0])
    mod_id_list = vector("int")()
    cluster_tree.Branch("EvtID", evt_id, "EvtID/I")
    cluster_tree.Branch("ModIDList", mod_id_list)

    # ===============================================================================

    # Step 3: Define per-event clustering function
    def cluster_event(edeps, event_id):
        """Cluster one event given its energy depositions."""
        edeps_full = [0.0] * nMod
        length = min(len(edeps), nMod)
        for i in range(length):
            edeps_full[i] = edeps[i]

        # 0: unclustered and available
        states = [0] * nMod

        # sort crystals by energy descending
        sorted_idx = sorted(range(nMod), key=lambda i: edeps_full[i], reverse=True)

        def is_lit_on(id):
            return edeps_full[id] >= energy_threshold

        def is_propagatable(cid):
            # print(f"NO.{cid}", end=" ")
            if states[cid] != 0 or not is_lit_on(cid):
                return False
            # print(f"n_neighbor = {len(neighbors[cid])}", end=" ")
            lit_nbors = sum(1 for nb in neighbors[cid] if is_lit_on(nb))
            # print(f"lit_nbors = {lit_nbors}")
            return lit_nbors >= n_neighbor_threshold

        seed_ptr = 0

        while True:
            # find next seed: state 0 and energy above threshold
            while seed_ptr < nMod and not (
                states[sorted_idx[seed_ptr]] == 0 and is_lit_on(sorted_idx[seed_ptr])
            ):
                seed_ptr += 1

            if seed_ptr >= nMod:
                break

            seed_id = sorted_idx[seed_ptr]
            # seed_ptr not incremented here; state change will skip it later

            states[seed_id] = 1
            wave = neighbors[seed_id]
            next_wave = []
            cluster_list = [seed_id]
            # level number: seed=1, neighbors=2, ...
            crystal_state = 1

            while wave:
                # print(wave)
                crystal_state += 1
                for cur_id in wave:
                    # assign cluster level
                    if is_propagatable(cur_id):
                        # print("yes")
                        states[cur_id] = crystal_state
                        cluster_list.append(cur_id)
                        # add neighbors to next wave (simple dedup)
                        for nb in neighbors[cur_id]:
                            if nb not in next_wave:
                                next_wave.append(nb)
                # else:
                # print("fuck")

                wave, next_wave = next_wave, []

            # save current cluster
            mod_id_list.clear()
            for cid in cluster_list:
                mod_id_list.push_back(cid)
            evt_id[0] = event_id
            cluster_tree.Fill()

    # ===============================================================================

    # Step 4: Loop over input events
    input_path = Path(data_dir) / "processed" / filename_edep_of_each_evt
    input_file = ROOT.TFile.Open(str(input_path), "READ")
    edep_tree = input_file.Get("EdepOfEachEvt")

    # set branch for Edeps (vector<float>)
    edeps_vec = vector("float")()
    edep_tree.SetBranchAddress("Edeps", edeps_vec)

    nEvents = edep_tree.GetEntries()
    for evt_id_loop in range(nEvents):
        edep_tree.GetEntry(evt_id_loop)
        # edeps_vec can be passed directly (matches C++ signature)
        cluster_event(edeps_vec, evt_id_loop)
        # print(f"========== Event NO.{evt_id_loop} ==========")

    input_file.Close()

    # ===============================================================================

    # Step 5: Write and close output
    output_file.cd()
    cluster_tree.Write()
    output_file.Close()
    print(f"Output written to {output_path}")
    return output_path

    # ===============================================================================


import scoring

# command line interface
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: [filename] [Eth=10] [Nth=2] [output_name] "
              "[expected_lambda] [--event-by-event] [--data-path path]")
        sys.exit(1)

    filename = sys.argv[1]
    E_th = float(sys.argv[2]) if len(sys.argv) >= 3 else 10.0
    N_th = float(sys.argv[3]) if len(sys.argv) >= 4 else 2.0
    output_name = sys.argv[4] if len(sys.argv) >= 5 else None
    exp_lambda = float(sys.argv[5]) if len(sys.argv) >= 6 else None

    event_by_event = "--event-by-event" in sys.argv
    data_path = None
    if "--data-path" in sys.argv:
        idx = sys.argv.index("--data-path")
        if idx + 1 < len(sys.argv):
            data_path = sys.argv[idx + 1]

    output_path = reconstruction(filename, E_th, N_th, output_path=output_name)
    scoring.run(output_path, expected_lambda=exp_lambda,
                event_by_event=event_by_event, data_path=data_path)
