import ROOT
from ROOT import std, vector
import os
import sys
from pathlib import Path
import array


def reconstruction(
    filename_edep_of_each_evt,
    energy_threshold_seed,
    energy_threshold_wave,
    n_neighbor_threshold,
    n_fallback_chance=1,
    neighbor_file_name="ecal_neighbor_info.root",
    data_dir=None,
    output_path=None,
):
    """
    Args:
        filename_edep_of_each_evt (str): input file name under data_dir/processed/
        energy_threshold_seed (float): minimum energy for a seed
        energy_threshold_wave (float): minimum energy to join the wavefront
        n_neighbor_threshold (int): minimum number of lit neighbors for a seed
        data_dir (str): data directory containing processed/, utilities/, results/
        output_path (str or None): output ROOT file path, defaults to auto-naming
    """
    # The wavefront threshold is meant to be the looser one: a high seed
    # threshold suppresses noise seeds while a low wave threshold lets the soft
    # shower halo join the cluster.
    if energy_threshold_seed <= energy_threshold_wave:
        print(
            f"[WARNING] energy_threshold_seed ({energy_threshold_seed}) should be "
            f"greater than energy_threshold_wave ({energy_threshold_wave}); "
            "the wavefront threshold is meant to be the looser one."
        )

    if data_dir is None:
        data_dir = os.environ.get("ECAL_CLUSTERING_DATA_DIR", ".")

    # Step 1: Load neighbor info
    neighbor_file = ROOT.TFile.Open(
        str(Path(data_dir) / "utilities" / neighbor_file_name), "READ"
    )
    neighbor_tree = neighbor_file.Get("ECALCrystalNeighbors")
    nMod = neighbor_tree.GetEntries()

    neighbors = []
    sec_neighbors = []
    # vector to receive branch data
    neighbor_vec = vector("int")()
    neighbor_tree.SetBranchAddress("neighbors", neighbor_vec)
    # 'second_neighbors' is optional: without it the gap-jumping fallback is a
    # no-op rather than a crash (SetBranchAddress on a missing branch throws)
    has_second = any(b.GetName() == "second_neighbors"
                     for b in neighbor_tree.GetListOfBranches())
    sec_neighbor_vec = vector("int")()
    if has_second:
        neighbor_tree.SetBranchAddress("second_neighbors", sec_neighbor_vec)
    else:
        print(f"[INFO] {neighbor_file_name} has no 'second_neighbors' branch; "
              "gap-jumping fallback disabled")
    for i in range(nMod):
        neighbor_tree.GetEntry(i)
        neighbors.append([neighbor_vec[j] for j in range(neighbor_vec.size())])
        sec_neighbors.append(
            [sec_neighbor_vec[j] for j in range(sec_neighbor_vec.size())]
        )
    neighbor_file.Close()
    print(f"Number of crystals: {nMod}")

    # ===============================================================================

    # Step 2: Prepare output
    if output_path is None:
        output_path = Path(
            "ecal_clusters"
            + "_Eth_s-"
            + str(energy_threshold_seed)
            + "_Eth_w-"
            + str(energy_threshold_wave)
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
            return edeps_full[id] >= energy_threshold_wave

        def is_seedable(cid):
            if states[cid] != 0 or edeps_full[cid] < energy_threshold_seed:
                return False
            lit_nbors = sum(1 for nb in neighbors[cid] if is_lit_on(nb))
            return lit_nbors >= n_neighbor_threshold

        def is_propagatable(cid):
            # must still be unclustered, otherwise an already-assigned crystal
            # keeps re-entering the wavefront and the loop never terminates
            return states[cid] == 0 and is_lit_on(cid)

        seed_ptr = 0

        while True:
            # find next seed: state 0 and energy above threshold
            while seed_ptr < nMod and not (
                states[sorted_idx[seed_ptr]] == 0 and is_seedable(sorted_idx[seed_ptr])
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
            fallback_count = 0

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
                # special case catching (directly fall back to second neighbors)
                if not next_wave and fallback_count < n_fallback_chance:
                    next_wave = sec_neighbors[cur_id]
                    fallback_count += 1
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
        print(
            "Usage: [filename] [Eth_seed=10] [Eth_wave=5] [Nth=3] [output_name] "
            "[expected_lambda] [--event-by-event] [--data-path path] [--data-dir path]"
        )
        sys.exit(1)

    filename = sys.argv[1]
    E_th_seed = float(sys.argv[2]) if len(sys.argv) >= 3 else 10.0
    E_th_wave = float(sys.argv[3]) if len(sys.argv) >= 4 else 5.0
    N_th = int(float(sys.argv[4])) if len(sys.argv) >= 5 else 3
    output_name = sys.argv[5] if len(sys.argv) >= 6 else None
    exp_lambda = float(sys.argv[6]) if len(sys.argv) >= 7 else None

    def _option(flag):
        if flag in sys.argv:
            idx = sys.argv.index(flag)
            if idx + 1 < len(sys.argv):
                return sys.argv[idx + 1]
        return None

    event_by_event = "--event-by-event" in sys.argv
    data_path = _option("--data-path")
    data_dir = _option("--data-dir")

    output_path = reconstruction(
        filename, E_th_seed, E_th_wave, N_th,
        data_dir=data_dir, output_path=output_name,
    )
    scoring.run(
        output_path,
        expected_lambda=exp_lambda,
        event_by_event=event_by_event,
        data_path=data_path,
    )
