import ROOT


def get_N_mod():
    file_read = ROOT.TFile.Open(
        "$ECAL_CLUSTERING_DATA_DIR/utilities/ecal_neighbor_info.root", "READ"
    )
    tree_read = file_read.Get("ECALCrystalNeighbors")
    n_entries = tree_read.GetEntries()
    return n_entries
