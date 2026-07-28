import ROOT

ROOT.EnableImplicitMT()
file_read = ROOT.TFile.Open("ecal_neighbor_info.root", "READ")
tree_read = file_read.Get("ECALCrystalNeighbors")

file_write = ROOT.TFile.Open("ecal_neighbor_info_added.root", "RECREATE")
tree_write = ROOT.TTree("ECALCrystalNeighbors", "ECALCrystalNeighbors")

vec_nbor_entry = ROOT.std.vector("int")()
tree_write.Branch("neighbors", vec_nbor_entry)

vec_sec_nbor_entry = ROOT.std.vector("int")()
tree_write.Branch("second_neighbors", vec_sec_nbor_entry)

n_entries = tree_read.GetEntries()

for modid in range(n_entries):
    tree_read.GetEntry(modid)
    nbor = list(tree_read.neighbors)

    sec_nbor_set = set()
    for nbor_id in nbor:
        tree_read.GetEntry(nbor_id)
        for sec_id in tree_read.neighbors:
            if sec_id != modid and sec_id not in nbor:
                sec_nbor_set.add(sec_id)

    sec_nbor_list = list(sec_nbor_set)

    vec_nbor_entry.clear()
    for x in nbor:
        vec_nbor_entry.push_back(x)

    vec_sec_nbor_entry.clear()
    for x in sec_nbor_list:
        vec_sec_nbor_entry.push_back(x)

    tree_write.Fill()

file_write.Write()
file_write.Close()
file_read.Close()
