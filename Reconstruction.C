#include "CryInfo.H"
#include "ROOT/RDataFrame.hxx"
#include <string>
#include <print>

// auto main(int argc, char *argv[]) -> int {
//     const auto filename_edep_of_each_evt{argv[1]};
//     const std::string E_th_str{argv[2]};
//     const auto energy_threshold{std::stof(E_th_str)};
//     auto df{ROOT::RDataFrame("EdepOfEachEvt", filename_edep_of_each_evt)};
//     df.Foreach([&](std::vector<float> edeps){ClusterEvent(edeps, energy_threshold);}, {"Edeps"});
//     return 0;
// }

auto Reconstruction(std::string filename_edep_of_each_evt, float energy_threshold = 1.5) -> int{
    ROOT::EnableImplicitMT(0);
    
    std::unique_ptr<TFile> file_neighbors{TFile::Open("ecal_neighbor_info.root", "READ")};
    auto tree_neighbors{static_cast<TTree*>(file_neighbors->Get("ECALCrystalNeighbors"))};
    std::vector<int> *nbor{nullptr};
    tree_neighbors->SetBranchAddress("neighbors", &nbor);
    const auto nMod{tree_neighbors->GetEntries()};

    std::unique_ptr<TFile> file_cluster{TFile::Open("ecal_clusters.root", "RECREATE")};
    std::unique_ptr<TTree> tree_cluster{new TTree("ECALClusters","ECALClusters")};
    tree_cluster->SetDirectory(file_cluster.get());
    int evt_id_cluster{};
    std::vector<int> mod_id_cluster{};
    tree_cluster->Branch("EvtID", &evt_id_cluster);
    tree_cluster->Branch("ModIDList", &mod_id_cluster);

    auto ClusterEvent = [&](std::vector<float> edep_of_cry){
        std::vector<int> wave;
        std::vector<int> next_wave;

        const auto nModHit{edep_of_cry.size()};

        // edep_of_cry is not raw data, use it differntly in actual program
        // std::println("Numbers of crystal hit in this file run: {}", nModHit);
        // std::println("Numbers of existent crystals: {}", nMod);
        auto crystal_info{CryInfo(nMod)};
g        for (int modID{0}; modID < nMod; ++modID){
            if (modID < nModHit){
                crystal_info.SetEdepOf(modID, edep_of_cry.at(modID));
            } else{
                crystal_info.SetEdepOf(modID, 0.0f);
            }
        }

        const auto &sort_idx{crystal_info.GetEdepSortedID()};

        int seed_edep_rank{0};
        for (int seed_count{0}; true; ++seed_count) {
            // exception for crystal with high edep but already clustered
            while (crystal_info.GetStateOf(sort_idx.at(seed_edep_rank)) != 0) {
                ++seed_edep_rank;
            }
            auto seed_id{sort_idx.at(seed_edep_rank)};
            if (crystal_info.GetEdepOf(seed_id) <= energy_threshold)
            break;

            // std::println("seed NO.{}", seed_count);
            mod_id_cluster.clear();

            wave.emplace_back(seed_id);
            int crystal_state{0};
            while (not wave.empty()) {
                // std::print("current wave: ");
                // PrintVec(wave);
                ++crystal_state;
                for (auto id_current : wave) {
                    // set current wave state
                    if (crystal_info.GetStateOf(id_current) == 0
                    and crystal_info.GetEdepOf(id_current) > energy_threshold) {
                        crystal_info.SetStateOf(id_current, crystal_state);
                        mod_id_cluster.emplace_back(id_current);
                        // std::println("SetStateOf ID {} to {}", id_current, crystal_state);
                        tree_neighbors->GetEntry(id_current);
                        // get next wave
                        for (auto id_neighbor : *nbor) {
                            // if condition just for avoiding duplicate elements
                            if (std::find(next_wave.begin(), next_wave.end(), id_neighbor)
                            == next_wave.end()) {
                                next_wave.emplace_back(id_neighbor);
                            }
                        }
                    }
                }
                wave.swap(next_wave);
                next_wave.clear();
            }
            tree_cluster->Fill();
        }
    };
    
    auto df{ROOT::RDataFrame("EdepOfEachEvt", filename_edep_of_each_evt)};
    evt_id_cluster = 0;
    df.Foreach(
    [&](std::vector<float> edeps){
        ClusterEvent(edeps);
        std::println("========== Event NO.{} ==========", evt_id_cluster);
        ++evt_id_cluster;
    }, {"Edeps"});
    
    // std::unique_ptr<TFile> file_edeps{TFile::Open(filename_edep_of_each_evt.c_str(), "READ")};
    // auto tree{static_cast<TTree*>(file_edeps->Get("EdepOfEachEvt"))};
    // std::vector<float> *edeps{nullptr};
    // tree->SetBranchAddress("Edeps", &edeps);
    // for (int evtID{0}; evtID < tree->GetEntries(); ++evtID){
    //     tree->GetEntry(evtID);
    //     std::println("========== Event NO.{} ==========", evtID);
    //     ClusterEvent(*edeps, energy_threshold);
    // }
    
    file_cluster->Write();
    return 0;
}
