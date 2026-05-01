#include "ClusterEvent.H"
#include "ROOT/RDataFrame.hxx"
#include <string>
auto main(int argc, char *argv[]) -> int {
  const auto filename_edep_of_each_evt{argv[1]};
  const std::string E_th_str{argv[2]};
  const auto energy_threshold{std::stof(E_th_str)};
  auto df{ROOT::RDataFrame("EdepOfEachEvt", filename_edep_of_each_evt)};
  df.Foreach([&](std::vector<float> edeps){ClusterEvent(edeps, energy_threshold);}, {"Edeps"});
  return 0;
}

auto Reconstruction(std::string filename_edep_of_each_evt, float energy_threshold) -> int{
    // auto df{ROOT::RDataFrame("EdepOfEachEvt", filename_edep_of_each_evt)};
    // df.Foreach([&](std::vector<float> edeps){ClusterEvent(edeps, energy_threshold);}, {"Edeps"});
    
    std::unique_ptr<TFile> file{TFile::Open(filename_edep_of_each_evt.c_str(), "READ")};
    auto tree{static_cast<TTree*>(file->Get("EdepOfEachEvt"))};
    std::vector<float> *edeps{nullptr};
    tree->SetBranchAddress("Edeps", &edeps);
    for (int evtID{0}; evtID < tree->GetEntries(); ++evtID){
        tree->GetEntry(evtID);
        ClusterEvent(*edeps, energy_threshold);
    }
    
    return 0;
}
