#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "TTree.h"
#include "TFile.h"

// Aggregate every ECALSimHit hit of an event into a single Edeps vector
// (one entry per crystal, duplicates on the same ModID summed).
//
// Run from the data/ directory; paths are relative to it:
//   root -l -b -q 'read_edep_of_each_evt.C()'
//   root -l -b -q 'read_edep_of_each_evt.C("raw/gamma_iso.root","processed/edep_of_each_evt_gamma_iso.root")'
int read_edep_of_each_evt(
    const std::string& src_file_name = "raw/ep_iso.root",
    const std::string& target_file_name = "processed/edep_of_each_evt_ep_iso.root")
{
    const auto src_tree_name{"G4Run0/ECALSimHit"};
    const auto target_tree_name{"EdepOfEachEvt"};

    ROOT::RDataFrame df(src_tree_name, src_file_name);

    if (*df.Count() == 0) {
        std::printf("[SKIP] %s: tree %s has 0 entries, no output written\n",
                    src_file_name.c_str(), src_tree_name);
        return 2;
    }

    auto nMod{ static_cast<short>(*df.Max("ModID")+1) };
    auto nEvt{static_cast<int>(*df.Max("EvtID")+1)};

    std::map<std::pair<int, short>, float> evt_mod_edep_map;
    df.Foreach([&](int evtID, short modID, float edep){
        evt_mod_edep_map[{evtID, modID}] += edep;
    }, {"EvtID","ModID","Edep"});

    auto file{TFile::Open(target_file_name.c_str(),"RECREATE")};
    auto tree = new TTree(target_tree_name, "Edep of each event");
    std::vector<float> edeps(nMod);
    tree->Branch("Edeps", &edeps);

    // The map is ordered by (evtID, modID), so walk it once instead of doing
    // nEvt * nMod lookups. Output is identical to the per-crystal find().
    auto it{evt_mod_edep_map.begin()};
    const int report_step{ std::max(1, nEvt / 10) };
    for (int evt{0}; evt < nEvt; ++evt){
        if (evt % report_step == 0){
            std::printf("Writing vector edeps of event%d\n", evt);
        }
        std::fill(edeps.begin(), edeps.end(), 0.0f);
        while (it != evt_mod_edep_map.end() && it->first.first == evt){
            edeps[it->first.second] = it->second;
            ++it;
        }
        tree->Fill();
    }

    file->Write();
    file->Close();

    std::printf("[OK] %s -> %s (%d events, %d modules)\n",
                src_file_name.c_str(), target_file_name.c_str(), nEvt, int(nMod));
    return 0;
}
