#include <stdio.h>
#include <vector>
#include <utility>
#include <map>
#include "ROOT/RDataFrame.hxx"
#include "TTree.h"
#include "TFile.h"

int read_edep_of_each_evt(){
    const auto src_tree_name{"G4Run0/ECALSimHit"};
    const auto src_file_name{"SimMACEPhaseI_20260401.root"};
    const auto target_tree_name{"EdepOfEachEvt"};
    const auto target_file_name{"edep_of_each_evt_20260401.root"};
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df(src_tree_name, src_file_name);
    auto nMod{ static_cast<short>(*df.Max("ModID")+1) };
    auto nEvt{static_cast<int>(*df.Max("EvtID")+1)};

    std::map<std::pair<int, short>, float> evt_mod_edep_map;
    df.Foreach([&](int evtID, short modID, float edep){
        evt_mod_edep_map[{evtID, modID}] += edep;
    }, {"EvtID","ModID","Edep"});

    auto file{TFile::Open(target_file_name,"RECREATE")};
    auto tree = new TTree(target_tree_name, "Edep of each event");
    std::vector<float> edeps(nMod);
    tree->Branch("Edeps", &edeps);

    for (int evt{0}; evt < nEvt; ++evt){
        printf("Writing vector edeps of event%d\n", evt);
        for (short mod{0}; mod < nMod; ++mod){
            auto it = evt_mod_edep_map.find({evt, mod});
            if (it == evt_mod_edep_map.end()){
                // printf("[Abnormal]: event%d, mod%d edep not found in map\n", evt, mod);
                edeps[mod] = 0.0f;
                continue;
            }
            // printf("[ok]: event%d, mod%d edep found in map\n", evt, mod);
            edeps[mod] = it->second;
        }
        tree->Fill();
        edeps.assign(nMod, .0f);
    }

    file->Write();
    file->Close();

    return 0;
}
