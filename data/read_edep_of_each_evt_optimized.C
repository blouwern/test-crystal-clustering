#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <map>

int read_edep_of_each_evt_optimized() {
    ROOT::EnableImplicitMT(18);
    ROOT::RDataFrame df("G4Run0/ECALSimHit", "SimMACEPhaseI_20260401.root");

    // Aggregate: build a map that sums Edep for each (EvtID, ModID)
    using Key = std::pair<int, short>;
    using Map = std::map<Key, float>;
    auto aggregator = [](Map &m, int evt, short mod, float edep){
        m[{evt, mod}] += edep;
    };
    auto merger = [](std::vector<Map>& maps){
        if (maps.empty()) return;

        Map& result = maps[0];

        for (size_t i = 1; i < maps.size(); ++i) {
            for (const auto& [key, value] : maps[i]) {
                if (result.count(key) > 0) {
                    throw std::runtime_error(
                    "Duplicate key found: (" + std::to_string(key.first) +
                    ", " + std::to_string(key.second) + ")"
                    );
                }
                result[key] += value;
            }
        }
    };

    auto summedMap = df.Aggregate(aggregator, merger, {"EvtID", "ModID", "Edep"}, Map{});

    // Retrieve the final map (this triggers the event loop)
    Map result = *summedMap;

    // Get dimensions
    int nEvt = *df.Max("EvtID") + 1;
    short nMod = *df.Max("ModID") + 1;

    // Output tree
    TFile *outFile = TFile::Open("edep_of_each_evt.root", "RECREATE");
    TTree *tree = new TTree("EdepOfEachEvt", "Edep of each event");
    std::vector<float> edeps(nMod);
    tree->Branch("Edeps", &edeps);

    // Fill sequentially
    for (int evt = 0; evt < nEvt; ++evt) {
        for (short mod = 0; mod < nMod; ++mod) {
            auto it = result.find({evt, mod});
            edeps[mod] = (it != result.end()) ? it->second : 0.0f;
        }
        tree->Fill();
    }

    outFile->Write();
    outFile->Close();
    return 0;
}
