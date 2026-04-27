#include <stdio.h>
#include <vector>
#include <unordered_map>
#include <mutex>
#include <functional>
#include "ROOT/RDataFrame.hxx"
#include "TTree.h"
#include "TFile.h"
#include "ROOT/TThreadedObject.hxx"

// Accumulator structure: for each event, a vector of sums per module
struct EvtModuleAccumulator {
    std::unordered_map<int, std::vector<float>> data; // key = EvtID

    // Add one hit to the accumulator (thread‑safe when used via TThreadedObject)
    void Add(int evtID, short modID, float edep, short nMod) {
        auto it = data.find(evtID);
        if (it == data.end()) {
            // First time we see this event: create a zero‑initialised vector
            auto& vec = data[evtID];
            vec.resize(nMod, 0.0f);
            vec[modID] = edep;
        } else {
            it->second[modID] += edep;
        }
    }

    // Merge another accumulator into this one (for combining thread‑local results)
    void Merge(const EvtModuleAccumulator& other) {
        for (const auto& [evt, otherVec] : other.data) {
            auto it = data.find(evt);
            if (it == data.end()) {
                data[evt] = otherVec;   // copy the whole vector
            } else {
                // Sum module‑wise
                auto& myVec = it->second;
                for (size_t i = 0; i < myVec.size(); ++i) {
                    myVec[i] += otherVec[i];
                }
            }
        }
    }
};

int read_edep_of_each_evt_acc() {
    ROOT::EnableImplicitMT(18);

    // 1. Open the input file and create the RDataFrame
    ROOT::RDataFrame df("G4Run0/ECALSimHit", "SimMACEPhaseI_20260401.root");

    // 2. Determine number of modules (ModID is 0‑based)
    short nMod = static_cast<short>(*df.Max("ModID")) + 1;

    // 3. Thread‑local accumulator: each thread gets its own copy
    ROOT::TThreadedObject<EvtModuleAccumulator> threadedAcc;

    // 4. Process all entries in parallel, filling the thread‑local accumulators
    //    No locks – each entry is handled by a single thread.
    df.Foreach([&](int evtID, short modID, float edep) {
        threadedAcc->Add(evtID, modID, edep, nMod);
    }, {"EvtID", "ModID", "Edep"});

    // 5. Merge all thread‑local accumulators into a single one
    EvtModuleAccumulator globalAcc;
    threadedAcc.Merge(static_cast<std::function<void(EvtModuleAccumulator&)>>([&](EvtModuleAccumulator& acc) {
        globalAcc.Merge(acc);
        // return globalAcc;
    }));
    // globalAcc.Merge(threadedAcc); 
    
    // 6. Write the per‑event module sums to a new TTree
    TFile* outFile = TFile::Open("edep_of_each_evt.root", "RECREATE");
    TTree* tree = new TTree("EdepOfEachEvt", "Edep of each event");

    std::vector<float> edeps(nMod);   // buffer for the current event
    tree->Branch("Edeps", &edeps);

    // Events may be processed out‑of‑order; sort them for output
    std::vector<int> evtIDs;
    evtIDs.reserve(globalAcc.data.size());
    for (const auto& p : globalAcc.data) evtIDs.push_back(p.first);
    std::sort(evtIDs.begin(), evtIDs.end());

    for (int evt : evtIDs) {
        const auto& sums = globalAcc.data[evt];
        // sums size must equal nMod; copy into the branch buffer
        std::copy(sums.begin(), sums.end(), edeps.begin());
        tree->Fill();
    }

    outFile->Write();
    outFile->Close();
    return 0;
}
