#include <stdio.h>
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "TTree.h"
#include "TFile.h"

int read_edep_of_each_evt(){
    ROOT::EnableImplicitMT(18);
    auto file{TFile::Open("edep_of_each_evt.root","RECREATE")};
    ROOT::RDataFrame df("G4Run0/ECALSimHit","SimMACEPhaseI_20260401.root");
    short nMod{ static_cast<short>(*df.Max("ModID")+1) };

    int nEvt{static_cast<int>(*df.Max("EvtID")+1)};
    short modID;
    float edep;
    std::vector<float> edeps(nMod);

    TTree* tree = new TTree("EdepOfEachEvt", "Edep of each event");
    tree->Branch("Edeps", &edeps);

    int evtID{ 0 };
    while (evtID < nEvt){
        // printf("Grabing edeps of crystals of event%d\n", evtID);
        modID = 0;
        while (modID < nMod){
            edep = *(df.Filter([&](int x){return x==evtID;},{"EvtID"})
                       .Filter([&](short x){return x==modID;},{"ModID"})
                       .Sum<float>("Edep"));
            edeps.emplace_back(edep);
            ++modID;
        }
        tree->Fill();
        ++evtID;
        edeps.clear();
    }
    file->Write();
    file->Close();

    return 0;
}
