#include <stdio.h>
#include "ROOT/RDataFrame.hxx"

int test_get_edep(){
    int evtID{ 0 };
    ROOT::RDataFrame df("G4Run0/ECALSimHit","SimMACEPhaseI_20260401.root");
    printf("Grabing edeps of crystals of event%d\n", evtID);
    auto df_evt{ df.Filter([&](int x){return x==evtID;},{"EvtID"}) };

    short nMod{ static_cast<short>(*df.Max("ModID")+1) };
    short modID{ 0 };
    float edep{ };
    while (modID < nMod){
        edep = *(df_evt.Filter([&](short x){return x==modID;},{"ModID"}).Sum<float>("Edep"));
        printf("Crystal NO.%d edep = %f\n",modID,edep);
        ++modID;
    }
    return 0;
}
