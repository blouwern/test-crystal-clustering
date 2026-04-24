#include "ClusterEvent.hxx"
#include "ROOT::RDataFrame.hxx"
auto Reconstruction(int argc, char *argv[]) -> int {
  const auto raw_file_name{argv[1]};
  auto df{ROOT::RDataFrame("G4Run0/ECALSimHit", raw_file_name)};

  return 0;
}
