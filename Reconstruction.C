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
    auto df{ROOT::RDataFrame("EdepOfEachEvt", filename_edep_of_each_evt)};
    df.Foreach([&](std::vector<float> edeps){ClusterEvent(edeps, energy_threshold);}, {"Edeps"});
    return 0;
}
