#pragma once
#include <algorithm>
#include <iostream>
#include <algorithm>
#include <vector>

#include "TFile.h"
#include "TTree.h"

class CryInfo {
    public:
    explicit CryInfo(int N) : edeps_(N, 0.0), states_(N, -1) {};

    auto GetEdepOf(int i) const -> const float { return edeps_.at(i); };

    auto GetStateOf(int i) const -> const int { return states_.at(i); };

    auto SetEdepOf(int modID, float edep) -> void {
        edeps_[modID] = edep;
        states_[modID] = 0;
    }

    auto SetStateOf(int modID, int state) -> void {
        if (state == 0 or state == -1) {
            std::cout << "[Warning]: Cannot set crystal state to " << state
            << "by function SetStateOf()!\n";
            return;
        }
        states_[modID] = state;
    }

    auto GetEdepSortedID() const -> std::vector<int> {
        std::vector<int> idx(states_.size(), -1);
        for (int n{0}; int &i : idx) {
            if (states_.at(n) == 0)
            i = n++;
        }
        std::sort(idx.begin(), idx.end(),
        [this](int a, int b) { return edeps_.at(a) > edeps_.at(b); });
        return idx;
    };

    private:
    std::vector<float> edeps_;
    // level: unset=-1, unclustered=0, seed=1, seed_neighbor=2, etc.
    std::vector<int> states_;
};

template <typename T>
auto PrintVec(const std::vector<T>& vec) -> void{
    std::cout << "{";
    for (auto i : vec){
        std::cout << i << " ";
    }
    std::cout << "}\n";
}

auto ClusterEvent(const std::vector<float> &edep_of_cry,
float energy_threshold) -> void {
    // PrintVec(edep_of_cry);
    std::vector<int> wave;
    std::vector<int> next_wave;

    std::unique_ptr<TFile> file{TFile::Open("ecal_neighbor_info.root", "READ")};
    auto tree{static_cast<TTree*>(file->Get("ECALCrystalNeighbors"))};
    std::vector<int> *nbor{nullptr};
    tree->SetBranchAddress("neighbors", &nbor);

    // edep_of_cry is not raw data, use it differntly in actual program
    auto crystal_info{CryInfo(edep_of_cry.size())};
    std::cout << "size of vector edep_of_cry: " << edep_of_cry.size() << "\n";
    for (int modID{0}; int edep : edep_of_cry) {
        crystal_info.SetEdepOf(modID++, edep);
    }

    const auto &sort_idx{crystal_info.GetEdepSortedID()};
    // std::cout << sort_idx[0] << "\n";

    int seed_edep_rank{0};
    for (int seed_count{0}; true; ++seed_count) {
        // exception for crystal with high edep but already clustered
        while (crystal_info.GetStateOf(sort_idx.at(seed_edep_rank)) != 0) {
            ++seed_edep_rank;
        }
        // std::cout << "seed_edep_rank of seed" << "[" <<seed_count << "]" << " is " << seed_edep_rank << "\n";
        auto seed_id{sort_idx.at(seed_edep_rank)};
        // std::cout << "seed_id of seed" << "[" << seed_count << "]" << " is " << seed_id << "\n";
        // std::cout << "edep of crystal seed_id" << "[" << seed_id << "]" << " is " << crystal_info.GetEdepOf(seed_id) << "\n";
        if (crystal_info.GetEdepOf(seed_id) < energy_threshold)
        break;

        std::cout << "fuck" << "\n";
        wave.emplace_back(seed_id);
        int crystal_state{0};
        while (not wave.empty()) {
            std::cout << "current wave: ";
            PrintVec(wave);
            ++crystal_state;
            for (auto id_current : wave) {
                // set current wave state
                if (crystal_info.GetStateOf(id_current) == 0
                and crystal_info.GetEdepOf(id_current) >= energy_threshold) {
                    crystal_info.SetStateOf(id_current, crystal_state);
                    std::cout << "SetStateOf ID " << id_current << " to " << crystal_state << "\n";
                    tree->GetEntry(id_current);
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
    }
}
