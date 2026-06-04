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
        std::vector<int> idx;
        for (size_t i = 0; i < states_.size(); ++i) {
            if (states_[i] == 0) idx.push_back(i);
        }
        std::sort(idx.begin(), idx.end(),
        [this](int a, int b) { return edeps_[a] > edeps_[b]; });
        return idx;
    }
    
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
