#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <iostream>
#include <utility>
#include <vector>

enum BlockMovement { PEAK, INTERVAL, MERGE, RANDOM };
enum Movement { REINSERTION, TWOOPT, SWAP };
enum Construction { GREEDY };

inline const char* toString(BlockMovement movement) {
    switch (movement) {
        case PEAK: return "PEAK";
        case INTERVAL: return "INTERVAL";
        case MERGE: return "MERGE";
        case RANDOM: return "RANDOM";
        default: return "UNKNOWN";
    }
}

class Solution {
public:
    std::vector<int> sol;
    std::vector<int> mE;
    std::vector<int> blocksCount;
    int cost;

    BlockMovement blockMovement;
    Movement movement;
    Construction construction;

    Solution() : cost(0) {}
    Solution(std::vector<int>& sol, int cost) : sol(sol), cost(cost) {}
    Solution(std::vector<int>&& sol, int cost) : sol(std::move(sol)), cost(cost) {}
};

inline std::ostream& operator<<(std::ostream& os, const Solution& s) {
    os << "Cost: " << s.cost
       << " | Block movement: " << toString(s.blockMovement);

    return os;
}

#endif