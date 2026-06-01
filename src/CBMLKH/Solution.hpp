#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <vector>

using namespace std;

enum Movement { REINSERTION, TWOOPT, SWAP };
enum Construction { GREEDY };

class Solution {
   public:
    vector<int> sol;
    vector<int> mE;
    vector<int> blocksCount;
    int cost;

    Movement movement;
    Construction construction;

    Solution() : cost(0) {}
    Solution(vector<int>& sol, int cost) : sol(sol), cost(cost) {}
    Solution(vector<int>&& sol, int cost) : sol(move(sol)), cost(cost) {}
};

#endif