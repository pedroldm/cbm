#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <vector>

using namespace std;

enum Movement { REINSERTION, TWOOPT, SWAP, ONEBLOCKM };
enum Construction { GREEDY, LKH, ONEBLOCK };

class Solution {
   public:
    vector<int> sol;
    vector<int> mE;
    int cost;

    Movement movement;
    Construction construction;

    Solution() : cost(0) {}
    Solution(const vector<int>& sol, int cost) : sol(sol), cost(cost) {}
    Solution(vector<int>&& sol, int cost) : sol(std::move(sol)), cost(cost) {}

    template <typename T>
    static int moveHelper(std::vector<T>& v, int from, int to) {
        if (from < 0 || from >= static_cast<int>(v.size()) || to < 0 || to > static_cast<int>(v.size())) {
            return -1;
        }

        T value = std::move(v[from]);
        v.erase(v.begin() + from);

        if (from < to) {
            to--;
        }

        v.insert(v.begin() + to, std::move(value));

        return to;
    }
};

#endif