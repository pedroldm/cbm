#ifndef SRFLP_HPP
#define SRFLP_HPP

#include "PTAPI/include/Problem.h"

using namespace std;

enum Movement {
    REINSERTION,
    TWOOPT,
    SWAP,
    ONEBLOCKM
};

enum Construction {
    GREEDY,
    LKH,
    ONEBLOCK
};

struct CBMSol : public solution {
    vector<int> sol;
    vector<int> mE;
    Movement movement;
    Construction construction;
    int cost;
};

inline ostream& operator<<(ostream& os, const CBMSol& s) {
    os << "Solution: [";
    for (size_t i = 0; i < s.sol.size(); ++i) {
        os << s.sol[i];
        if (i + 1 < s.sol.size()) {
            os << ", ";
        }
    }
    os << "]" << endl;
    os << "Cost: " << s.cost << endl;
    return os;
}

template <typename T>
int moveHelper(std::vector<T>& v, int from, int to) {
    if (from < 0 || from >= v.size() || to < 0 || to > v.size()) 
        return -1;

    T value = v[from];
    v.erase(v.begin() + from);

    if (from < to) 
        to--;

    v.insert(v.begin() + to, value);

    return to;
}

#endif