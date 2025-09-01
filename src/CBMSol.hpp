#ifndef SRFLP_HPP
#define SRFLP_HPP

#include "PTAPI/include/Problem.h"

using namespace std;

enum Movement {
    REINSERTION,
    TWOOPT,
    SWAP,
    ADJACENTSWAP,
    BLOCK_INSERTION
};

enum Construction {
    GREEDY,
    LKH
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

#endif