#ifndef SRFLP_HPP
#define SRFLP_HPP

#include "PTAPI/include/Problem.h"

using namespace std;

enum Movement {
    REINSERTION,
    TWOOPT,
    SWAP
};

struct CBMSol : public solution {
    vector<int> sol;
    tuple<int, int, int, int, int> mE; /* j-1, j+1, k-1, k, e */
    Movement movement;
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
    os << "mE: " << get<0>(s.mE) << get<1>(s.mE) << get<2>(s.mE) << get<3>(s.mE) << get<4>(s.mE) << endl;
    return os;
}

#endif