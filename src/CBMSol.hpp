#ifndef SRFLP_HPP
#define SRFLP_HPP

#include "PTAPI/include/Problem.h"

struct CBMSol : public solution {
    std::vector<int> sol;
    int cost;
};

inline std::ostream& operator<<(std::ostream& os, const CBMSol& s) {
    os << "Solution: [";
    for (size_t i = 0; i < s.sol.size(); ++i) {
        os << s.sol[i];
        if (i + 1 < s.sol.size()) {
            os << ", ";
        }
    }
    os << "]\nCost: " << s.cost;
    return os;
}

#endif