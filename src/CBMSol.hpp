#ifndef SRFLP_HPP
#define SRFLP_HPP

#include "PTAPI/include/Problem.h"

struct CBMSol : public solution {
    std::vector<int> sol;
    int cost;
};

#endif