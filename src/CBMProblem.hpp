#ifndef CBMProblem_HPP
#define CBMProblem_HPP

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "CBMSol.hpp"
#include "PTAPI/include/Problem.h"

class CBMProblem : public Problem<CBMSol> {
   public:
    int n, m;
    std::vector<std::vector<bool>> binaryMatrix;

    int movementType;
    double maxTemp;

    std::random_device rng_device;
    std::mt19937 mersenne_engine;

    CBMProblem(std::string filename, int movementType,
               double maxTempProportion);
    CBMSol construction();
    CBMSol neighbor(CBMSol sol);
    int evaluate(CBMSol& sol);
};

#endif