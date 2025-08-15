#ifndef CBMProblem_HPP
#define CBMProblem_HPP

#include <algorithm>
#include <bitset>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <utility>

#include "CBMSol.hpp"
#include "PTAPI/include/Problem.h"

#define COLUMNS 10000

using namespace std;

class CBMProblem : public Problem<CBMSol> {
   public:
    int l, c;
    vector<vector<int>> W;
    vector<bitset<COLUMNS + 2>> B;
    vector<bitset<COLUMNS>> binaryMatrix;

    int movementType;
    double maxTemp;
    double constructionBias;

    random_device rng_device;
    mt19937 mersenne_engine;

    CBMProblem(string filename, int movementType, double constructionBias);
    CBMSol construction();
    CBMSol neighbor(CBMSol sol);
    int evaluate(CBMSol sol);
    void printMatrix(const CBMSol* s = nullptr);
    int nextInsertion(int& curr, unordered_set<int>& out);
    int calculateHamming(int& curr, int& candidate);
    void computeW();
    int deltaEvaluate(CBMSol s);
};

#endif