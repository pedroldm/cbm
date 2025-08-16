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
    vector<vector<int>> diffMatrix;
    vector<vector<int>> onesToZeros;
    vector<vector<int>> zerosToOnes;
    vector<bitset<COLUMNS>> binaryMatrix;

    int movementType;
    double maxTemp;
    double constructionBias;

    random_device rng_device;
    mt19937 mersenne_engine;

    CBMProblem(string filename, int movementType, double constructionBias);
    CBMSol construction();
    CBMSol neighbor(CBMSol sol);
    int completeEval(CBMSol& s);
    int nextInsertion(int& curr, unordered_set<int>& out);
    void computeMatrixes();
    int deltaEval(CBMSol& s);
    int evaluate(CBMSol& s);
    void printS(CBMSol& s);

    template <typename T>
    void printMatrix(const vector<vector<T>>& matrix) {
        for (const auto& row : matrix) {
            for (const auto& elem : row) {
                cout << elem << " ";
            }
            cout << "\n";
        }
    }
};

#endif