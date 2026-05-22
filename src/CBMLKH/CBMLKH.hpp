#include <vector>

#include "Solution.hpp"

#ifndef CBMLKH_HPP
#define CBMLKH_HPP

using namespace std;

class CBMLKH {
    int threads;

    int l, c;

    vector<vector<int>> binaryMatrix;
    vector<vector<int>> diffMatrix;
    vector<vector<int>> onesToZeros;
    vector<vector<int>> zerosToOnes;
    vector<vector<int>> onesToOnes;
    vector<vector<int>> tspMatrix;

    CBMLKH(int l, int c, int threads);

    int deltaEval(Solution& s);
    int completeEval(Solution& s);
    void computeMatrixes();
};

#endif