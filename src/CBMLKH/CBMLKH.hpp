#include <algorithm>
#include <fstream>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "Config.hpp"
#include "Solution.hpp"

#ifndef CBMLKH_HPP
#define CBMLKH_HPP

using namespace std;

class CBMLKH {
    string filePath;
    int threads;

    int l, c;

    vector<vector<bool>> binaryMatrix;
    vector<vector<int>> diffMatrix;
    vector<vector<int>> onesToZeros;
    vector<vector<int>> zerosToOnes;
    vector<vector<int>> onesToOnes;
    vector<vector<int>> tspMatrix;

    random_device rng_device;
    mt19937 mersenne_engine;
    float constructionBias;

   public:
    CBMLKH(const Config& cfg);

    int deltaEval(Solution& s);
    int completeEval(Solution& s);
    Solution greedyConstruction();
    int nextInsertion(int& current, unordered_set<int>& remaining);
    void computeMatrixes();
    vector<int> count1BlocksPerColumn(const Solution& s);
};

#endif