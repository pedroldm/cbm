#include <algorithm>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "Config.hpp"
#include "Solution.hpp"

#ifndef CBMLKH_HPP
#define CBMLKH_HPP

static constexpr const char* LKH_EXEC_PATH = "/home/pdamasceno/cbm/src/LKH3";

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
    vector<int> onesSum;
    filesystem::path currPath;
    filesystem::path lkhPath;
    string instanceName;
    int lkhMaxTime;
    bool lkhCache = false;

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
    void countBlocksPerColumn(Solution& s);
    void toTSP(const string& suffix = "", const vector<int>& subset = {});
    void initialTour(const string& suffix = "", const vector<int>& subset = {});
    void runLKH(string suffix);
    vector<int> fromTSP(const string& suffix = "", const vector<int>& subset = {});
};

#endif