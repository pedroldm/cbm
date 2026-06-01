#include <algorithm>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Config.hpp"
#include "Solution.hpp"
#include "Trajectory.hpp"
#include "Validator.hpp"

#ifndef CBMLKH_HPP
#define CBMLKH_HPP

static constexpr const char* LKH_EXEC_PATH = "/home/pdamasceno/cbm/src/LKH3";

using namespace std;

struct DenseSegment {
    int start;
    int end;
    int sum;
    double avg;
    double score;
};

class CBMLKH {
   public:
    Config cfg;
    Validator validator;

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
    bool lkhCache = false;

    random_device rng_device;
    mt19937 mersenne_engine;
    float constructionBias;

    CBMLKH(const Config& cfg);

    void run();
    Trajectory LKHILS(Solution& initial);
    Solution ILSNeighbor(Solution s);
    int deltaEval(Solution& s);
    int completeEval(Solution& s);
    Solution greedyConstruction();
    int nextInsertion(int& current, unordered_set<int>& remaining);
    void computeMatrixes();
    void countBlocksPerColumn(Solution& s);
    vector<DenseSegment> findDenseSegments(Solution s, int minSize, int maxSize, double minScore);
};

#endif