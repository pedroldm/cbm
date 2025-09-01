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
#include <filesystem>
#include <utility>
#include <mutex>
#include <atomic>

#include "CBMSol.hpp"
#include "PTAPI/include/Problem.h"
#include <sys/syscall.h>
#include <unistd.h>

#define COLUMNS 1000

using namespace std;
namespace fs = filesystem;

class CBMProblem : public Problem<CBMSol> {
   public:
    int l, c;
    vector<CBMSol> initialSolutions;
    vector<vector<int>> diffMatrix;
    vector<vector<int>> onesToZeros;
    vector<vector<int>> zerosToOnes;
    vector<vector<int>> onesToOnes;
    vector<bitset<COLUMNS>> binaryMatrix;
    filesystem::path currPath;
    filesystem::path tspPath;
    filesystem::path solPath;
    filesystem::path inputTourPath;
    filesystem::path outputTourPath;
    filesystem::path parPath;
    filesystem::path lkhPath;
    string instanceName;

    mutex ISMutex;
    atomic<int> ISCount;

    int movementType;
    double maxTemp;
    double constructionBias;
    int maxBlockSize;
    int threads;

    random_device rng_device;
    mt19937 mersenne_engine;

    CBMProblem(string filename, int movementType, double constructionBias, int maxBlockSize, int threads);
    CBMSol construction();
    CBMSol greedyConstruction();
    vector<int> lkhConstruction();
    CBMSol neighbor(CBMSol sol);
    int completeEval(CBMSol& s);
    int nextInsertion(int& curr, unordered_set<int>& out);
    void computeMatrixes();
    int deltaEval(CBMSol& s);
    int evaluate(CBMSol& s);
    void printS(CBMSol& s);
    void toTSP(string sufix = "");
    vector<int> fromTSP(string sufix = "");
    void initialTour(string sufix = "");
    void runLKH(string sufix = "");

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