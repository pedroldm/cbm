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
#include <omp.h>
#include <atomic>

#include "CBMSol.hpp"
#include "PTAPI/include/Problem.h"
#include <sys/syscall.h>
#include <unistd.h>

#define COLUMNS 500

using namespace std;
namespace fs = filesystem;

class CBMProblem : public Problem<CBMSol> {
   public:
    int l, c;
    vector<CBMSol> initialSolutions;
    vector<tuple<double, int>> diffVec;
    vector<vector<int>> lkhInitialSolutions;
    vector<double> selectionWeights;
    vector<vector<int>> diffMatrix;
    vector<vector<int>> onesToZeros;
    vector<vector<int>> zerosToOnes;
    vector<vector<int>> onesToOnes;
    vector<vector<int>> tspMatrix;
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
    double selectionBias;
    double totalWeight;
    double constructionBias;
    int maxBlockSize;
    int threads;
    int lkhS;
    int lkhMaxTime;
    bool lkhCache;

    int constructionMethod;

    random_device rng_device;
    mt19937 mersenne_engine;

    CBMProblem(string filename, int movementType, int constructionMethod, double constructionBias, double selectionBias, int maxBlockSize, int threads, int lkhS, int lkhMaxTime, bool lkhCache);
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
    void createLKHInitialS();
    int biasedSelection();
    CBMSol oneBlockSearch(CBMSol& s, vector<int> lines = {});
    vector<pair<int, int>> findOneBlocks(CBMSol& s, int line);
    bool moveOneBlockColumns(CBMSol& s, int& currBest, pair<int, int> b1, pair<int, int> b2, bool returnAnyway);
    CBMSol oneBlockMovement(CBMSol& s);
    void moveOneBlockRandomly(CBMSol& s, pair<int, int> b1, pair<int, int> b2);

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