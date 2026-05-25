#ifndef LKHWrapper_HPP
#define LKHWrapper_HPP

#include <algorithm>
#include <climits>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

class LKHWrapper {
    inline static string lkhPath = "/home/pedroldm/MSc/cbm/src/LKH3/LKH";
    inline static string tmpDir = "/tmp/LKH/";

    int c, l;
    vector<vector<int>> tspMatrix;

   public:
    LKHWrapper(int c, int l, const vector<vector<int>>& tspMatrix);
    vector<int> run(const vector<int>& slice, string instanceName, int maxTime);
    void writeTSP(const vector<int>& slice, string instanceName, string tspFile, long execId);
    void writeInitialTour(const vector<int>& slice, string instanceName, string tourFile, long execId);
    void writePar(string parFile, string tspFile, string InitialTourFile, string resultTourFile, int maxTime);
    void runLKH(string parFile);
    long getExecutionId();
    vector<int> getResultTour(string resultTourFile);
};

#endif