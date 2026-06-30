#include "LKHWrapper.hpp"

#include <algorithm>
#include <climits>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace fs = std::filesystem;

LKHWrapper::LKHWrapper(const ColumnStore& columns) : columns(columns) {
    if (!fs::exists(tmpDir)) {
        fs::create_directories(tmpDir);
    }
}

vector<int> LKHWrapper::run(const vector<int>& slice, string instanceName, int maxTime) {
    long execId = getExecutionId();
    string tspFile = tmpDir + instanceName + "_tsp_" + to_string(execId) + ".tsp";
    string parFile = tmpDir + instanceName + "_par_" + to_string(execId) + ".par";
    string resultTourFile = tmpDir + instanceName + "_result_tour_" + to_string(execId) + ".tsp";

    writeTSP(slice, instanceName, tspFile, execId);
    writePar(parFile, tspFile, resultTourFile, maxTime);
    runLKH(parFile);

    vector<int> resultTour = getResultTour(resultTourFile);
    vector<int> mappedTour(resultTour.size());
    for (size_t i = 0; i < resultTour.size(); i++) {
        mappedTour[i] = slice[resultTour[i]];
    }

    return mappedTour;
}

void LKHWrapper::writePar(string parFile, string tspFile, string resultTourFile, int maxTime) {
    ofstream out(parFile);
    if (!out.is_open()) throw runtime_error("Error opening file: " + parFile);

    out << "PROBLEM_FILE = " << tspFile << endl;
    out << "TOUR_FILE = " << resultTourFile << endl;
    out << "MOVE_TYPE = 5" << endl;
    out << "PATCHING_C = 3" << endl;
    out << "PATCHING_A = 2" << endl;
    out << "RUNS = 1" << endl;
    out << "TIME_LIMIT = " << maxTime << endl;
}

void LKHWrapper::writeTSP(const vector<int>& slice, string instanceName, string tspFile, long execId) {
    ofstream out(tspFile);
    if (!out.is_open()) throw runtime_error("Error opening file: " + tspFile);

    out << "NAME : " << instanceName << "_tsp_" << execId << endl;
    out << "TYPE : TSP" << endl;
    out << "DIMENSION : " << slice.size() + 1 << endl;
    out << "EDGE_WEIGHT_TYPE : EXPLICIT" << endl;
    out << "EDGE_WEIGHT_FORMAT : FULL_MATRIX" << endl;
    out << "EDGE_WEIGHT_SECTION" << endl;

    // Node 0 is a dummy depot: its distance to a column equals that column's
    // 1-count, which turns the cycle LKH solves into an open Hamiltonian path.
    out << 0;
    for (size_t j = 0; j < slice.size(); j++) out << " " << columns.onesCount(slice[j]);
    out << endl;

    for (size_t i = 0; i < slice.size(); i++) {
        out << columns.onesCount(slice[i]);
        for (size_t j = 0; j < slice.size(); j++) {
            out << " " << columns.hamming(slice[i], slice[j]);
        }
        out << endl;
    }
}

vector<int> LKHWrapper::getResultTour(string resultTourFile) {
    ifstream sol(resultTourFile);
    if (!sol) throw runtime_error("Error opening solution file: " + resultTourFile);

    vector<int> tour;
    string line;
    bool inTourSection = false;

    while (getline(sol, line)) {
        if (!inTourSection) {
            if (line == "TOUR_SECTION") inTourSection = true;
            continue;
        }
        istringstream iss(line);
        int node;
        while (iss >> node) {
            if (node == -1) {
                inTourSection = false;
                break;
            }
            tour.push_back(node - 2);
        }
    }

    auto depotIt = find(tour.begin(), tour.end(), -1);
    if (depotIt != tour.end()) {
        rotate(tour.begin(), depotIt, tour.end());
        tour.erase(tour.begin());
    }

    return tour;
}

void LKHWrapper::runLKH(string parFile) {
    if (!fs::exists(lkhPath)) throw runtime_error("LKH executable not found: " + lkhPath);
    if (!fs::exists(parFile)) throw runtime_error("PAR file not found: " + parFile);
    string command = lkhPath + " " + parFile + " > /dev/null 2>&1";
    int ret = system(command.c_str());
    if (ret != 0) cerr << "LKH returned error code: " << ret << endl;
}

long LKHWrapper::getExecutionId() {
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<long> dist(1, LONG_MAX);
    return dist(gen);
}

void LKHWrapper::clearTmpDir() {
    if (!fs::exists(tmpDir)) {
        return;
    }

    for (const auto& entry : fs::directory_iterator(tmpDir)) {
        try {
            fs::remove_all(entry.path());
        } catch (const fs::filesystem_error& e) {
            cerr << "Failed to remove " << entry.path() << ": " << e.what() << endl;
        }
    }
}
