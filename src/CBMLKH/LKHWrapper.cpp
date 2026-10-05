#include "LKHWrapper.hpp"

#include <unistd.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "../common/cbm_lkh_params.h"
#include "../common/cbm_spawn.h"

using namespace std;

namespace fs = std::filesystem;

void LKHWrapper::configure(const string& executablePath, const string& tmpParent) {
    if (access(executablePath.c_str(), X_OK) != 0) throw runtime_error("LKH executable not found or not executable: " + executablePath);
    lkhPath = executablePath;

    fs::create_directories(tmpParent);
    string pattern = (fs::path(tmpParent) / "cbmlkh-XXXXXX").string();
    if (!mkdtemp(pattern.data())) throw runtime_error("Cannot create a scratch directory under " + tmpParent);
    workDir = pattern;
}

void LKHWrapper::cleanup() {
    if (workDir.empty()) return;
    error_code ignored;
    fs::remove_all(workDir, ignored);
    workDir.clear();
}

LKHWrapper::LKHWrapper(const ColumnStore& columns) : columns(columns) {}

LKHResult LKHWrapper::run(const vector<int>& slice, int maxTime, uint32_t seed) const {
    if (workDir.empty()) throw runtime_error("LKHWrapper::configure() was not called");

    const string base = workDir + "/call" + to_string(nextCallId++);
    const string tspFile = base + ".tsp", parFile = base + ".par", tourFile = base + ".tour", logFile = base + ".log";

    writeTSP(slice, tspFile);
    {
        ofstream par(parFile);
        par << "PROBLEM_FILE = " << tspFile << "\n"
            << "TOUR_FILE = " << tourFile << "\n"
            << "MOVE_TYPE = " << CBM_LKH_MOVE_TYPE << "\n"
            << "PATCHING_C = " << CBM_LKH_PATCHING_C << "\n"
            << "PATCHING_A = " << CBM_LKH_PATCHING_A << "\n"
            << "RUNS = 1\n"
            << "TIME_LIMIT = " << maxTime << "\n"
            << "SEED = " << seed << "\n";
        if (!par) throw runtime_error("Error writing " + parFile);
    }

    string exe = lkhPath, parArg = parFile;
    char* argv[] = {exe.data(), parArg.data(), nullptr};
    const int status = cbm_run(argv, logFile.c_str());
    if (status != 0) {
        // The scratch directory is removed at exit, so carry the log in the error.
        ifstream log(logFile);
        string tail, line;
        while (getline(log, line)) tail = tail.size() > 2000 ? line + "\n" : tail + line + "\n";
        throw runtime_error("LKH failed with status " + to_string(status) + "; output:\n" + tail);
    }

    LKHResult result;
    // TRACE_LEVEL defaults to 1, at which LKH reports a binding TIME_LIMIT.
    ifstream log(logFile);
    for (string line; getline(log, line);) {
        if (line.find("Time limit exceeded") != string::npos) {
            result.timeLimitHit = true;
            break;
        }
    }

    vector<int> local = readTour(tourFile);
    if (local.size() != slice.size())
        throw runtime_error("LKH returned " + to_string(local.size()) + " columns, expected " + to_string(slice.size()));
    vector<bool> seen(slice.size(), false);
    result.tour.reserve(local.size());
    for (int idx : local) {
        if (idx < 0 || idx >= static_cast<int>(slice.size()) || seen[idx]) throw runtime_error("LKH returned an invalid tour in " + tourFile);
        seen[idx] = true;
        result.tour.push_back(slice[idx]);  // local node index -> global column id
    }

    for (const string& f : {tspFile, parFile, tourFile, logFile}) fs::remove(f);
    return result;
}

void LKHWrapper::writeTSP(const vector<int>& slice, const string& tspFile) const {
    ofstream out(tspFile);
    if (!out.is_open()) throw runtime_error("Error opening file: " + tspFile);

    out << "NAME : " << fs::path(tspFile).stem().string() << "\n";
    out << "TYPE : TSP\n";
    out << "DIMENSION : " << slice.size() + 1 << "\n";
    out << "EDGE_WEIGHT_TYPE : EXPLICIT\n";
    out << "EDGE_WEIGHT_FORMAT : FULL_MATRIX\n";
    out << "EDGE_WEIGHT_SECTION\n";

    // Node 0 is a dummy depot: its distance to a column equals that column's
    // 1-count, which turns the cycle LKH solves into an open Hamiltonian path.
    out << 0;
    for (size_t j = 0; j < slice.size(); j++) out << " " << columns.onesCount(slice[j]);
    out << "\n";

    for (size_t i = 0; i < slice.size(); i++) {
        out << columns.onesCount(slice[i]);
        for (size_t j = 0; j < slice.size(); j++) out << " " << columns.hamming(slice[i], slice[j]);
        out << "\n";
    }
    out << "EOF\n";
    if (!out) throw runtime_error("Error writing " + tspFile);
}

// Parses the TOUR_SECTION into 0-based column positions of the slice, starting
// right after the depot (LKH node 1).
vector<int> LKHWrapper::readTour(const string& tourFile) {
    ifstream sol(tourFile);
    if (!sol) throw runtime_error("Error opening solution file: " + tourFile);

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
    if (depotIt == tour.end()) throw runtime_error("LKH tour misses the depot: " + tourFile);
    rotate(tour.begin(), depotIt, tour.end());
    tour.erase(tour.begin());
    return tour;
}
