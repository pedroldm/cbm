// Standalone LKH baseline for CBM: solve the whole instance as one open-path TSP.
//
// The TSP is the one LKHWrapper builds for a sub-segment, applied to every
// column: a depot node whose distance to a column is its 1-count, plus Hamming
// distances between columns, so a tour's length is twice the number of 1-blocks
// of the column order it induces. It is written as UPPER_ROW, which holds the
// same symmetric matrix in half the disk space of FULL_MATRIX.
//
// Usage:
//   lkh_standalone --filePath=<instance> --lkhPath=<LKH> --seed=<n>
//                  [--timeLimit=<s>] [--outputPath=<json>] [--workDir=<dir>]
//                  [--runs=1]
//
// MOVE_TYPE / PATCHING_C / PATCHING_A are fixed project-wide (cbm_lkh_params.h).

#include <unistd.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../IO/json.hpp"
#include "../common/cbm_lkh_params.h"
#include "../common/cbm_spawn.h"
#include "ColumnStore.hpp"

using namespace std;
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace {

struct Options {
    string filePath, outputPath, lkhPath, workDir;
    unsigned long seed = 1;
    double timeLimit = 7200.0;
    int runs = 1;
};

double secondsSince(chrono::steady_clock::time_point start) { return chrono::duration<double>(chrono::steady_clock::now() - start).count(); }

Options parseArgs(int argc, char* argv[]) {
    Options opt;
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        auto take = [&](const string& key, auto setter) {
            if (arg.rfind(key, 0) != 0) return false;
            try {
                setter(arg.substr(key.size()));
            } catch (const exception&) {
                throw runtime_error("invalid value in " + arg);
            }
            return true;
        };
        bool known =
            take("--filePath=", [&](const string& v) { opt.filePath = v; }) || take("--outputPath=", [&](const string& v) { opt.outputPath = v; }) ||
            take("--lkhPath=", [&](const string& v) { opt.lkhPath = v; }) || take("--workDir=", [&](const string& v) { opt.workDir = v; }) ||
            take("--seed=", [&](const string& v) { opt.seed = stoul(v); }) ||
            take("--timeLimit=", [&](const string& v) { opt.timeLimit = stod(v); }) ||
            take("--runs=", [&](const string& v) { opt.runs = stoi(v); });
        if (!known) throw runtime_error("unknown option " + arg);
    }
    if (opt.filePath.empty()) throw runtime_error("missing --filePath=<instance>");
    if (opt.lkhPath.empty()) {
        if (const char* env = getenv("LKH_PATH")) opt.lkhPath = env;
    }
    if (opt.lkhPath.empty() || access(opt.lkhPath.c_str(), X_OK) != 0)
        throw runtime_error("LKH not executable: pass --lkhPath=<binary> (or set LKH_PATH)");
    // LKH reads SEED as unsigned; 0 is legal there but kept out for symmetry with Linkern.
    if (opt.seed == 0 || opt.seed > 2147483647UL) throw runtime_error("--seed must be in [1, 2^31 - 1]");
    if (opt.timeLimit <= 0 || opt.runs < 1) throw runtime_error("--timeLimit must be > 0 and --runs >= 1");
    return opt;
}

void writeTSP(const ColumnStore& columns, const string& tspFile) {
    const int n = columns.cols();
    ofstream out(tspFile);
    out << "NAME : cbm\nTYPE : TSP\nDIMENSION : " << n + 1 << "\n"
        << "EDGE_WEIGHT_TYPE : EXPLICIT\nEDGE_WEIGHT_FORMAT : UPPER_ROW\nEDGE_WEIGHT_SECTION\n";
    for (int j = 0; j < n; j++) out << columns.onesCount(j) << (j + 1 < n ? ' ' : '\n');  // depot row
    for (int i = 0; i + 1 < n; i++)
        for (int j = i + 1; j < n; j++) out << columns.hamming(i, j) << (j + 1 < n ? ' ' : '\n');
    out << "EOF\n";
    if (!out) throw runtime_error("cannot write " + tspFile + " (disk full?)");
}

// 0-based column order of the open path (depot removed, path starts after it).
vector<int> readTour(const string& tourFile) {
    ifstream in(tourFile);
    if (!in) throw runtime_error("LKH wrote no tour file " + tourFile);
    vector<int> nodes;
    string line;
    bool inSection = false;
    while (getline(in, line)) {
        if (!inSection) {
            inSection = line == "TOUR_SECTION";
            continue;
        }
        istringstream iss(line);
        int node;
        bool done = false;
        while (iss >> node) {
            if (node == -1) {
                done = true;
                break;
            }
            nodes.push_back(node - 2);  // depot (node 1) -> -1
        }
        if (done) break;
    }
    vector<int> path;
    size_t depot = 0;
    while (depot < nodes.size() && nodes[depot] != -1) depot++;
    if (depot == nodes.size()) throw runtime_error("LKH tour misses the depot");
    for (size_t k = 1; k < nodes.size(); k++) path.push_back(nodes[(depot + k) % nodes.size()]);
    return path;
}

struct LogSummary {
    bool timeLimitHit = false;
    json preprocessingTime = nullptr, runTime = nullptr, timeToBest = nullptr;
    long improvements = 0;
};

// Pulls timings out of LKH's TRACE_LEVEL 1 output.
LogSummary parseLog(const string& logFile) {
    LogSummary s;
    ifstream in(logFile);
    string line;
    while (getline(in, line)) {
        double t;
        long long cost;
        int k;
        if (line.find("Time limit exceeded") != string::npos) s.timeLimitHit = true;
        if (sscanf(line.c_str(), "Preprocessing time = %lf", &t) == 1) s.preprocessingTime = t;
        if (sscanf(line.c_str(), "Run %d: Cost = %lld, Time = %lf", &k, &cost, &t) == 3) s.runTime = t;
        if (sscanf(line.c_str(), "* %d: Cost = %lld, Time = %lf", &k, &cost, &t) == 3) {
            s.timeToBest = t;
            s.improvements++;
        }
    }
    return s;
}

void writeAtomically(const string& path, const string& content) {
    const string tmp = path + ".tmp";
    {
        ofstream out(tmp);
        out << content;
        out.flush();
        if (!out) throw runtime_error("cannot write " + tmp);
    }
    fs::rename(tmp, path);
}

}  // namespace

int main(int argc, char* argv[]) {
    const auto start = chrono::steady_clock::now();
    string ownedWorkDir;
    int status = 0;
    try {
        Options opt = parseArgs(argc, argv);
        if (opt.workDir.empty()) {
            string pattern = (fs::temp_directory_path() / "lkh-standalone-XXXXXX").string();
            if (!mkdtemp(pattern.data())) throw runtime_error("cannot create a work directory");
            opt.workDir = ownedWorkDir = pattern;
        } else {
            fs::create_directories(opt.workDir);
        }

        ColumnStore columns(opt.filePath);
        const string tspFile = opt.workDir + "/problem.tsp", parFile = opt.workDir + "/problem.par";
        const string tourFile = opt.workDir + "/problem.tour", logFile = opt.workDir + "/lkh.log";

        writeTSP(columns, tspFile);
        const double buildSeconds = secondsSince(start);
        {
            ofstream par(parFile);
            par << "PROBLEM_FILE = " << tspFile << "\nTOUR_FILE = " << tourFile << "\n"
                << "MOVE_TYPE = " << CBM_LKH_MOVE_TYPE << "\nPATCHING_C = " << CBM_LKH_PATCHING_C << "\nPATCHING_A = " << CBM_LKH_PATCHING_A << "\n"
                << "RUNS = " << opt.runs << "\nTIME_LIMIT = " << opt.timeLimit << "\nSEED = " << opt.seed << "\nTRACE_LEVEL = 1\n";
        }

        string exe = opt.lkhPath, parArg = parFile;
        char* lkhArgv[] = {exe.data(), parArg.data(), nullptr};
        const auto lkhStart = chrono::steady_clock::now();
        const int rc = cbm_run(lkhArgv, logFile.c_str());
        const double lkhSeconds = secondsSince(lkhStart);
        fs::remove(tspFile);  // the largest file by far; drop it as soon as possible
        if (rc != 0) throw runtime_error("LKH exited with status " + to_string(rc) + "; see " + logFile);

        vector<int> path = readTour(tourFile);
        const int n = columns.cols();
        vector<bool> seen(n, false);
        if (static_cast<int>(path.size()) != n) throw runtime_error("LKH tour has " + to_string(path.size()) + " columns, expected " + to_string(n));
        for (int c : path) {
            if (c < 0 || c >= n || seen[c]) throw runtime_error("LKH tour is not a permutation of the columns");
            seen[c] = true;
        }

        long blocks = columns.onesCount(path[0]);
        for (int k = 1; k < n; k++) blocks += columns.zerosToOnes(path[k - 1], path[k]);

        LogSummary log = parseLog(logFile);
        vector<int> permutation(path.size());
        for (size_t k = 0; k < path.size(); k++) permutation[k] = path[k] + 1;

        json out = {
            {"method", "LKH"},
            {"instance", opt.filePath},
            {"rows", columns.rows()},
            {"cols", n},
            {"seed", opt.seed},
            {"parameters",
             {{"time_limit_s", opt.timeLimit},
              {"runs", opt.runs},
              {"move_type", CBM_LKH_MOVE_TYPE},
              {"patching_c", CBM_LKH_PATCHING_C},
              {"patching_a", CBM_LKH_PATCHING_A}}},
            {"best_blocks", blocks},
            {"time_limit_reached", log.timeLimitHit},
            {"stop_reason", log.timeLimitHit ? "time_limit" : "trials"},
            {"elapsed_s", secondsSince(start)},
            {"tsp_build_time_s", buildSeconds},
            {"lkh_time_s", lkhSeconds},
            {"lkh_preprocessing_time_s", log.preprocessingTime},
            {"lkh_run_time_s", log.runTime},
            {"lkh_time_to_best_s", log.timeToBest},
            {"lkh_improvements", log.improvements},
            {"permutation", permutation},
        };
        if (opt.outputPath.empty())
            cout << out.dump(2) << endl;
        else
            writeAtomically(opt.outputPath, out.dump(2) + "\n");
        fs::remove(tourFile);
        fs::remove(parFile);
    } catch (const exception& e) {
        cerr << "lkh_standalone error: " << e.what() << endl;
        status = 1;
    }
    if (!ownedWorkDir.empty() && status == 0) fs::remove_all(ownedWorkDir);
    return status;
}
