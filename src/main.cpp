#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "CBMProblem.hpp"
#include "IO/json.hpp"
#include "PTAPI/include/PT.h"

using json = nlohmann::json;
using Clock = chrono::high_resolution_clock;
using Duration = chrono::duration<double>;

// ─── Config ──────────────────────────────────────────────────────────────────

struct Config {
    float tempMin = 0.05f;
    float tempMax = 2.0f;
    int tempL = 4;
    int tempD = 4;
    int tempUpdate = 3;
    int upType = 1;
    float MKL = 400;
    int PTL = 2000;
    double constructionBias = 2.5;
    int constructionMethod = 1;
    double selectionBias = 1.0;
    int movementType = 4;
    int maxBlockSize = 3;
    int threads = 1;
    int lkhS = 0;
    int lkhMaxTime = 120;
    bool irace = false;
    bool lkhCache = false;
    string filePath;
    string tspPath = "./instances/tsp/";
};

// ─── Argument parsing ─────────────────────────────────────────────────────────

template <typename T>
T parseValue(const string& raw) {
    T val;
    istringstream(raw) >> val;
    return val;
}

bool parseBool(const string& key, const string& value) {
    if (value == "true" || value == "1") return true;
    if (value == "false" || value == "0") return false;
    throw invalid_argument("Invalid value for --" + key + " (expected true/false)");
}

// Returns the part after "prefix" if arg starts with it, else "".
string stripPrefix(const string& arg, const string& prefix) {
    if (arg.rfind(prefix, 0) == 0) return arg.substr(prefix.size());
    return {};
}

Config parseArgs(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];

        // Macro-like helper: try to match a flag and assign a parsed value.
        auto tryParse = [&]<typename T>(const string& flag, T& field) -> bool {
            string val = stripPrefix(arg, "--" + flag + "=");
            if (val.empty()) return false;
            field = parseValue<T>(val);
            return true;
        };

        if (tryParse("tempMin", cfg.tempMin)) {
        } else if (tryParse("tempMax", cfg.tempMax)) {
        } else if (tryParse("tempL", cfg.tempL)) {
        } else if (tryParse("tempD", cfg.tempD)) {
        } else if (tryParse("tempUpdate", cfg.tempUpdate)) {
        } else if (tryParse("upType", cfg.upType)) {
        } else if (tryParse("MKL", cfg.MKL)) {
        } else if (tryParse("PTL", cfg.PTL)) {
        } else if (tryParse("constructionBias", cfg.constructionBias)) {
        } else if (tryParse("constructionMethod", cfg.constructionMethod)) {
        } else if (tryParse("selectionBias", cfg.selectionBias)) {
        } else if (tryParse("movementType", cfg.movementType)) {
        } else if (tryParse("maxBlockSize", cfg.maxBlockSize)) {
        } else if (tryParse("threads", cfg.threads)) {
        } else if (tryParse("lkhS", cfg.lkhS)) {
        } else if (tryParse("lkhMaxTime", cfg.lkhMaxTime)) {
        } else if (tryParse("filePath", cfg.filePath)) {
        } else if (string v = stripPrefix(arg, "--irace="); !v.empty())
            cfg.irace = parseBool("irace", v);
        else if (string v = stripPrefix(arg, "--lkhCache="); !v.empty())
            cfg.lkhCache = parseBool("lkhCache", v);
        else
            throw runtime_error("Unknown argument: " + arg);
    }
    return cfg;
}

// ─── Output helpers ───────────────────────────────────────────────────────────

string constructionLabel(int construction) {
    switch (construction) {
        case GREEDY:
            return "GREEDY";
        case ONEBLOCK:
            return "ONEBLOCK";
        default:
            return "LKH";
    }
}

json buildSolutionMatrix(const CBMSol& s, const CBMProblem& prob) {
    vector<vector<int>> matrix;
    matrix.reserve(prob.l);
    for (int row = 0; row < prob.l; ++row) {
        vector<int> rowVec;
        rowVec.reserve(prob.c);
        for (int col = 0; col < prob.c; ++col) rowVec.push_back(prob.binaryMatrix[row][s.sol[col]]);
        matrix.push_back(move(rowVec));
    }
    return matrix;
}

json buildInitialSolutions(const CBMProblem& prob) {
    json initSols = json::array();
    for (const auto& sol : prob.initialSolutions) {
        initSols.push_back({
            {"cost", sol.cost},
            {"solution", sol.sol},
            {"construction", constructionLabel(sol.construction)},
        });
    }
    return initSols;
}

void jsonOutput(const CBMSol& s, const CBMProblem& prob, Duration execution, Duration lkh) {
    json j = {
        {"final_solution",
         {
             {"cost", s.cost},
             {"solution", s.sol},
             {"matrix", buildSolutionMatrix(s, prob)},
         }},
        {"instance",
         {
             {"columns", prob.c},
             {"lines", prob.l},
         }},
        {"initial_solutions", buildInitialSolutions(prob)},
        {"execution_time_seconds", execution.count()},
        {"lkh_time_seconds", lkh.count()},
    };
    cout << j.dump(4) << "\n";
}

void iraceOutput(CBMSol& s, CBMProblem& prob, Duration execution) {
    cout << s.cost + (execution.count() / 10000.0) << "\n";
    prob.printS(s);
}

// ─── Debug helper ─────────────────────────────────────────────────────────────

void debugDeltaEval(CBMProblem& prob, int iterations) {
    for (int i = 0; i < iterations; ++i) {
        CBMSol s = prob.construction();
        prob.completeEval(s);

        CBMSol neighbor = prob.neighbor(s);
        prob.deltaEval(neighbor);
        const int deltaCost = neighbor.cost;

        prob.completeEval(neighbor);
        const int completeCost = neighbor.cost;

        if (deltaCost != completeCost) {
            ostringstream oss;
            oss << "Delta evaluation mismatch at iteration " << i << "\nInitial cost:  " << s.cost << "\nDelta cost:    " << deltaCost
                << "\nComplete cost: " << completeCost << "\nDifference:    " << (deltaCost - completeCost);
            throw runtime_error(oss.str());
        }
    }
}

// ─── main ─────────────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {
    const Config cfg = parseArgs(argc, argv);

    CBMProblem prob(cfg.filePath, cfg.movementType, cfg.constructionMethod, cfg.constructionBias, cfg.selectionBias, cfg.maxBlockSize, cfg.threads,
                    cfg.lkhS, cfg.lkhMaxTime, cfg.lkhCache);

    PT<CBMSol> algo(cfg.tempMin, cfg.tempMax, cfg.tempL, cfg.MKL, cfg.PTL, cfg.tempD, cfg.upType, max(cfg.PTL / cfg.tempUpdate, 1));

    const auto t0 = Clock::now();
    prob.createLKHInitialS();
    const auto t1 = Clock::now();
    CBMSol sol = algo.start(cfg.threads, &prob);
    const auto t2 = Clock::now();

    const Duration lkh = t1 - t0;
    const Duration execution = t2 - t0;

    if (cfg.irace)
        iraceOutput(sol, prob, execution);
    else
        jsonOutput(sol, prob, execution, lkh);

    return 0;
}