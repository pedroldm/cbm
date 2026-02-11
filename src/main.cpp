#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "CBMProblem.hpp"
#include "IO/json.hpp"
#include "PTAPI/include/PT.h"

using json = nlohmann::json;

void jsonOutput(CBMSol& s, CBMProblem& prob, PT<CBMSol>& algo,
                chrono::duration<double> execution,
                chrono::duration<double> lkh);
void debugDeltaEval(CBMProblem& prob, PT<CBMSol>& algo, int iterations);

int main(int argc, char* argv[]) {
    float tempMin = 0.05f;
    float tempMax = 2.0f;
    double constructionBias = 2.5;
    int constructionMethod = 1;
    double selectionBias = 1.0;
    int tempL = 4;
    int lkhS = 0;
    float MKL = 400;
    int PTL = 2000;
    int tempD = 4;
    int upType = 1;
    int tempUpdate = 3;
    int movementType = 4;
    int threads = 1;
    int lkhMaxTime = 120;
    int maxBlockSize = 3;
    bool irace = false;
    bool lkhCache = false;
    string filePath;
    string tspPath = "./instances/tsp/";

    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg.find("--tempMin=") == 0) {
            istringstream(arg.substr(10)) >> tempMin;
        } else if (arg.find("--tempL=") == 0) {
            istringstream(arg.substr(8)) >> tempL;
        } else if (arg.find("--MKL=") == 0) {
            istringstream(arg.substr(6)) >> MKL;
        } else if (arg.find("--PTL=") == 0) {
            istringstream(arg.substr(6)) >> PTL;
        } else if (arg.find("--selectionBias=") == 0) {
            istringstream(arg.substr(16)) >> selectionBias;
        } else if (arg.find("--constructionMethod=") == 0) {
            istringstream(arg.substr(21)) >> constructionMethod;
        } else if (arg.find("--tempD=") == 0) {
            istringstream(arg.substr(8)) >> tempD;
        } else if (arg.find("--upType=") == 0) {
            istringstream(arg.substr(9)) >> upType;
        } else if (arg.find("--tempUpdate=") == 0) {
            istringstream(arg.substr(13)) >> tempUpdate;
        } else if (arg.find("--tempMax=") == 0) {
            istringstream(arg.substr(10)) >> tempMax;
        } else if (arg.find("--lkhMaxTime=") == 0) {
            istringstream(arg.substr(13)) >> lkhMaxTime;
        } else if (arg.find("--lkhS=") == 0) {
            istringstream(arg.substr(7)) >> lkhS;
        } else if (arg.find("--maxBlockSize=") == 0) {
            istringstream(arg.substr(15)) >> maxBlockSize;
        } else if (arg.find("--threads=") == 0) {
            istringstream(arg.substr(10)) >> threads;
        } else if (arg.find("--movementType=") == 0) {
            istringstream(arg.substr(15)) >> movementType;
        } else if (arg.find("--constructionBias=") == 0) {
            istringstream(arg.substr(19)) >> constructionBias;
        } else if (arg.find("--irace=") == 0) {
            string value = arg.substr(8);
            if (value == "true" || value == "1")
                irace = true;
            else if (value == "false" || value == "0")
                irace = false;
            else
                throw invalid_argument(
                    "Invalid value for --irace (expected true/false)");
        } else if (arg.find("--lkhCache=") == 0) {
            string value = arg.substr(11);
            if (value == "true" || value == "1")
                lkhCache = true;
            else if (value == "false" || value == "0")
                lkhCache = false;
            else
                throw invalid_argument(
                    "Invalid value for --lkhCache (expected true/false)");
        } else if (arg.find("--filePath=") == 0) {
            istringstream(arg.substr(11)) >> filePath;
        } else {
            throw runtime_error("Unkown argument: " + arg);
        }
    }

    CBMProblem* prob =
        new CBMProblem(filePath, movementType, constructionMethod, constructionBias, selectionBias,
                       maxBlockSize, threads, lkhS, lkhMaxTime, lkhCache);
    PT<CBMSol> algo(tempMin, tempMax, tempL, MKL, PTL, tempD, upType,
                    max(PTL / tempUpdate, 1));

    auto start = chrono::high_resolution_clock::now();
    prob->createLKHInitialS();
    auto lkhFinish = chrono::high_resolution_clock::now();
    CBMSol sol = algo.start(threads, prob);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> execution = end - start;
    chrono::duration<double> lkh = lkhFinish - start;

    if (irace) {
        cout << sol.cost + (execution.count() / 10000.0) << endl;
        prob->printS(sol);
    } else
        jsonOutput(sol, *prob, algo, execution, lkh);

    /*
    debugDeltaEval(*prob, algo, 100000);
    */

    delete prob;
    return 0;
}

void debugDeltaEval(CBMProblem& prob, PT<CBMSol>& algo, int iterations) {
    for (int i = 0 ; i < iterations ; i++) {
        CBMSol s = prob.construction();
        prob.completeEval(s);

        cout << "Iteration: " << i << endl;
        cout << "Initial Cost: " << s.cost << endl;

        CBMSol neighbor = prob.neighbor(s);

        prob.deltaEval(neighbor);
        int deltaCost = neighbor.cost;
        cout << "Movement: " << neighbor.movement << endl;
        cout << "Delta Cost: " << deltaCost << endl;

        prob.completeEval(neighbor);
        int completeCost = neighbor.cost;
        cout << "Complete Cost: " << completeCost << endl;

        if (deltaCost != completeCost) {
            std::ostringstream oss;
            oss << "Delta evaluation mismatch at iteration " << i
                << "\nInitial cost: " << s.cost
                << "\nDelta cost: " << deltaCost
                << "\nComplete cost: " << completeCost
                << "\nDifference: " << (deltaCost - completeCost);

            throw std::runtime_error(oss.str());
        }
    }
}


void jsonOutput(CBMSol& s, CBMProblem& prob, PT<CBMSol>& algo,
                chrono::duration<double> execution,
                chrono::duration<double> lkh) {
    json j;

    // Final solution
    j["final_solution"]["cost"] = s.cost;
    j["final_solution"]["solution"] = s.sol;

    j["instance"]["columns"] = prob.c;
    j["instance"]["lines"] = prob.l;

    vector<vector<int>> matrix;
    for (int row = 0; row < prob.l; row++) {
        vector<int> rowVec;
        for (int col = 0; col < prob.c; col++) {
            rowVec.push_back(prob.binaryMatrix[row][s.sol[col]]);
        }
        matrix.push_back(rowVec);
    }
    j["final_solution"]["matrix"] = matrix;

    // Initial solutions
    json initSols = json::array();
    for (const auto& sol : prob.initialSolutions) {
        json js;
        js["cost"] = sol.cost;
        js["solution"] = sol.sol;

        if (sol.construction == GREEDY)
            js["construction"] = "GREEDY";
        else if (sol.construction == ONEBLOCK)
            js["construction"] = "ONEBLOCK";
        else
            js["construction"] = "LKH";

        initSols.push_back(js);
    }
    j["initial_solutions"] = initSols;

    // Times
    j["execution_time_seconds"] = execution.count();
    j["lkh_time_seconds"] = lkh.count();

    cout << j.dump(4) << endl;
}