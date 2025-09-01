#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <iomanip>

#include "CBMProblem.hpp"
#include "PTAPI/include/PT.h"
#include "IO/json.hpp"

using json = nlohmann::json;

void validateDeltaF(CBMProblem* prob);
void jsonOutput(CBMSol& s, CBMProblem& prob, PT<CBMSol>& algo);

int main(int argc, char* argv[]) {
    float tempMin = 0.05f;
    float tempMax = 2.0f;
    double constructionBias = 2.5;
    int tempL = 4;
    float MKL = 400;
    int PTL = 2000;
    int tempD = 4;
    int upType = 1;
    int tempUpdate = 3;
    int movementType = 4;
    int threads = 1;
    int maxBlockSize = 3;
    bool irace = false;
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
        } else if (arg.find("--tempD=") == 0) {
            istringstream(arg.substr(8)) >> tempD;
        } else if (arg.find("--upType=") == 0) {
            istringstream(arg.substr(9)) >> upType;
        } else if (arg.find("--tempUpdate=") == 0) {
            istringstream(arg.substr(13)) >> tempUpdate;
        } else if (arg.find("--tempMax=") == 0) {
            istringstream(arg.substr(10)) >> tempMax;
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
                throw invalid_argument("Invalid value for --irace (expected true/false)");
        } else if (arg.find("--filePath=") == 0) {
            istringstream(arg.substr(11)) >> filePath;
        } else {
            throw runtime_error("Unkown argument: " + arg);
        }
    }

    CBMProblem* prob = new CBMProblem(filePath, movementType, constructionBias, maxBlockSize, threads);
    PT<CBMSol> algo(tempMin, tempMax, tempL, MKL, PTL, tempD, upType, max(PTL / tempUpdate, 1));
    auto start = chrono::high_resolution_clock::now();
    CBMSol sol = algo.start(threads, prob);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;

    if(irace) {
        cout << sol.cost + (elapsed.count() / 10000.0) << endl;
    }
    else
        jsonOutput(sol, *prob, algo);

    delete prob;
    return 0;
}

void jsonOutput(CBMSol& s, CBMProblem& prob, PT<CBMSol>& algo) {
    json j;

    // Final solution
    j["final_solution"]["cost"] = s.cost;
    j["final_solution"]["solution"] = s.sol;

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
        js["construction"] = (sol.construction == GREEDY ? "GREEDY" : "LKH");
        initSols.push_back(js);
    }
    j["initial_solutions"] = initSols;

    cout << j.dump(4) << endl;
}

void validateDeltaF(CBMProblem* prob) {
    for(int i = 0 ; i < 1000 ; i++) {
        CBMSol sol = prob->construction();
        int cost = prob->evaluate(sol);
        CBMSol newS = prob->neighbor(sol);
        newS.cost = 0;
        int newCost = prob->evaluate(newS);
        cout << "\nNovo Custo (n²): " << newCost << endl;
        newS.cost = cost;
        int deltaCost = prob->deltaEval(newS);
        cout << "Novo Custo (1): " << deltaCost << endl;
        if(newCost != deltaCost) {
            newS.cost = cost;
            int deltaCost = prob->deltaEval(newS);
            throw runtime_error("Unmatch: " + to_string(newCost) + " != " + to_string(deltaCost));
        } else {
            cout << newCost << " = " << deltaCost << endl;;
        }
    }
}