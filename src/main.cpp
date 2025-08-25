#include <iostream>
#include <sstream>
#include <string>

#include "CBMProblem.hpp"


void validateDeltaF(CBMProblem* prob);

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
    int movementType = 2;
    int threads = 4;
    string filePath;

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
        } else if (arg.find("--threads=") == 0) {
            istringstream(arg.substr(10)) >> threads;
        } else if (arg.find("--movementType=") == 0) {
            istringstream(arg.substr(15)) >> movementType;
        } else if (arg.find("--constructionBias=") == 0) {
            istringstream(arg.substr(19)) >> constructionBias;
        } else if (arg.find("--filePath=") == 0) {
            istringstream(arg.substr(11)) >> filePath;
        } else {
            throw runtime_error("Unkown argument: " + arg);
        }
    }

    CBMProblem* prob = new CBMProblem(filePath, movementType, constructionBias);
    CBMProblem->run();

    delete prob;

    return 0;
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