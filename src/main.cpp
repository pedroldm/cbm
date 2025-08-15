#include <iostream>
#include <sstream>
#include <string>

#include "CBMProblem.hpp"
#include "PTAPI/include/PT.h"

int main(int argc, char* argv[]) {
    float tempMin = 0.3f;
    float tempMax = 2.0f;
    double constructionBias = 2.5;
    int tempL = 12;
    float MKL = 400;
    int PTL = 2000;
    int tempD = 4;
    int upType = 1;
    int tempUpdate = 3;
    int movementType = 3;
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
            istringstream(arg.substr(20)) >> tempMax;
        } else if (arg.find("--movementType=") == 0) {
            istringstream(arg.substr(15)) >> movementType;
        } else if (arg.find("--filePath=") == 0) {
            istringstream(arg.substr(11)) >> filePath;
        } else {
            throw runtime_error("Unkown argument: " + arg);
        }
    }

    CBMProblem* prob = new CBMProblem(filePath, movementType, constructionBias);
    PT<CBMSol> algo(tempMin, tempMax, tempL, MKL, PTL, tempD, upType, max(PTL / tempUpdate, 1));
    
    CBMSol sol = prob->construction();
    int cost = prob->evaluate(sol);
    
    CBMSol newS = prob->neighbor(sol);

    int newCost = prob->evaluate(newS);
    newS.cost = cost;
    int deltaCost = prob->deltaEvaluate(newS);
    //CBMSol sol = algo.start(1, prob);
    //cout << sol << endl;
    //prob->printMatrix();
    //cout << endl << "# ------------- #" << endl;
    //prob->printMatrix(&sol);
    delete prob;

    return 0;
}