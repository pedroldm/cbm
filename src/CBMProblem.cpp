#include "CBMProblem.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

CBMProblem::CBMProblem(std::string filename, int movementType,
                       double maxTempProportion)
    : mersenne_engine(rng_device()) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    input >> this->n;
    input >> this->m;

    this->binaryMatrix.resize(m, std::vector<bool>(n + 1, false));

    for (int block = 0; block < m; ++block) {
        int k;
        input >> k;
        for (int x = 0; x < k; ++x) {
            int elementId;
            input >> elementId;
            this->binaryMatrix[block][elementId] = true;
        }
    }

    this->movementType = movementType;
    this->maxTemp = this->n * maxTempProportion;
}

CBMSol CBMProblem::construction() {
    CBMSol ss;
    ss.cost = 0;
    ss.sol.resize(n);
    std::iota(ss.sol.begin(), ss.sol.end(), 0);
    std::shuffle(ss.sol.begin(), ss.sol.end(), this->mersenne_engine);
    return ss;
}

CBMSol CBMProblem::neighbor(CBMSol s) {
    std::uniform_int_distribution<> dist(0, this->n - 1);
    int index = dist(this->mersenne_engine);
    int newIndex = dist(this->mersenne_engine);

    switch (this->movementType) {
        case 1: /* swap */
            std::swap(s.sol[index], s.sol[newIndex]);
            break;
        case 2: /* 2-opt */
            if (index > newIndex) std::swap(index, newIndex);
            std::reverse(s.sol.begin() + index, s.sol.begin() + newIndex + 1);
            break;
        case 3: /* re-insertion */
            if (index < newIndex) newIndex -= 1;
            int element = s.sol[index];
            s.sol.erase(s.sol.begin() + index);
            s.sol.insert(s.sol.begin() + newIndex, element);
    }

    return s;
}

int CBMProblem::evaluate(CBMSol& s) {

}