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

    input >> this->l;
    input >> this->c;

    this->binaryMatrix.resize(this->l, std::vector<bool>(this->c, false));

    int n, e;
    for(int i = 0 ; i < this->l ; i++) {
        input >> n;
        for(int j = 0 ; j < n ; j++) {
            input >> e;
            this->binaryMatrix[i][e - 1] = true;
        }
    }

    this->movementType = movementType;
    this->maxTemp = this->l * maxTempProportion;
}

CBMSol CBMProblem::construction() {
    CBMSol ss;
    ss.cost = 0;
    ss.sol.resize(this->c);
    std::iota(ss.sol.begin(), ss.sol.end(), 0);
    std::shuffle(ss.sol.begin(), ss.sol.end(), this->mersenne_engine);
    return ss;
}

CBMSol CBMProblem::neighbor(CBMSol s) {
    std::uniform_int_distribution<> dist(0, this->c - 1);
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
    s.cost = 0;

    for (int row = 0; row < this->l; row++) {
        if (this->binaryMatrix[row][s.sol[0]]) {
            s.cost++;
        }
    }
    for (int col = 1 ; col < this->c ; col++) {
        for(int row = 0 ; row < this->l ; row++) {
            if(this->binaryMatrix[row][s.sol[col]] && !this->binaryMatrix[row][s.sol[col - 1]])
                s.cost++;
        }
    }

    return s.cost;
}

void CBMProblem::printMatrix(const CBMSol* s) {
    if(!s) {
           for (int row = 0; row < this->l; row++) {
            for (int col = 0; col < this->c; col++) {
                std::cout << this->binaryMatrix[row][col] << " ";
            }
            std::cout << "\n";
        }
    }
    else {
        for(int row = 0 ; row < this->l ; row ++) {
            for (int col = 0; col < this->c; col++) {
                std::cout << this->binaryMatrix[row][s->sol[col]] << " ";
            }
            std::cout << "\n";
        }
    }
}