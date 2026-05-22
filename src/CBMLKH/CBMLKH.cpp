#include "CBMLKH.hpp"

CBMLKH::CBMLKH(int l, int c, int threads) : l(l), c(c) {
    this->diffMatrix.resize(c, vector<int>(c, 0));
    this->tspMatrix.resize(c + 1, vector<int>(c + 1, 0));
    this->onesToZeros.resize(c, vector<int>(c, 0));
    this->zerosToOnes.resize(c, vector<int>(c, 0));
    this->onesToOnes.resize(c, vector<int>(c, 0));
    this->threads = threads;

    this->computeMatrixes();
}

int CBMLKH::deltaEval(Solution& s) {
    auto insertionCost = [&](int L, int R, int e) -> int {
        if (L != -1 && R != -1) return (this->diffMatrix[L][e] + this->diffMatrix[e][R] - this->diffMatrix[L][R]) / 2;
        if (L == -1) return this->onesToZeros[e][R];
        return this->onesToZeros[e][L];
    };

    auto sharedOnes = [&](int i, int j) -> int { return (i == -1 || j == -1) ? 0 : this->onesToOnes[i][j]; };
    auto edgeDiff = [&](int i, int j) -> int { return (i != -1 && j != -1) ? this->diffMatrix[i][j] : 0; };

    switch (s.movement) {
        case REINSERTION: {
            int removedCost = insertionCost(s.mE[0], s.mE[1], s.mE[4]);
            int insertedCost = insertionCost(s.mE[2], s.mE[3], s.mE[4]);
            s.cost += -removedCost + insertedCost;
            break;
        }
        case SWAP: {
            int beforeA = insertionCost(s.mE[0], s.mE[1], s.mE[4]);
            int beforeB = insertionCost(s.mE[2], s.mE[3], s.mE[5]);
            int afterA = insertionCost(s.mE[2], s.mE[3], s.mE[4]);
            int afterB = insertionCost(s.mE[0], s.mE[1], s.mE[5]);
            s.cost += -(beforeA + beforeB) + (afterA + afterB);
            break;
        }
        case TWOOPT: {
            int before = sharedOnes(s.mE[1], s.mE[0]) + sharedOnes(s.mE[2], s.mE[3]);
            int after = sharedOnes(s.mE[2], s.mE[0]) + sharedOnes(s.mE[1], s.mE[3]);
            s.cost += before - after;
            break;
        }
    }
    return s.cost;
}

int CBMLKH::completeEval(Solution& s) {
    s.cost = 0;

    for (int row = 0; row < this->l; row++)
        if (this->binaryMatrix[row][s.sol[0]]) s.cost++;

    for (int col = 1; col < this->c; col++)
        for (int row = 0; row < this->l; row++)
            if (this->binaryMatrix[row][s.sol[col]] && !this->binaryMatrix[row][s.sol[col - 1]]) s.cost++;

    return s.cost;
}

void CBMLKH::computeMatrixes() {
#pragma omp parallel for num_threads(this->threads) collapse(2)
    for (int i = 0; i < c; i++) {
        for (int j = 0; j < c; j++) {
            if (i == j) continue;

            int o2z = 0, z2o = 0, o2o = 0;
            for (int row = 0; row < l; row++) {
                bool bi = binaryMatrix[row][i];
                bool bj = binaryMatrix[row][j];
                if (bi && !bj)
                    o2z++;
                else if (!bi && bj)
                    z2o++;
                else if (bi && bj)
                    o2o++;
            }
            this->onesToZeros[i][j] = o2z;
            this->zerosToOnes[i][j] = z2o;
            this->onesToOnes[i][j] = o2o;
            this->diffMatrix[i][j] = o2z + z2o;
            this->tspMatrix[i + 1][j + 1] = this->diffMatrix[i][j];
        }
    }

    this->tspMatrix[0][0] = 0;
    for (int i = 1; i <= this->c; i++) {
        int onesCount = 0;
        for (int row = 0; row < this->l; row++)
            if (binaryMatrix[row][i - 1]) onesCount++;
        this->tspMatrix[0][i] = onesCount;
        this->tspMatrix[i][0] = onesCount;
    }
}