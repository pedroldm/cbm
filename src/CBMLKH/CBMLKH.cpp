#include "CBMLKH.hpp"

#include <filesystem>
#include <iostream>
#include <sstream>

namespace fs = std::filesystem;

CBMLKH::CBMLKH(const Config& cfg) : cfg(cfg), validator(cfg.instancePath), l(0), c(0), rng_device(), mersenne_engine(rng_device()) {
    ifstream input(cfg.instancePath);
    if (!input.is_open()) throw runtime_error("Error opening file: " + cfg.instancePath);

    input >> this->l >> this->c;

    this->binaryMatrix.resize(this->l);
    for (int row = 0; row < this->l; row++) {
        this->binaryMatrix[row].resize(this->c, 0);
    }
    this->diffMatrix.resize(c, vector<int>(c, 0));
    this->tspMatrix.resize(c + 1, vector<int>(c + 1, 0));
    this->onesToZeros.resize(c, vector<int>(c, 0));
    this->zerosToOnes.resize(c, vector<int>(c, 0));
    this->onesToOnes.resize(c, vector<int>(c, 0));
    this->onesSum.resize(c, 0);

    auto pos = cfg.instancePath.find_last_of("/\\");
    this->instanceName = (pos == string::npos) ? cfg.instancePath : cfg.instancePath.substr(pos + 1);
    this->currPath = fs::current_path();
    this->lkhPath = this->currPath / "src" / "LKH3" / "LKH";
    this->lkhCache = false;

    int nonZeroCount, col;
    for (int row = 0; row < this->l; row++) {
        input >> nonZeroCount;
        for (int j = 0; j < nonZeroCount; j++) {
            input >> col;
            this->binaryMatrix[row][col - 1] = true;
        }
    }

    this->computeMatrixes();
}

void CBMLKH::run() {
    vector<Solution> initialSolutions(this->cfg.threads);

#pragma omp parallel for num_threads(this->cfg.threads)
    for (int i = 0; i < this->cfg.threads; i++) {
        initialSolutions[i] = this->greedyConstruction();
        Trajectory t = this->LKHILS(initialSolutions[i]);
    }
}

Trajectory CBMLKH::LKHILS(Solution& initial) {
    Trajectory trajectory(initial);

    auto start = chrono::steady_clock::now();
    auto deadline = start + chrono::seconds(cfg.lkhMaxTime);

    for (int i = 0; i < cfg.maxIterations; i++) {
        auto now = chrono::steady_clock::now();
        if (now >= deadline) break;
        Solution neighbor = this->ILSNeighbor(trajectory.currentSolution);
        if (neighbor.cost < trajectory.bestSolution.cost) {
            validator.validate(neighbor.sol, neighbor.cost);
            trajectory.record(neighbor, i, chrono::duration_cast<chrono::milliseconds>(now - start).count());
        }
    }

    return trajectory;
}

Solution CBMLKH::ILSNeighbor(const Solution s) {}

int CBMLKH::deltaEval(Solution& s) {
    auto insertionCost = [&](int L, int R, int e) -> int {
        if (L != -1 && R != -1) return (this->diffMatrix[L][e] + this->diffMatrix[e][R] - this->diffMatrix[L][R]) / 2;
        if (L == -1) return this->onesToZeros[e][R];
        return this->onesToZeros[e][L];
    };

    auto sharedOnes = [&](int i, int j) -> int { return (i == -1 || j == -1) ? 0 : this->onesToOnes[i][j]; };

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
    s.cost = this->onesSum[s.sol[0]];

    for (int col = 1; col < this->c; col++) {
        s.cost += this->zerosToOnes[s.sol[col - 1]][s.sol[col]];
    }

    return s.cost;
}

Solution CBMLKH::greedyConstruction() {
    uniform_int_distribution<> colDist(0, this->c - 1);
    unordered_set<int> remaining;
    Solution s;

    s.cost = 0;
    s.sol.resize(this->c);
    s.blocksCount.resize(this->c);
    int current = colDist(this->mersenne_engine);

    for (int i = 0; i < this->c; i++) remaining.insert(i);

    int pos = 0;
    remaining.erase(current);
    s.sol[pos++] = current;

    while (!remaining.empty()) {
        current = this->nextInsertion(current, remaining);
        s.sol[pos++] = current;
        remaining.erase(current);
    }

    return s;
}

int CBMLKH::nextInsertion(int& current, unordered_set<int>& remaining) {
    vector<tuple<int, int>> candidates;
    candidates.reserve(remaining.size());
    for (int i : remaining) candidates.push_back({this->l - this->diffMatrix[current][i], i});

    sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return get<0>(a) > get<0>(b); });

    vector<double> weights(candidates.size());
    for (size_t rank = 0; rank < candidates.size(); ++rank) weights[rank] = 1.0 / pow(rank + 1, this->constructionBias);

    double totalWeight = accumulate(weights.begin(), weights.end(), 0.0);
    uniform_real_distribution<double> dist(0.0, totalWeight);
    double roll = dist(this->mersenne_engine);
    double cumulative = 0.0;
    size_t chosen = 0;

    for (; chosen < weights.size(); ++chosen) {
        cumulative += weights[chosen];
        if (roll < cumulative) break;
    }

    return get<1>(candidates[chosen]);
}

void CBMLKH::computeMatrixes() {
#pragma omp parallel for num_threads(this->cfg.threads) collapse(2)
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
        this->onesSum[i - 1] = onesCount;
    }
}

void CBMLKH::countBlocksPerColumn(Solution& s) {
    s.blocksCount[0] = this->onesSum[s.sol[0]];
    for (int i = 1; i < this->c; i++) s.blocksCount[i] = this->zerosToOnes[s.sol[i - 1]][s.sol[i]];
}

vector<DenseSegment> CBMLKH::findDenseSegments(Solution s, int minSize, int maxSize, double minScore) {
    vector<DenseSegment> segments;
    vector<int> prefix(s.blocksCount.size() + 1, 0);

    for (int i = 0; i < (int)s.blocksCount.size(); i++) prefix[i + 1] = prefix[i] + s.blocksCount[i];
    for (int l = 0; l < (int)s.blocksCount.size(); l++) {
        for (int r = l + minSize - 1; r < (int)s.blocksCount.size() && (r - l + 1) <= maxSize; r++) {
            int size = r - l + 1;
            int sum = prefix[r + 1] - prefix[l];
            double avg = (double)sum / size;
            double score = avg * size;
            if (score >= minScore) {
                segments.push_back({l, r, sum, avg, score});
            }
        }
    }

    sort(segments.begin(), segments.end(), [](const auto& a, const auto& b) { return a.score > b.score; });

    return segments;
}