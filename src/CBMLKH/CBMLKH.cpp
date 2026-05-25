#include "CBMLKH.hpp"

#include <filesystem>
#include <iostream>
#include <sstream>

namespace fs = std::filesystem;

CBMLKH::CBMLKH(const Config& cfg)
    : filePath(cfg.instancePath),
      threads(cfg.threads),
      l(0),
      c(0),
      rng_device(),
      mersenne_engine(rng_device()),
      constructionBias(cfg.constructionBias) {
    ifstream input(filePath);
    if (!input.is_open()) throw runtime_error("Error opening file: " + filePath);

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
    this->threads = threads;

    auto pos = filePath.find_last_of("/\\");
    this->instanceName = (pos == string::npos) ? filePath : filePath.substr(pos + 1);
    this->currPath = fs::current_path();
    this->lkhPath = this->currPath / "src" / "LKH3" / "LKH";
    this->lkhMaxTime = cfg.lkhMaxTime;
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

void CBMLKH::toTSP(const string& suffix, const vector<int>& subset) {
    fs::path tspPath = fs::path("/tmp") / (this->instanceName + suffix + ".tsp");
    fs::path parPath = fs::path("/tmp") / (this->instanceName + suffix + ".par");

    ofstream tsp(tspPath), par(parPath);
    if (!tsp || !par) throw runtime_error("Error creating TSP/PAR files");

    int dim = subset.empty() ? this->c : static_cast<int>(subset.size());

    tsp << "NAME : " << this->instanceName << "\n"
        << "TYPE : TSP\n"
        << "DIMENSION : " << dim + 1 << "\n"
        << "EDGE_WEIGHT_TYPE : EXPLICIT\n"
        << "EDGE_WEIGHT_FORMAT : FULL_MATRIX\n"
        << "EDGE_WEIGHT_SECTION\n";

    for (int i = 0; i <= dim; ++i) {
        for (int j = 0; j <= dim; ++j) {
            int val = 0;
            int orig_i = (i == 0) ? -1 : (subset.empty() ? (i - 1) : subset[i - 1]);
            int orig_j = (j == 0) ? -1 : (subset.empty() ? (j - 1) : subset[j - 1]);
            if (orig_i == -1 && orig_j == -1)
                val = 0;
            else if (orig_i == -1)
                val = this->onesSum[orig_j];
            else if (orig_j == -1)
                val = this->onesSum[orig_i];
            else
                val = this->diffMatrix[orig_i][orig_j];
            tsp << val << (j == dim ? "\n" : " ");
        }
    }
    tsp << "EOF\n";

    string solFile = (fs::path("/tmp") / (this->instanceName + suffix + ".sol")).string();
    string initialTour = (fs::path("/tmp") / (this->instanceName + suffix + "_initial_.tour")).string();

    par << "PROBLEM_FILE = " << tspPath.string() << "\n"
        << "OUTPUT_TOUR_FILE = " << solFile << "\n"
        << "RUNS = 1\n"
        << "INITIAL_TOUR_FILE = " << initialTour << "\n"
        << "TIME_LIMIT = " << this->lkhMaxTime << "\n";
}

void CBMLKH::initialTour(const string& suffix, const vector<int>& subset) {
    fs::path tourPath = fs::path("/tmp") / (this->instanceName + suffix + "_initial_.tour");
    ofstream tour(tourPath);
    if (!tour) throw runtime_error("Error creating TOUR file");

    int dim = subset.empty() ? this->c : static_cast<int>(subset.size());

    tour << "TOUR_SECTION\n"
         << "1\n";
    for (int i = 0; i < dim; i++) tour << (i + 2) << "\n";
    tour << "-1\nEOF";
}

void CBMLKH::runLKH(string suffix) {
    fs::path parPath = fs::path("/tmp") / (this->instanceName + suffix + ".par");

    fs::path lkhExec = fs::path(LKH_EXEC_PATH);
    if (!fs::exists(lkhExec)) throw runtime_error("LKH executable not found: " + lkhExec.string());
    if (!fs::exists(parPath)) throw runtime_error("PAR file not found: " + parPath.string());

    string command = lkhExec.string() + " " + parPath.string() + " > /dev/null 2>&1";
    int ret = system(command.c_str());
    if (ret != 0) cerr << "LKH returned error code: " << ret << endl;
}

vector<int> CBMLKH::fromTSP(const string& suffix, const vector<int>& subset) {
    fs::path solPath = fs::path("/tmp") / (this->instanceName + suffix + ".sol");

    ifstream sol(solPath);
    if (!sol) throw runtime_error("Error opening solution file: " + solPath.string());

    vector<int> tour;
    string line;
    bool inTourSection = false;

    while (getline(sol, line)) {
        if (!inTourSection) {
            if (line == "TOUR_SECTION") inTourSection = true;
            continue;
        }
        istringstream iss(line);
        int node;
        while (iss >> node) {
            if (node == -1) {
                inTourSection = false;
                break;
            }
            int idx = node - 2;
            if (idx == -1)
                tour.push_back(-1);
            else if (subset.empty())
                tour.push_back(idx);
            else
                tour.push_back(subset[idx]);
        }
    }

    auto depotIt = find(tour.begin(), tour.end(), -1);
    if (depotIt != tour.end()) {
        rotate(tour.begin(), depotIt, tour.end());
        tour.erase(tour.begin());
    }

    return tour;
}

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
    s.cost = this->tspMatrix[0][s.sol[0] + 1];

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
        this->onesSum[i - 1] = onesCount;
    }
}

void CBMLKH::countBlocksPerColumn(Solution& s) {
    s.blocksCount[0] = this->onesSum[s.sol[0]];
    for (int i = 1; i < this->c; i++) s.blocksCount[i] = this->zerosToOnes[s.sol[i - 1]][s.sol[i]];
}