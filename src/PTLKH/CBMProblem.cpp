#include "CBMProblem.hpp"

// ─────────────────────────────────────────────
//  Constructor
// ─────────────────────────────────────────────

CBMProblem::CBMProblem(string filename, int movementType, int constructionMethod, double constructionBias, double selectionBias, int maxBlockSize,
                       int threads, int lkhS, int lkhMaxTime, bool lkhCache)
    : mersenne_engine(rng_device()),
      constructionMethod(constructionMethod),
      constructionBias(constructionBias),
      selectionBias(selectionBias),
      maxBlockSize(maxBlockSize),
      threads(threads),
      ISCount(lkhS),
      lkhMaxTime(lkhMaxTime),
      lkhCache(lkhCache) {
    ifstream input(filename);
    if (!input.is_open()) throw runtime_error("Error opening file: " + filename);

    auto pos = filename.find_last_of("/\\");
    this->lkhS = lkhS;
    this->instanceName = (pos == string::npos) ? filename : filename.substr(pos + 1);
    this->currPath = filesystem::current_path();
    this->lkhPath = this->currPath / "src" / "LKH3" / "LKH";

    input >> this->l >> this->c;

    // Allocate all matrices and vectors up-front
    this->binaryMatrix.resize(this->l);
    this->diffMatrix.resize(this->c, vector<int>(this->c, 0));
    this->tspMatrix.resize(this->c + 1, vector<int>(this->c + 1, 0));
    this->onesToZeros.resize(this->c, vector<int>(this->c, 0));
    this->zerosToOnes.resize(this->c, vector<int>(this->c, 0));
    this->onesToOnes.resize(this->c, vector<int>(this->c, 0));
    this->diffVec.resize(this->c);
    this->selectionWeights.resize(this->c);

    // Read sparse binary matrix (rows given as lists of 1-indexed column indices)
    int nonZeroCount, col;
    for (int row = 0; row < this->l; row++) {
        input >> nonZeroCount;
        for (int j = 0; j < nonZeroCount; j++) {
            input >> col;
            this->binaryMatrix[row][col - 1] = true;
        }
    }

    this->computeMatrixes();
    this->movementType = movementType;
}

// ─────────────────────────────────────────────
//  Selection
// ─────────────────────────────────────────────

// Picks a column index using weighted-random selection (higher diff → more likely).
int CBMProblem::biasedSelection() {
    uniform_real_distribution<double> dist(0.0, this->totalWeight);
    double roll = dist(this->mersenne_engine);

    double cumulative = 0.0;
    size_t idx = 0;
    for (; idx < this->c; ++idx) {
        cumulative += this->selectionWeights[idx];
        if (roll < cumulative) break;
    }

    auto [score, column] = this->diffVec[idx];
    return column;
}

// ─────────────────────────────────────────────
//  Neighbor generation
// ─────────────────────────────────────────────

CBMSol CBMProblem::neighbor(CBMSol s) {
    uniform_int_distribution<> colDist(0, this->c - 1);
    uniform_int_distribution<> rowDist(0, this->l - 1);
    uniform_int_distribution<> moveDist(1, 3);

    // Pick the primary column (biased or uniform)
    int index = (this->selectionBias > 1) ? this->biasedSelection() : colDist(this->mersenne_engine);
    int newIndex = colDist(this->mersenne_engine);

    // Helpers to safely read adjacent elements in the permutation
    auto leftOf = [&](int i) { return (i > 0) ? s.sol[i - 1] : -1; };
    auto rightOf = [&](int i) { return (i + 1 < this->c) ? s.sol[i + 1] : -1; };

    int movement = (this->movementType < 4) ? this->movementType : moveDist(this->mersenne_engine);

    switch (movement) {
        // ── Case 1: swap two columns (adjacent → treated as 2-opt, non-adjacent → SWAP)
        case 1:
            if (index > newIndex) swap(index, newIndex);
            if (abs(newIndex - index) > 1) {
                s.movement = SWAP;
                s.mE = {leftOf(index), rightOf(index), leftOf(newIndex), rightOf(newIndex), s.sol[index], s.sol[newIndex]};
                swap(s.sol[index], s.sol[newIndex]);
            } else {
                s.movement = TWOOPT;
                s.mE = {leftOf(index), s.sol[index], s.sol[newIndex], rightOf(newIndex)};
                swap(s.sol[index], s.sol[newIndex]);
            }
            break;

        // ── Case 2: reverse a segment (classic 2-opt)
        case 2:
            if (index > newIndex) swap(index, newIndex);
            s.movement = TWOOPT;
            s.mE = {leftOf(index), s.sol[index], s.sol[newIndex], rightOf(newIndex)};
            reverse(s.sol.begin() + index, s.sol.begin() + newIndex + 1);
            break;

        // ── Case 3: remove a column and reinsert it at another position
        case 3: {
            s.movement = REINSERTION;
            int elem = s.sol[index];
            int origLeft = leftOf(index);
            int origRight = rightOf(index);
            s.sol.erase(s.sol.begin() + index);
            s.sol.insert(s.sol.begin() + newIndex, elem);
            int destLeft = leftOf(newIndex);
            int destRight = rightOf(newIndex);
            s.mE = {origLeft, origRight, destLeft, destRight, elem};
            break;
        }

        // ── Case 4: move a contiguous "one-block" to another position
        case 4:
            this->oneBlockMovement(s);
            s.movement = ONEBLOCKM;
            break;
    }

    return s;
}

// ─────────────────────────────────────────────
//  Evaluation
// ─────────────────────────────────────────────

int CBMProblem::evaluate(CBMSol& s) {
    // Use fast incremental eval when cost is already known, full scan otherwise
    return s.cost ? this->deltaEval(s) : this->completeEval(s);
}

// Full O(l·c) evaluation: count 0→1 transitions across all rows.
int CBMProblem::completeEval(CBMSol& s) {
    s.cost = 0;

    // First column contributes a block-start for every row where it is 1
    for (int row = 0; row < this->l; row++)
        if (this->binaryMatrix[row][s.sol[0]]) s.cost++;

    // Subsequent columns: count transitions from 0 to 1
    for (int col = 1; col < this->c; col++)
        for (int row = 0; row < this->l; row++)
            if (this->binaryMatrix[row][s.sol[col]] && !this->binaryMatrix[row][s.sol[col - 1]]) s.cost++;

    return s.cost;
}

// Incremental evaluation: only recomputes cost for the edges touched by the last move.
int CBMProblem::deltaEval(CBMSol& s) {
    // Cost of inserting `e` between neighbours L and R
    auto insertionCost = [&](int L, int R, int e) -> int {
        if (L != -1 && R != -1) return (this->diffMatrix[L][e] + this->diffMatrix[e][R] - this->diffMatrix[L][R]) / 2;
        if (L == -1) return this->onesToZeros[e][R];
        /* R == -1 */ return this->onesToZeros[e][L];
    };

    // Shared-ones count for adjacent pair (used by 2-opt)
    auto sharedOnes = [&](int i, int j) -> int { return (i == -1 || j == -1) ? 0 : this->onesToOnes[i][j]; };

    // Diff between two columns (0 when either is a boundary sentinel)
    auto edgeDiff = [&](int i, int j) -> int { return (i != -1 && j != -1) ? this->diffMatrix[i][j] : 0; };

    switch (s.movement) {
        case REINSERTION: {
            // mE = { origLeft, origRight, destLeft, destRight, element }
            int removedCost = insertionCost(s.mE[0], s.mE[1], s.mE[4]);
            int insertedCost = insertionCost(s.mE[2], s.mE[3], s.mE[4]);
            s.cost += -removedCost + insertedCost;
            break;
        }
        case SWAP: {
            // mE = { leftA, rightA, leftB, rightB, elemA, elemB }
            int beforeA = insertionCost(s.mE[0], s.mE[1], s.mE[4]);
            int beforeB = insertionCost(s.mE[2], s.mE[3], s.mE[5]);
            int afterA = insertionCost(s.mE[2], s.mE[3], s.mE[4]);
            int afterB = insertionCost(s.mE[0], s.mE[1], s.mE[5]);
            s.cost += -(beforeA + beforeB) + (afterA + afterB);
            break;
        }
        case TWOOPT: {
            // mE = { leftEdge, elemA, elemB, rightEdge }
            // Reversing the segment changes which pairs share ones at the boundary edges
            int before = sharedOnes(s.mE[1], s.mE[0]) + sharedOnes(s.mE[2], s.mE[3]);
            int after = sharedOnes(s.mE[2], s.mE[0]) + sharedOnes(s.mE[1], s.mE[3]);
            s.cost += before - after;
            break;
        }
    }
    return s.cost;
}

// ─────────────────────────────────────────────
//  Matrix pre-computation
// ─────────────────────────────────────────────

void CBMProblem::computeMatrixes() {
// Fill pairwise relationship counts between every pair of columns
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

    // Depot node (index 0) edge weights = number of ones in each column
    this->tspMatrix[0][0] = 0;
    for (int i = 1; i <= this->c; i++) {
        int onesCount = 0;
        for (int row = 0; row < this->l; row++)
            if (binaryMatrix[row][i - 1]) onesCount++;
        this->tspMatrix[0][i] = onesCount;
        this->tspMatrix[i][0] = onesCount;
    }

// Build biased-selection weights: weight ∝ totalDiff ^ selectionBias
#pragma omp parallel for num_threads(this->threads)
    for (int i = 0; i < c; i++) {
        int totalDiff = accumulate(this->diffMatrix[i].begin(), this->diffMatrix[i].end(), 0);
        double weight = pow(static_cast<double>(totalDiff), this->selectionBias);
        this->selectionWeights[i] = weight;
        this->diffVec[i] = make_tuple(weight / this->c, i);
    }
    sort(this->diffVec.begin(), this->diffVec.end());
    this->totalWeight = accumulate(this->selectionWeights.begin(), this->selectionWeights.end(), 0.0);
}

// ─────────────────────────────────────────────
//  Construction
// ─────────────────────────────────────────────

// Generates lkhS solutions by running LKH and stores them in lkhInitialSolutions.
void CBMProblem::createLKHInitialS() {
#pragma omp parallel for
    for (int i = 0; i < this->lkhS; i++) {
        string threadSuffix = to_string(omp_get_thread_num());
        if (!this->lkhCache) {
            this->toTSP(threadSuffix);
            this->initialTour(threadSuffix);
            this->runLKH(threadSuffix);
        }
        {
            lock_guard<mutex> lock(this->ISMutex);
            this->lkhInitialSolutions.push_back(this->fromTSP(threadSuffix));
        }
    }
}

CBMSol CBMProblem::construction() {
    CBMSol s;
    s.cost = 0;

    if (this->ISCount) {
        // Consume one pre-computed LKH solution
        this->ISCount--;
        lock_guard<mutex> lock(this->ISMutex);
        s.sol = this->lkhInitialSolutions.back();
        s.construction = LKH;
        this->lkhInitialSolutions.pop_back();
        this->evaluate(s);
        this->initialSolutions.push_back(s);
    } else {
        switch (constructionMethod) {
            case 1:
                s = this->greedyConstruction();
                this->evaluate(s);
                s.construction = GREEDY;
                break;
            case 2:
                // Random permutation improved by one-block local search
                s.sol.resize(this->c);
                iota(s.sol.begin(), s.sol.end(), 0);
                shuffle(s.sol.begin(), s.sol.end(), this->mersenne_engine);
                this->evaluate(s);
                s = this->oneBlockSearch(s);
                s.construction = ONEBLOCK;
                break;
        }
        lock_guard<mutex> lock(this->ISMutex);
        this->initialSolutions.push_back(s);
    }

    return s;
}

// Greedy nearest-neighbour construction with biased randomisation.
CBMSol CBMProblem::greedyConstruction() {
    uniform_int_distribution<> colDist(0, this->c - 1);
    unordered_set<int> remaining;
    CBMSol s;

    s.cost = 0;
    s.sol.resize(this->c);
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

// Picks the next column to append using biased-random selection over sorted candidates.
int CBMProblem::nextInsertion(int& current, unordered_set<int>& remaining) {
    // Build similarity list: more similar (fewer diffs) = ranked higher
    vector<tuple<int, int>> candidates;
    candidates.reserve(remaining.size());
    for (int i : remaining) candidates.push_back({this->l - this->diffMatrix[current][i], i});

    sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return get<0>(a) > get<0>(b); });

    // Weight by rank: rank-1 gets weight 1/1^bias, rank-2 gets 1/2^bias, etc.
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

// ─────────────────────────────────────────────
//  One-block local search
// ─────────────────────────────────────────────

// Iterative improvement: for each row, tries moving elements of one-blocks
// closer to other one-blocks until no improving move is found.
CBMSol CBMProblem::oneBlockSearch(CBMSol& s, vector<int> rows) {
    bool searchAllRows = rows.empty();
    if (searchAllRows) {
        rows.resize(this->c);
        iota(rows.begin(), rows.end(), 0);
        shuffle(rows.begin(), rows.end(), this->mersenne_engine);
    }

    int best = s.cost;
    bool improved = true;

    while (improved) {
        improved = false;
        for (size_t ri = 0; ri < rows.size() && !improved; ri++) {
            vector<pair<int, int>> blocks = this->findOneBlocks(s, rows[ri]);
            int numBlocks = blocks.size();

            for (int i = 0; i < numBlocks && !improved; i++) {
                for (int j = 0; j < numBlocks && !improved; j++) {
                    if (i == j) continue;
                    bool isLastPair = searchAllRows && (i == numBlocks - 1) && (j == numBlocks - 1);
                    improved = this->moveOneBlockColumns(s, best, blocks[i], blocks[j], isLastPair);
                }
            }
        }
    }

    return s;
}

// Randomly moves one element from one block to inside another (used as a neighbourhood move).
CBMSol CBMProblem::oneBlockMovement(CBMSol& s) {
    uniform_int_distribution<> rowDist(0, this->l - 1);
    const int maxAttempts = 10;

    int row = rowDist(this->mersenne_engine);
    vector<pair<int, int>> blocks = this->findOneBlocks(s, row);

    // Retry until we find a row with at least two blocks to move between
    for (int attempt = 0; attempt < maxAttempts && blocks.size() < 2; attempt++) {
        row = rowDist(this->mersenne_engine);
        blocks = this->findOneBlocks(s, row);
    }
    if (blocks.size() < 2) return s;

    // Randomly pick two distinct blocks
    vector<int> indices(blocks.size());
    iota(indices.begin(), indices.end(), 0);
    shuffle(indices.begin(), indices.end(), this->mersenne_engine);

    moveOneBlockRandomly(s, blocks[indices[0]], blocks[indices[1]]);
    return s;
}

// Picks a random element from b1 and reinserts it at a random position inside b2.
void CBMProblem::moveOneBlockRandomly(CBMSol& s, pair<int, int> b1, pair<int, int> b2) {
    auto leftOf = [&](int i) { return (i > 0) ? s.sol[i - 1] : -1; };
    auto rightOf = [&](int i) { return (i + 1 < this->c) ? s.sol[i + 1] : -1; };

    uniform_int_distribution<int> srcDist(b1.first, b1.second);
    uniform_int_distribution<int> dstDist(b2.first, b2.second + 1);
    int from = srcDist(mersenne_engine);
    int to = dstDist(mersenne_engine);

    int elem = s.sol[from];
    int origLeft = leftOf(from);
    int origRight = rightOf(from);
    int dest = moveHelper(s.sol, from, to);
    int destLeft = leftOf(dest);
    int destRight = rightOf(dest);

    s.mE = {origLeft, origRight, destLeft, destRight, elem};
    s.movement = REINSERTION;
    this->deltaEval(s);
}

// Tries every (source, destination) pair between two blocks; accepts the first improvement.
bool CBMProblem::moveOneBlockColumns(CBMSol& s, int& best, pair<int, int> srcBlock, pair<int, int> dstBlock, bool acceptAnyway) {
    auto leftOf = [&](int i) { return (i > 0) ? s.sol[i - 1] : -1; };
    auto rightOf = [&](int i) { return (i + 1 < this->c) ? s.sol[i + 1] : -1; };

    int savedCost = s.cost;
    vector<int> savedSol = s.sol;

    for (int i = srcBlock.first; i <= srcBlock.second; i++) {
        for (int j = dstBlock.first; j <= dstBlock.second + 1; j++) {
            int elem = s.sol[i];
            int origLeft = leftOf(i);
            int origRight = rightOf(i);
            int dest = moveHelper(s.sol, i, j);
            int destLeft = leftOf(dest);
            int destRight = rightOf(dest);

            s.mE = {origLeft, origRight, destLeft, destRight, elem};
            s.movement = REINSERTION;
            this->deltaEval(s);

            if (s.cost < best) {
                best = s.cost;
                return true;  // first-improvement: accept immediately
            }

            // No improvement — roll back
            bool isVeryLast = acceptAnyway && (i == srcBlock.second) && (j == dstBlock.second);
            if (isVeryLast) return false;
            s.sol = savedSol;
            s.cost = savedCost;
        }
    }

    return false;
}

// Returns all contiguous runs of 1s in row `row` under the current column permutation.
vector<pair<int, int>> CBMProblem::findOneBlocks(CBMSol& s, int row) {
    vector<pair<int, int>> blocks;

    for (int i = 0; i < this->c; i++) {
        if (this->binaryMatrix[row][s.sol[i]]) {
            int start = i;
            while (i < this->c && this->binaryMatrix[row][s.sol[i]]) i++;
            blocks.push_back({start, i - 1});
        }
    }

    return blocks;
}

// ─────────────────────────────────────────────
//  LKH interface
// ─────────────────────────────────────────────

// Writes the CBM instance as a TSP FULL_MATRIX problem that LKH can read.
void CBMProblem::toTSP(string suffix) {
    filesystem::path tspPath = this->currPath / "instances" / "tsp" / (this->instanceName + suffix + ".tsp");
    filesystem::path parPath = this->currPath / "instances" / "tsp" / (this->instanceName + suffix + ".par");

    ofstream tsp(tspPath), par(parPath);
    if (!tsp || !par) throw runtime_error("Error creating TSP/PAR files");

    tsp << "NAME : " << this->instanceName << "\n"
        << "TYPE : TSP\n"
        << "DIMENSION : " << this->c + 1 << "\n"
        << "EDGE_WEIGHT_TYPE : EXPLICIT\n"
        << "EDGE_WEIGHT_FORMAT : FULL_MATRIX\n"
        << "EDGE_WEIGHT_SECTION\n";

    for (int i = 0; i <= this->c; ++i) {
        for (int j = 0; j <= this->c; ++j) tsp << this->tspMatrix[i][j] << (j == this->c - 1 ? "\n" : " ");
    }
    tsp << "EOF\n";

    string solFile = (this->currPath / "instances" / "tsp" / (this->instanceName + suffix + ".sol")).string();
    string initialTour = (this->currPath / "instances" / "tsp" / (this->instanceName + suffix + "_initial_.tour")).string();

    par << "PROBLEM_FILE = " << tspPath.string() << "\n"
        << "OUTPUT_TOUR_FILE = " << solFile << "\n"
        << "RUNS = 1\n"
        << "INITIAL_TOUR_FILE = " << initialTour << "\n"
        << "TIME_LIMIT = " << this->lkhMaxTime << "\n";
}

// Writes a greedy initial tour file that LKH uses as its starting point.
void CBMProblem::initialTour(string suffix) {
    filesystem::path tourPath = this->currPath / "instances" / "tsp" / (this->instanceName + suffix + "_initial_.tour");
    ofstream tour(tourPath);
    if (!tour) throw runtime_error("Error creating TOUR file");

    CBMSol s = this->greedyConstruction();
    tour << "TOUR_SECTION\n"
         << "1\n";
    for (int i = 0; i < this->c; i++) tour << s.sol[i] + 2 << "\n";  // +2: shift from 0-indexed to TSP node numbering (depot = 1)
    tour << "-1\nEOF";
}

// Invokes the LKH binary on the generated parameter file.
void CBMProblem::runLKH(string suffix) {
    filesystem::path parPath = this->currPath / "instances" / "tsp" / (this->instanceName + suffix + ".par");

    if (!fs::exists(this->lkhPath)) throw runtime_error("LKH executable not found: " + lkhPath.string());
    if (!fs::exists(parPath)) throw runtime_error("PAR file not found: " + parPath.string());

    string command = lkhPath.string() + " " + parPath.string() + " > /dev/null 2>&1";
    int ret = system(command.c_str());
    if (ret != 0) cerr << "LKH returned error code: " << ret << endl;
}

// Reads the LKH output tour and converts it back to a 0-indexed column permutation.
vector<int> CBMProblem::fromTSP(string suffix) {
    filesystem::path solPath = this->lkhCache ? this->currPath / "instances" / "tsp" / "cache" / (this->instanceName + suffix + ".sol")
                                              : this->currPath / "instances" / "tsp" / (this->instanceName + suffix + ".sol");

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
            tour.push_back(node - 2);  // -2: undo TSP node numbering, depot (node 1) becomes -1
        }
    }

    // If the depot (-1) ended up inside the tour, rotate it to the front and remove it
    auto depotIt = find(tour.begin(), tour.end(), -1);
    if (depotIt != tour.end()) {
        rotate(tour.begin(), depotIt, tour.end());
        tour.erase(tour.begin());
    }

    return tour;
}

// ─────────────────────────────────────────────
//  Debug
// ─────────────────────────────────────────────

void CBMProblem::printS(CBMSol& s) {
    for (int row = 0; row < this->l; row++) {
        for (int col = 0; col < this->c; col++) cout << this->binaryMatrix[row][s.sol[col]] << " ";
        cout << "\n";
    }
}