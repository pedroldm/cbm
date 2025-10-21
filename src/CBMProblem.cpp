#include "CBMProblem.hpp"

CBMProblem::CBMProblem(string filename,
                       int movementType,
                       int constructionMethod,
                       double constructionBias,
                       double selectionBias,
                       int maxBlockSize,
                       int threads,
                       int lkhS,
                       int lkhMaxTime)
    : mersenne_engine(rng_device()),
      constructionMethod(constructionMethod),
      constructionBias(constructionBias),
      selectionBias(selectionBias),
      maxBlockSize(maxBlockSize),
      threads(threads),
      ISCount(lkhS),
      lkhMaxTime(lkhMaxTime){    
    ifstream input(filename);
    if (!input.is_open()) {
        throw runtime_error("Error opening file: " + filename);
    }
    auto pos = filename.find_last_of("/\\");
    this->lkhS = lkhS;
    this->instanceName = (pos == string::npos) 
        ? filename 
        : filename.substr(pos + 1);
    this->currPath = filesystem::current_path();
    this->lkhPath = this->currPath / "src" / "LKH3" / "LKH";

    input >> this->l;
    input >> this->c;

    this->binaryMatrix.resize(this->l);
    this->diffMatrix.resize(this->c, vector<int>(this->c, 0));
    this->tspMatrix.resize(this->c + 1, vector<int>(this->c + 1, 0));
    this->onesToZeros.resize(this->c, vector<int>(this->c, 0));
    this->zerosToOnes.resize(this->c, vector<int>(this->c, 0));
    this->onesToOnes.resize(this->c, vector<int>(this->c, 0));
    this->diffVec.resize(this->c);
    this->selectionWeights.resize(this->c);

    int n, e;
    for (int i = 0; i < this->l; i++) {
        input >> n;
        for (int j = 0; j < n; j++) {
            input >> e;
            this->binaryMatrix[i][e - 1] = true;
        }
    }

    this->computeMatrixes();
    this->movementType = movementType;
}

int CBMProblem::biasedSelection() {
    uniform_real_distribution<double> dist(0.0, this->totalWeight);
    double rnd = dist(this->mersenne_engine);

    double cum_sum = 0.0;
    size_t chosen_idx = 0;
    for (; chosen_idx < this->c; ++chosen_idx) {
        cum_sum += this->selectionWeights[chosen_idx];
        if (rnd < cum_sum) break;
    }

    auto [finalDiff, chosen] = this->diffVec[chosen_idx];
    return chosen;
}

CBMSol CBMProblem::neighbor(CBMSol s) {
    uniform_int_distribution<> dist(0, this->c - 1);
    uniform_int_distribution<> distL(0, this->l - 1);
    uniform_int_distribution<> movementDist(1, 3);
    
    int index = (this->selectionBias > 1) ? this->biasedSelection() : dist(this->mersenne_engine);
    int newIndex = dist(this->mersenne_engine);

    auto getLeft = [&](int idx) { return (idx > 0) ? s.sol[idx - 1] : -1; };
    auto getRight = [&](int idx) { return (idx + 1 < this->c) ? s.sol[idx + 1] : -1; };

    int movement = (this->movementType < 4) ? this->movementType : movementDist(this->mersenne_engine);

    switch (movement) {
        case 1:
            if (index > newIndex) swap(index, newIndex);
            if(abs(newIndex - index) > 1) {
                s.movement = SWAP;
                s.mE = {getLeft(index), getRight(index), getLeft(newIndex), getRight(newIndex), s.sol[index], s.sol[newIndex]};
                swap(s.sol[index], s.sol[newIndex]);
            }
            else {
                s.movement = TWOOPT;
                s.mE = {getLeft(index), s.sol[index], s.sol[newIndex], getRight(newIndex)};
                swap(s.sol[index], s.sol[newIndex]);
            }
            break;
        case 2:
            if (index > newIndex) swap(index, newIndex);
            s.movement = TWOOPT;
            s.mE = {getLeft(index), s.sol[index], s.sol[newIndex], getRight(newIndex)};
            reverse(s.sol.begin() + index, s.sol.begin() + newIndex + 1);
            break;
        case 3: {
            s.movement = REINSERTION;
            int element = s.sol[index];
            int oL = getLeft(index);
            int oR = getRight(index);
            s.sol.erase(s.sol.begin() + index);
            s.sol.insert(s.sol.begin() + newIndex, element);
            int nL = getLeft(newIndex);
            int nR = getRight(newIndex);
            s.mE = {oL, oR, nL, nR, element};
            break;
        }
        case 4: {
            int line = distL(this->mersenne_engine);
            this->oneBlockSearch(s, {line});
            s.movement = ONEBLOCKM;
            break;
        }
    }

    return s;
}

int CBMProblem::deltaEval(CBMSol& s) {
    auto lambda = [&](int L, int R, int e) -> int {
        if (L != -1 && R != -1) {
            return (this->diffMatrix[L][e] +
                    this->diffMatrix[e][R] -
                    this->diffMatrix[L][R]) / 2;
        } else if (L == -1) {
            return this->onesToZeros[e][R];
        } else if (R == -1) {
            return this->onesToZeros[e][L];
        }
    };

    auto andLambda = [&](int i, int j) -> int {
        if(i == -1 || j == -1)
            return 0;
        return this->onesToOnes[i][j];
    };

    auto blockLambda = [&](int i, int j) -> int {
        return (i != -1 && j != -1) ? this->diffMatrix[i][j] : 0;
    };

    switch(s.movement) {
        case REINSERTION: {
            int prev = lambda(s.mE[0], s.mE[1], s.mE[4]);
            int after = lambda(s.mE[2], s.mE[3], s.mE[4]);
            s.cost += (-prev + after);
            break;
        }
        case SWAP: {
            int prevLeft = lambda(s.mE[0], s.mE[1], s.mE[4]);
            int prevRight = lambda(s.mE[2], s.mE[3], s.mE[5]);
            int afterLeft = lambda(s.mE[2], s.mE[3], s.mE[4]);
            int afterRight = lambda(s.mE[0], s.mE[1], s.mE[5]);
            s.cost += -(prevLeft + prevRight) + (afterLeft + afterRight);
            break;
        }
        case TWOOPT: {
            int prevLeft = andLambda(s.mE[1], s.mE[0]);
            int prevRight = andLambda(s.mE[2], s.mE[3]);
            int afterLeft = andLambda(s.mE[2], s.mE[0]);
            int afterRight = andLambda(s.mE[1], s.mE[3]);
            s.cost += (prevLeft + prevRight) - (afterLeft + afterRight);
            break;
        }
    }
    return s.cost;
}

void CBMProblem::computeMatrixes() {
    #pragma omp parallel for num_threads(this->threads) collapse(2)
    for (int i = 0; i < c; i++) {
        for (int j = 0; j < c; j++) {
            if (i == j) continue;
            int o2z = 0, z2o = 0, o2o = 0;
            for (int row = 0; row < l; row++) {
                if (binaryMatrix[row][i] && !binaryMatrix[row][j]) o2z++;
                else if (!binaryMatrix[row][i] && binaryMatrix[row][j]) z2o++;
                else if (binaryMatrix[row][i] && binaryMatrix[row][j]) o2o++;
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
            if (binaryMatrix[row][i-1])
                onesCount++;
        this->tspMatrix[0][i] = onesCount;
        this->tspMatrix[i][0] = onesCount;
    }

    #pragma omp parallel for num_threads(this->threads)
    for(int i = 0 ; i < c ; i++) {
        int diffCount = accumulate(this->diffMatrix[i].begin(), this->diffMatrix[i].end(), 0);
        double weight = pow(static_cast<double>(diffCount), this->selectionBias);
        double biasedCount = weight / this->c;
        this->selectionWeights[i] = weight;
        this->diffVec[i] = make_tuple(biasedCount, i);
    }
    sort(this->diffVec.begin(), this->diffVec.end());
    this->totalWeight = accumulate(this->selectionWeights.begin(), this->selectionWeights.end(), 0.0);
}

int CBMProblem::evaluate(CBMSol& s) {
    return (s.cost) ? this->deltaEval(s) : this->completeEval(s);
}

int CBMProblem::completeEval(CBMSol& s) {
    s.cost = 0;

    for (int row = 0; row < this->l; row++) {
        if (this->binaryMatrix[row][s.sol[0]]) {
            s.cost++;
        }
    }
    for (int col = 1; col < this->c; col++) {
        for (int row = 0; row < this->l; row++) {
            if (this->binaryMatrix[row][s.sol[col]] &&
                !this->binaryMatrix[row][s.sol[col - 1]])
                s.cost++;
        }
    }

    return s.cost;
}

void CBMProblem::createLKHInitialS() {
    #pragma omp parallel for
    for (int i = 0; i < this->lkhS; i++) {
        string tid = to_string(omp_get_thread_num());
        this->toTSP(tid);
        this->initialTour(tid);
        this->runLKH(tid);
        vector<int> s = this->fromTSP(tid);
        {
            lock_guard<mutex> lock(this->ISMutex);
            this->lkhInitialSolutions.push_back(s);
        }
    }
}

CBMSol CBMProblem::construction() {
    CBMSol s;
    s.cost = 0;
    if(this->ISCount) {
        this->ISCount--;
        {
            lock_guard<mutex> lock(this->ISMutex);
            s.sol = this->lkhInitialSolutions.back();
            this->evaluate(s);
            s.construction = LKH;
            this->lkhInitialSolutions.pop_back();
            this->initialSolutions.push_back(s);
        }
    } else {
        switch(constructionMethod) {
            case 1:
                s = this->greedyConstruction();
                this->evaluate(s);
                s.construction = GREEDY;
                break;
            case 2:
                s.sol.resize(this->c);
                iota(s.sol.begin(), s.sol.end(), 0);
                shuffle(s.sol.begin(), s.sol.end(), this->mersenne_engine);
                this->evaluate(s);
                s = this->oneBlockSearch(s);
                s.construction = ONEBLOCK;
                break;
        }
        {
            lock_guard<mutex> lock(this->ISMutex);
            this->initialSolutions.push_back(s);
        }
    }
    return s;
}

CBMSol CBMProblem::greedyConstruction() {
    uniform_int_distribution<> dist(0, this->c - 1);
    unordered_set<int> out;
    CBMSol ss;

    ss.cost = 0;
    ss.sol.resize(this->c);
    int curr = dist(this->mersenne_engine);

    for (int i = 0; i < this->c; i++) 
        out.insert(i);

    int idx = 0;
    out.erase(curr);
    ss.sol[idx++] = curr;
    while (!out.empty()) {
        int nI = this->nextInsertion(curr, out);
        ss.sol[idx++] = nI;
        out.erase(nI);
        curr = nI;
    }

    return ss;
}

int CBMProblem::nextInsertion(int& curr, unordered_set<int>& out) {
    vector<tuple<int, int>> similarity;
    for (int i = 0; i < this->c; i++)
        if (out.count(i))
            similarity.push_back({this->l - this->diffMatrix[curr][i], i});

    sort(similarity.begin(), similarity.end(),
         [](const tuple<int, int>& a, const tuple<int, int>& b) {
             return get<0>(a) > get<0>(b);
         });

    vector<double> weights(similarity.size());
    for (size_t i = 0; i < similarity.size(); ++i) {
        weights[i] = 1.0 / pow(i + 1, this->constructionBias);
    }

    double total_weight = accumulate(weights.begin(), weights.end(), 0.0);
    uniform_real_distribution<double> dist(0.0, total_weight);
    double rnd = dist(this->mersenne_engine);
    double cum_sum = 0.0;
    size_t chosen_idx = 0;

    for (; chosen_idx < weights.size(); ++chosen_idx) {
        cum_sum += weights[chosen_idx];
        if (rnd < cum_sum) break;
    }

    return get<1>(similarity[chosen_idx]);
}

void CBMProblem::printS(CBMSol& s) {
    for (int row = 0; row < this->l; row++) {
        for (int col = 0; col < this->c; col++) {
            cout << this->binaryMatrix[row][s.sol[col]] << " ";
        }
        cout << endl;
    }
}

CBMSol CBMProblem::oneBlockSearch(CBMSol& s, vector<int> lines) {
    bool returnAnyway = false;
    if (lines.empty()) {
        lines.resize(this->c);
        iota(lines.begin(), lines.end(), 0);
        shuffle(lines.begin(), lines.end(), this->mersenne_engine);
        returnAnyway = true;
    }

    int currBest = s.cost;
    bool improved = true;
    while(improved) {
        improved = false;
        for (size_t li = 0; li < lines.size() && !improved; li++) {
            vector<pair<int, int>> oneBlocks = this->findOneBlocks(s, lines[li]);
            int size = oneBlocks.size();
            for(int i = 0 ; i < size && !improved ; i++) {
                for(int j = 0 ; j < size && !improved ; j++) {
                    if(i == j)
                        continue;
                    improved = this->moveOneBlockColumns(s, currBest, oneBlocks[i], oneBlocks[j], (returnAnyway && i == size - 1 && j == size - 1));
                }
            }
        }
    }

    return s;
}

bool CBMProblem::moveOneBlockColumns(CBMSol& s, int& currBest, pair<int, int> b1, pair<int, int> b2, bool returnAnyway) {
    auto getLeft = [&](int idx) { return (idx > 0) ? s.sol[idx - 1] : -1; };
    auto getRight = [&](int idx) { return (idx + 1 < this->c) ? s.sol[idx + 1] : -1; };
    int currCost = s.cost;
    vector<int> originalS = s.sol;

    for(int i = b1.first ; i <= b1.second ; i++) {
        for (int j = b2.first; j <= b2.second + 1; ++j) {
            int element = s.sol[i];
            int oL = getLeft(i);
            int oR = getRight(i);
            int to = moveHelper(s.sol, i, j);
            int nL = getLeft(to);
            int nR = getRight(to);
            s.mE = {oL, oR, nL, nR, element};
            s.movement = REINSERTION;
            this->deltaEval(s);
            if(currBest > s.cost){
                currBest = s.cost;
                return true;
            } else {
                if(returnAnyway && i == b1.second && j == b2.second)
                    return false;
                s.sol = originalS;
                s.cost = currCost;
            }
        }
    }

    return false;
}

vector<pair<int, int>> CBMProblem::findOneBlocks(CBMSol& s, int line) {
    vector<pair<int, int>> pairs;
    int start, end;

    for(int i = 0 ; i < this->c ; i++) {
        if(this->binaryMatrix[line][s.sol[i]]) {
            start = i;
            while(i < this->c && this->binaryMatrix[line][s.sol[i]])
                i++;
            end = i - 1;
            pairs.push_back({start, end});
        }
    }

    return pairs;
}

void CBMProblem::toTSP(string sufix) {
    filesystem::path tspPath = this->currPath / "instances" / "tsp" / (this->instanceName + sufix + ".tsp");
    filesystem::path parPath = this->currPath / "instances" / "tsp" / (this->instanceName + sufix + ".par");
    ofstream tsp(tspPath);
    ofstream par(parPath);
    if (!tsp || !par) {
        throw runtime_error("Error creating TSP/PAR files");
    }

    tsp << "NAME : " << this->instanceName << "\n";
    tsp << "TYPE : TSP\n";
    tsp << "DIMENSION : " << this->c + 1 << "\n";
    tsp << "EDGE_WEIGHT_TYPE : EXPLICIT\n";
    tsp << "EDGE_WEIGHT_FORMAT : FULL_MATRIX\n";
    tsp << "EDGE_WEIGHT_SECTION\n";

    for (int i = 0; i <= this->c; ++i) {
        for (int j = 0; j <= this->c; ++j) {
            tsp << this->tspMatrix[i][j] << (j == this->c - 1 ? "\n" : " ");
        }
    }
    tsp << "EOF\n";

    par << "PROBLEM_FILE = " << tspPath.string() << "\n";
    par << "OUTPUT_TOUR_FILE = " << (this->currPath / "instances" / "tsp" / (this->instanceName + sufix + ".sol")).string() << "\n";
    par << "RUNS = 1\n";
    par << "INITIAL_TOUR_FILE = " << (this->currPath / "instances" / "tsp" / (this->instanceName + sufix + "_initial_.tour")).string() << endl;
    par << "TIME_LIMIT = " << this->lkhMaxTime << endl;
}

void CBMProblem::runLKH(string sufix) {
    filesystem::path parPath = this->currPath / "instances" / "tsp" / (this->instanceName + sufix + ".par");
    if (!fs::exists(this->lkhPath)) {
        throw runtime_error("LKH executable not found: " + lkhPath.string());
    }
    if (!fs::exists(parPath)) {
        throw runtime_error("PAR file not found: " + parPath.string());
    }

    string command = lkhPath.string() + " " + parPath.string() + " > /dev/null 2>&1";

    int ret = system(command.c_str());
    if (ret != 0) {
        cerr << "LKH returned error code: " << ret << endl;
    }
}

void CBMProblem::initialTour(string sufix) {
    filesystem::path inputTourPath = this->currPath / "instances" / "tsp" / (this->instanceName + sufix + "_initial_.tour");
    ofstream tour(inputTourPath);
    if (!tour) {
        throw runtime_error("Error creating TOUR file");
    }
    tour << "TOUR_SECTION" << endl;
    CBMSol s = this->greedyConstruction();
    tour << "1" << endl;
    for(int i = 0 ; i < this->c ; i++)
        tour << s.sol[i] + 2 << endl;
    tour << "-1" << endl << "EOF";
}

vector<int> CBMProblem::fromTSP(string sufix) {
    filesystem::path solPath = this->currPath / "instances" / "tsp" / (this->instanceName + sufix + ".sol");
    ifstream sol(solPath);
    if (!sol) {
        throw runtime_error("Error opening solution file: " + solPath.string());
    }

    vector<int> tour;
    string line;
    bool inTourSection = false;

    while (getline(sol, line)) {
        if (!inTourSection) {
            if (line == "TOUR_SECTION") {
                inTourSection = true;
            }
        } else {
            istringstream iss(line);
            int node;
            while (iss >> node) {
                if (node == -1) {
                    inTourSection = false;
                    break;
                }
                tour.push_back(node - 2);
            }
        }
    }
    auto it = find(tour.begin(), tour.end(), -1);
    if (it != tour.end()) {
        rotate(tour.begin(), it, tour.end());
        tour.erase(tour.begin());
    }
    return tour;
}

