#include "CBMProblem.hpp"

CBMProblem::CBMProblem(string filename, int movementType, double constructionBias): mersenne_engine(rng_device()), constructionBias(constructionBias) {
    ifstream input(filename);
    if (!input.is_open()) {
        throw runtime_error("Error opening file: " + filename);
    }

    input >> this->l;
    input >> this->c;

    this->binaryMatrix.resize(this->l);
    this->diffMatrix.resize(this->c, vector<int>(this->c, 0));
    this->onesToZeros.resize(this->c, vector<int>(this->c, 0));
    this->zerosToOnes.resize(this->c, vector<int>(this->c, 0));
    this->onesToOnes.resize(this->c, vector<int>(this->c, 0));

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

CBMSol CBMProblem::neighbor(CBMSol s) {
    uniform_int_distribution<> dist(0, this->c - 1);
    uniform_int_distribution<> movementDist(1, 3);
    int index = dist(this->mersenne_engine);
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
        case 3:
            s.movement = REINSERTION;
            int element = s.sol[index];
            int oL = getLeft(index);
            int oR = getRight(index);
            s.sol.erase(s.sol.begin() + index);
            s.sol.insert(s.sol.begin() + newIndex, element);
            int nL = getLeft(newIndex);
            int nR = getRight(newIndex);
            s.mE = {oL, oR, nL, nR, element};
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

    auto twoOptLambda = [&](int i, int j) -> int {
        if(i == -1 || j == -1)
            return 0;
        return this->onesToOnes[i][j];
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
            int prevLeft = twoOptLambda(s.mE[1], s.mE[0]);
            int prevRight = twoOptLambda(s.mE[2], s.mE[3]);
            int afterLeft = twoOptLambda(s.mE[2], s.mE[0]);
            int afterRight = twoOptLambda(s.mE[1], s.mE[3]);
            s.cost += (prevLeft + prevRight) - (afterLeft + afterRight);
            break;
        }
    }
    return s.cost;
}

void CBMProblem::computeMatrixes() {
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
            this->diffMatrix[i][j]  = o2z + z2o;
        }
    }
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

CBMSol CBMProblem::construction() {
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