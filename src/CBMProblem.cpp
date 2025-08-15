#include "CBMProblem.hpp"

CBMProblem::CBMProblem(string filename, int movementType, double constructionBias): mersenne_engine(rng_device()), constructionBias(constructionBias) {
    ifstream input(filename);
    if (!input.is_open()) {
        throw runtime_error("Error opening file: " + filename);
    }

    input >> this->l;
    input >> this->c;

    this->binaryMatrix.resize(this->l);
    this->B.resize(this->l);
    this->W.resize(this->c + 2, vector<int>(this->c + 2, false));

    int n, e;
    for (int i = 0; i < this->l; i++) {
        input >> n;
        for (int j = 0; j < n; j++) {
            input >> e;
            this->binaryMatrix[i][e - 1] = true;
        }
    }
    for (int i = 0; i < this->l; i++) {
        for (int j = 0; j < this->c; j++) {
            this->B[i][j + 1] = this->binaryMatrix[i][j];
        }
    }

    this->computeW();
    this->movementType = movementType;
}

CBMSol CBMProblem::neighbor(CBMSol s) {
    uniform_int_distribution<> dist(0, this->c - 1);
    int index = dist(this->mersenne_engine);
    int newIndex = dist(this->mersenne_engine);

    auto getLeft = [&](int idx) { return (idx > 0) ? s.sol[idx - 1] : -1; };
    auto getRight = [&](int idx) { return (idx + 1 < this->c) ? s.sol[idx + 1] : -1; };

    switch (this->movementType) {
        case 1:
            s.movement = SWAP;
            swap(s.sol[index], s.sol[newIndex]);
            s.mE = {getLeft(index), getRight(index), getLeft(newIndex), getRight(newIndex), s.sol[index]};
            break;
        case 2:
            s.movement = TWOOPT;
            if (index > newIndex) swap(index, newIndex);
            reverse(s.sol.begin() + index, s.sol.begin() + newIndex + 1);
            s.mE = {getLeft(index), getRight(index), getLeft(newIndex), getRight(newIndex), s.sol[index]};
            break;
        case 3:
            s.movement = REINSERTION;
            if (index < newIndex) newIndex -= 1;
            int element = s.sol[index];
            s.sol.erase(s.sol.begin() + index);
            s.sol.insert(s.sol.begin() + newIndex, element);
            s.mE = {getLeft(index), getRight(index), getLeft(newIndex), getRight(newIndex), s.sol[index]};
    }

    return s;
}

int CBMProblem::deltaEvaluate(CBMSol s) {
    switch(s.movement) {
        case REINSERTION: {
            int prev = this->W[get<0>(s.mE)][get<4>(s.mE)] +
                        this->W[get<4>(s.mE)][get<1>(s.mE)] - 
                        this->W[get<0>(s.mE)][get<1>(s.mE)];
            int after = this->W[get<2>(s.mE)][get<4>(s.mE)] +
                        this->W[get<4>(s.mE)][get<3>(s.mE)] - 
                        this->W[get<2>(s.mE)][get<3>(s.mE)];
            s.cost += -(prev - after) / 2;
            break;
        }
        case SWAP: {
            break;
        }
        case TWOOPT: {
            break;
        }
    }
    return s.cost;
}

int CBMProblem::evaluate(CBMSol s) {
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
            similarity.push_back({this->W[curr][i], i});

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

void CBMProblem::computeW() {
    int similarity;
    for(int i = 0 ; i <= this->c ; i++) {
        for(int j = i + 1 ; j <= this->c ; j++) {
            similarity = this->calculateHamming(i, j);
            this->W[i][j] = similarity;
            this->W[j][i] = similarity;
        }
    }
}

int CBMProblem::calculateHamming(int& curr, int& candidate) {
    int diffs = 0;
    for (int i = 0; i < l; ++i)
        diffs += (B[i][curr] ^ B[i][candidate]);
    return diffs;
}

void CBMProblem::printMatrix(const CBMSol* s) {
    if (!s) {
        for (int row = 0; row < this->l; row++) {
            for (int col = 0; col < this->c; col++) {
                cout << this->binaryMatrix[row][col] << " ";
            }
            cout << "\n";
        }
    } else {
        for (int row = 0; row < this->l; row++) {
            for (int col = 0; col < this->c; col++) {
                cout << this->binaryMatrix[row][s->sol[col]] << " ";
            }
            cout << "\n";
        }
    }
}