#include "LKHWrapper.hpp"

LKHWrapper::LKHWrapper(int c, int l, const vector<vector<int>>& tspMatrix) : c(c), l(l), tspMatrix(tspMatrix) {
    if (!fs::exists(tmpDir)) {
        fs::create_directories(tmpDir);
    }
}

vector<int> LKHWrapper::run(const vector<int>& slice, string instanceName, int maxTime) {
    long execId = getExecutionId();
    string tspFile = LKHWrapper::tmpDir + instanceName + "_tsp_" + to_string(execId) + ".tsp";
    string parFile = LKHWrapper::tmpDir + instanceName + "_par_" + to_string(execId) + ".par";
    string resultTourFile = LKHWrapper::tmpDir + instanceName + "_result_tour_" + to_string(execId) + ".tsp";
    string InitialTourFile = LKHWrapper::tmpDir + instanceName + "_initial_tour_" + to_string(execId) + ".tsp";

    writeTSP(slice, instanceName, tspFile, execId);
    writeInitialTour(slice, instanceName, InitialTourFile, execId);
    writePar(parFile, tspFile, InitialTourFile, resultTourFile, maxTime);
    runLKH(parFile);

    vector<int> resultTour = getResultTour(resultTourFile);
    vector<int> mappedTour(resultTour.size());
    for (size_t i = 0; i < resultTour.size(); i++) {
        mappedTour[i] = slice[resultTour[i]];
    }

    return mappedTour;
}

void LKHWrapper::writePar(string parFile, string tspFile, string InitialTourFile, string resultTourFile, int maxTime) {
    ofstream out(parFile);
    if (!out.is_open()) throw runtime_error("Error opening file: " + parFile);

    out << "PROBLEM_FILE = " << tspFile << endl;
    out << "INITIAL_TOUR_FILE = " << InitialTourFile << endl;
    out << "TOUR_FILE = " << resultTourFile << endl;
    out << "MOVE_TYPE = 5" << endl;
    out << "PATCHING_C = 3" << endl;
    out << "PATCHING_A = 2" << endl;
    out << "RUNS = 1" << endl;
    out << "TIME_LIMIT = " << maxTime << endl;

    out.close();
}

void LKHWrapper::writeInitialTour(const vector<int>& slice, string instanceName, string tourFile, long execId) {
    ofstream out(tourFile);
    if (!out.is_open()) throw runtime_error("Error opening file: " + tourFile);

    out << "NAME : " << instanceName << "_tour_" << execId << endl;
    out << "TYPE : TOUR" << endl;
    out << "DIMENSION : " << slice.size() + 1 << endl;
    out << "TOUR_SECTION" << endl;

    out << 1 << endl;
    for (size_t i = 0; i < slice.size(); i++) {
        out << slice[i] + 2 << endl;
    }
    out << -1 << endl;
    out << "EOF" << endl;

    out.close();
}

void LKHWrapper::writeTSP(const vector<int>& slice, string instanceName, string tspFile, long execId) {
    ofstream out(tspFile);
    if (!out.is_open()) throw runtime_error("Error opening file: " + tspFile);

    out << "NAME : " << instanceName << "_tsp_" << execId << endl;
    out << "TYPE : TSP" << endl;
    out << "DIMENSION : " << slice.size() + 1 << endl;
    out << "EDGE_WEIGHT_TYPE : EXPLICIT" << endl;
    out << "EDGE_WEIGHT_FORMAT : FULL_MATRIX" << endl;
    out << "EDGE_WEIGHT_SECTION" << endl;

    out << tspMatrix[0][0];
    for (size_t j = 0; j < slice.size(); j++) out << " " << tspMatrix[0][slice[j] + 1];
    out << endl;

    for (size_t i = 0; i < slice.size(); i++) {
        out << tspMatrix[slice[i] + 1][0];
        for (size_t j = 0; j < slice.size(); j++) {
            out << " " << tspMatrix[slice[i] + 1][slice[j] + 1];
        }
        out << endl;
    }

    out.close();
}

vector<int> LKHWrapper::getResultTour(string resultTourFile) {
    ifstream sol(resultTourFile);
    if (!sol) throw runtime_error("Error opening solution file: " + resultTourFile);

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
            tour.push_back(node - 2);
        }
    }

    auto depotIt = find(tour.begin(), tour.end(), -1);
    if (depotIt != tour.end()) {
        rotate(tour.begin(), depotIt, tour.end());
        tour.erase(tour.begin());
    }

    return tour;
}

void LKHWrapper::runLKH(string parFile) {
    if (!fs::exists(LKHWrapper::lkhPath)) throw runtime_error("LKH executable not found: " + lkhPath);
    if (!fs::exists(parFile)) throw runtime_error("PAR file not found: " + parFile);
    string command = LKHWrapper::lkhPath + " " + parFile + " > /dev/null 2>&1";
    int ret = system(command.c_str());
    if (ret != 0) cerr << "LKH returned error code: " << ret << endl;
}

long LKHWrapper::getExecutionId() {
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<long> dist(1, LONG_MAX);
    return dist(gen);
}