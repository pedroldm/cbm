// ─────────────────────────────────────────────────────────────────────────────
//  ILS for Consecutive Block Minimization
//
//  Usage:
//    ils --filePath=<instance> --algorithm=<linkern|lkh> --solverPath=<binary>
//        [--outputPath=<json>] [--seed=<n>] [--timeLimit=<s>] [--workDir=<dir>]
//        [--iterations=<NBITER>] [--phi=<phi>] [--V=<V>] [--runs=<LKH RUNS>]
//
//  Input format (1-indexed columns):
//    First line : l (rows)  c (cols)
//    Next l lines: n  e1 e2 ... en
//
//  Reproducibility: the perturbation stream and every TSP call derive their
//  seeds from --seed (call k of the run gets cbm_derive_seed(seed, k)), and
//  all scratch files live in a private --workDir. A failing TSP solver is an
//  error, not a silent fallback to a greedy tour.
// ─────────────────────────────────────────────────────────────────────────────

#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../src/common/cbm_lkh_params.h"
#include "../src/common/cbm_seed.h"
#include "../src/common/cbm_spawn.h"

using namespace std;
namespace fs = std::filesystem;

enum class Solver { LKH, LINKERN };

struct Options {
    string filePath, outputPath, solverPath, workDir;
    Solver solver = Solver::LINKERN;
    uint64_t seed = 1;
    double timeLimit = 7200.0;
    int iterations = 24;
    double phi = 0.2;
    int V = 1;
    int runs = 1;
};

static std::chrono::steady_clock::time_point g_start;
static double elapsedSeconds() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - g_start).count(); }

// ─────────────────────────────────────────────────────────────────────────────
//  TSP solver interface (one instance per run; calls are sequential)
// ─────────────────────────────────────────────────────────────────────────────
class TSPSolver {
   public:
    long calls = 0, timeLimitHits = 0;
    double seconds = 0.0;
    vector<uint32_t> seeds;

    TSPSolver(const Options& opt) : opt(opt) {}

    // Solves `dist` and returns a tour rotated to start at node 0. `callIndex`
    // is the logical identity of the call (0 = initial tour, k = iteration k)
    // and alone determines the solver seed.
    vector<int> solve(const vector<vector<int>>& dist, int callIndex) {
        const int n = static_cast<int>(dist.size());
        const string base = opt.workDir + "/call" + to_string(callIndex);
        const string tspFile = base + ".tsp", tourFile = base + ".tour", parFile = base + ".par", logFile = base + ".log";
        const uint32_t seed = cbm_derive_seed(opt.seed, static_cast<uint64_t>(callIndex));
        const double remaining = opt.timeLimit - elapsedSeconds();
        const int limit = remaining > 1.0 ? static_cast<int>(ceil(remaining)) : 1;
        seeds.push_back(seed);

        writeTSP(tspFile, dist);
        vector<string> args;
        if (opt.solver == Solver::LKH) {
            ofstream par(parFile);
            par << "PROBLEM_FILE = " << tspFile << "\n"
                << "OUTPUT_TOUR_FILE = " << tourFile << "\n"
                << "MOVE_TYPE = " << CBM_LKH_MOVE_TYPE << "\n"
                << "PATCHING_C = " << CBM_LKH_PATCHING_C << "\n"
                << "PATCHING_A = " << CBM_LKH_PATCHING_A << "\n"
                << "RUNS = " << opt.runs << "\n"
                << "TIME_LIMIT = " << limit << "\n"
                << "SEED = " << seed << "\n";
            if (!par) throw runtime_error("cannot write " + parFile);
            args = {opt.solverPath, parFile};
        } else {
            args = {opt.solverPath, "-Q", "-s", to_string(seed), "-t", to_string(limit), "-o", tourFile, tspFile};
        }

        vector<char*> argv;
        for (string& a : args) argv.push_back(a.data());
        argv.push_back(nullptr);
        const double start = elapsedSeconds();
        const int rc = cbm_run(argv.data(), logFile.c_str());
        const double took = elapsedSeconds() - start;
        calls++;
        seconds += took;
        if (rc != 0) throw runtime_error("TSP solver failed on call " + to_string(callIndex) + " (status " + to_string(rc) + "); see " + logFile);
        if (hitTimeLimit(logFile, took, limit)) timeLimitHits++;

        vector<int> tour = opt.solver == Solver::LKH ? parseTourLKH(tourFile) : parseTourLinkern(tourFile, n);
        if (static_cast<int>(tour.size()) != n)
            throw runtime_error("TSP solver returned a tour of " + to_string(tour.size()) + " nodes, expected " + to_string(n));
        auto it = find(tour.begin(), tour.end(), 0);
        if (it == tour.end()) throw runtime_error("TSP tour does not visit node 0");
        rotate(tour.begin(), it, tour.end());

        for (const string& f : {tspFile, tourFile, parFile, logFile}) fs::remove(f);
        return tour;
    }

   private:
    const Options& opt;

    static void writeTSP(const string& tspFile, const vector<vector<int>>& dist) {
        const int n = static_cast<int>(dist.size());
        ofstream f(tspFile);
        f << "NAME : tsp\n"
          << "TYPE : TSP\n"
          << "DIMENSION : " << n << "\n"
          << "EDGE_WEIGHT_TYPE : EXPLICIT\n"
          << "EDGE_WEIGHT_FORMAT : FULL_MATRIX\n"
          << "EDGE_WEIGHT_SECTION\n";
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) f << dist[i][j] << (j + 1 < n ? " " : "\n");
        f << "EOF\n";
        if (!f) throw runtime_error("cannot write " + tspFile);
    }

    // LKH reports a binding TIME_LIMIT explicitly; Linkern does not, so its
    // check falls back to the elapsed time.
    bool hitTimeLimit(const string& logFile, double took, int limit) const {
        if (opt.solver == Solver::LINKERN) return took >= 0.99 * limit;
        ifstream log(logFile);
        string line;
        while (getline(log, line))
            if (line.find("Time limit exceeded") != string::npos) return true;
        return false;
    }

    static vector<int> parseTourLKH(const string& tourFile) {
        ifstream f(tourFile);
        string line;
        bool inSection = false;
        vector<int> tour;
        while (getline(f, line)) {
            if (line.find("TOUR_SECTION") != string::npos) {
                inSection = true;
                continue;
            }
            if (!inSection) continue;
            istringstream iss(line);
            int v;
            while (iss >> v) {
                if (v == -1) return tour;
                tour.push_back(v - 1);  // LKH is 1-based
            }
        }
        return tour;
    }

    // Linkern edge-list format: "<n> <m>" then one "<a> <b> <w>" line per tour
    // edge, 0-based, starting at an arbitrary node.
    static vector<int> parseTourLinkern(const string& tourFile, int expectedNodes) {
        ifstream f(tourFile);
        int nNodes, nEdges;
        if (!(f >> nNodes >> nEdges) || nNodes != expectedNodes) return {};
        vector<int> tour;
        tour.reserve(nNodes);
        for (int i = 0; i < nNodes; i++) {
            int a, b, w;
            if (!(f >> a >> b >> w)) return {};
            tour.push_back(a);
        }
        return tour;
    }
};

// =============================================================================
struct HistoryEntry {
    int iteration;
    int blocks;
    double elapsed;
};

struct CBM {
    int l, c;
    vector<vector<bool>> A;
    vector<vector<int>> W;
    TSPSolver& tsp;
    long rejectedTours = 0;

    explicit CBM(TSPSolver& tsp) : tsp(tsp) {}

    void read(istream& in) {
        if (!(in >> l >> c) || l <= 0 || c <= 0) throw runtime_error("malformed instance header");
        A.assign(l, vector<bool>(c, false));
        for (int i = 0; i < l; i++) {
            int n;
            if (!(in >> n)) throw runtime_error("malformed instance row " + to_string(i + 1));
            for (int j = 0; j < n; j++) {
                int e;
                if (!(in >> e) || e < 1 || e > c) throw runtime_error("malformed instance row " + to_string(i + 1));
                A[i][e - 1] = true;
            }
        }
    }

    // Hamming distances between columns, with dummy all-zero columns 0 and c+1.
    void buildW() {
        int N = c + 2;
        W.assign(N, vector<int>(N, 0));
        auto getB = [&](int row, int col) -> int {
            if (col == 0 || col == c + 1) return 0;
            return A[row][col - 1] ? 1 : 0;
        };
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++) {
                int d = 0;
                for (int r = 0; r < l; r++) {
                    int bi = getB(r, i), bj = getB(r, j);
                    d += bi * (1 - bj) + (1 - bi) * bj;
                }
                W[i][j] = W[j][i] = d;
            }
    }

    int countBlocks(const vector<int>& pi) const {
        int blocks = 0;
        for (int r = 0; r < l; r++) {
            bool prev = false;
            for (int col : pi) {
                bool cur = A[r][col - 1];
                if (cur && !prev) blocks++;
                prev = cur;
            }
        }
        return blocks;
    }

    bool isPermutation(const vector<int>& pi) const {
        if (static_cast<int>(pi.size()) != c) return false;
        vector<bool> seen(c + 1, false);
        for (int col : pi) {
            if (col < 1 || col > c || seen[col]) return false;
            seen[col] = true;
        }
        return true;
    }

    void bestInsertion(vector<int>& pi, int& nb1) {
        int n = (int)pi.size();
        bool improved = true;
        while (improved) {
            improved = false;

            for (int j = 1; j < n && !improved; j++) {
                int pj_1 = pi[j - 1];
                int pj = pi[j];
                int pj1 = (j == n - 1) ? c + 1 : pi[j + 1];
                for (int k = 0; k <= j - 2 && !improved; k++) {
                    int pk_1 = (k == 0) ? 0 : pi[k - 1];
                    int pk = pi[k];
                    int d = W[pj_1][pj] + W[pj][pj1] - W[pj_1][pj1] - W[pk_1][pj] - W[pj][pk] + W[pk_1][pk];
                    if (d > 0) {
                        nb1 -= d / 2;
                        for (int i = j; i > k; i--) pi[i] = pi[i - 1];
                        pi[k] = pj;
                        improved = true;
                    }
                }
            }

            for (int j = 0; j < n - 1 && !improved; j++) {
                int pj_1 = (j == 0) ? 0 : pi[j - 1];
                int pj = pi[j];
                int pj1 = pi[j + 1];
                for (int k = j + 2; k <= n && !improved; k++) {
                    int pk_1 = pi[k - 1];
                    int pk = (k == n) ? c + 1 : pi[k];
                    int d = W[pj_1][pj] + W[pj][pj1] - W[pj_1][pj1] - W[pk_1][pj] - W[pj][pk] + W[pk_1][pk];
                    if (d > 0) {
                        nb1 -= d / 2;
                        for (int i = j; i < k - 1; i++) pi[i] = pi[i + 1];
                        pi[k - 1] = pj;
                        improved = true;
                    }
                }
            }
        }
    }

    void submatricesMoving(vector<int>& pi, int& nb1, int p, int callIndex) {
        int n = (int)pi.size();
        if (p <= 0 || n < 4) return;

        vector<pair<int, int>> adjDist;
        for (int j = 0; j < n - 1; j++) adjDist.push_back({W[pi[j]][pi[j + 1]], j});
        sort(adjDist.rbegin(), adjDist.rend());

        vector<bool> cut(n - 1, false);
        int chosen = 0;
        for (auto& [d, j] : adjDist) {
            if (chosen >= p) break;
            if (j == 0 || j == n - 2) continue;
            if (j > 0 && cut[j - 1]) continue;
            if (j < n - 2 && cut[j + 1]) continue;
            cut[j] = true;
            chosen++;
        }
        if (chosen == 0) return;

        vector<vector<int>> subs;
        int start = 0;
        for (int j = 0; j < n - 1; j++) {
            if (cut[j]) {
                subs.push_back({pi.begin() + start, pi.begin() + j + 1});
                start = j + 1;
            }
        }
        subs.push_back({pi.begin() + start, pi.end()});
        int nsubs = (int)subs.size();

        int maxW = 0;
        for (auto& row : W)
            for (int v : row) maxW = max(maxW, v);
        int BIGM = maxW * (2 * nsubs + 2) + 1;

        // City 0 is the depot; cities 2s+1 / 2s+2 are the ends of sub-matrix s,
        // joined by a zero-cost (forced) edge.
        int nCities = 1 + 2 * nsubs;
        vector<int> cityCol(nCities);
        cityCol[0] = 0;
        for (int s = 0; s < nsubs; s++) {
            cityCol[1 + 2 * s] = subs[s].front();
            cityCol[1 + 2 * s + 1] = subs[s].back();
        }

        vector<vector<int>> D(nCities, vector<int>(nCities, 0));
        for (int i = 0; i < nCities; i++)
            for (int j = 0; j < nCities; j++) {
                if (i == j) continue;
                bool forced = (i > 0 && j > 0 && (i - 1) / 2 == (j - 1) / 2);
                D[i][j] = forced ? 0 : BIGM + W[cityCol[i]][cityCol[j]];
            }

        vector<int> tour = tsp.solve(D, callIndex);

        vector<int> newPi;
        newPi.reserve(n);
        for (int step = 1; step < nCities; step++) {
            int idx = tour[step % nCities];
            int s = (idx - 1) / 2;
            if (idx % 2 == 1) {
                for (int col : subs[s]) newPi.push_back(col);
            } else {
                for (int i = (int)subs[s].size() - 1; i >= 0; i--) newPi.push_back(subs[s][i]);
            }
            step++;  // skip the partner end of the same sub-matrix
        }

        // A tour that breaks a forced edge decodes into duplicates/omissions;
        // such a move is rejected rather than scored.
        if (!isPermutation(newPi)) {
            rejectedTours++;
            return;
        }
        int nb = countBlocks(newPi);
        if (nb < nb1) {
            nb1 = nb;
            pi = newPi;
        }
    }

    void twoOpt(vector<int>& pi, const vector<vector<int>>& Wp) {
        int n = (int)pi.size();
        bool improved = true;
        while (improved) {
            improved = false;
            for (int i = 0; i < n - 1 && !improved; i++)
                for (int j = i + 2; j < n && !improved; j++) {
                    int a = (i == 0) ? 0 : pi[i - 1];
                    int b = pi[i];
                    int cc = pi[j];
                    int d = (j == n - 1) ? c + 1 : pi[j + 1];
                    if (Wp[a][cc] + Wp[b][d] < Wp[a][b] + Wp[cc][d]) {
                        reverse(pi.begin() + i, pi.begin() + j + 1);
                        improved = true;
                    }
                }
        }
    }

    void perturb(vector<int>& pi, double phi, int V, mt19937& rng) {
        int N = c + 2;
        vector<vector<int>> Wp = W;
        uniform_real_distribution<double> uni(0.0, 1.0);
        for (int i = 0; i < N; i++)
            for (int j = i + 1; j < N; j++) {
                double x = uni(rng);
                int e = Wp[i][j];
                if (x <= phi && e >= V)
                    e -= V;
                else if (x > 1.0 - phi)
                    e += V;
                Wp[i][j] = Wp[j][i] = e;
            }
        twoOpt(pi, Wp);
    }

    vector<int> initialSolution() {
        vector<int> tour = tsp.solve(W, 0);
        vector<int> pi;
        pi.reserve(c);
        for (int node : tour)
            if (node != 0 && node != c + 1) pi.push_back(node);
        if (!isPermutation(pi)) throw runtime_error("initial TSP tour is not a permutation of the columns");
        return pi;
    }
};

struct Result {
    vector<int> bestPi;
    int bestBlocks = 0, initialBlocks = 0, iterationsDone = 0;
    double timeToBest = 0.0;
    string stopReason = "iterations";
    vector<HistoryEntry> history;
};

static Result solve(CBM& cbm, const Options& opt) {
    Result res;
    mt19937 rng(cbm_derive_seed(opt.seed, UINT64_MAX));  // perturbation stream, disjoint from TSP call indices
    cbm.buildW();

    vector<int> piStar = cbm.initialSolution();
    int bStar = cbm.countBlocks(piStar);
    res.initialBlocks = bStar;
    res.timeToBest = elapsedSeconds();
    res.history.push_back({0, bStar, res.timeToBest});

    vector<int> pi = piStar;
    int nb1 = bStar;
    const int c = cbm.c;

    auto getP = [&](int iter) -> int {
        if (iter < 6) return max(1, c / 5);
        if (iter < 12) return max(1, c / 10);
        return max(1, c / 20);
    };

    for (int iter = 0; iter < opt.iterations; iter++) {
        // Soft time limit, checked between iterations.
        if (elapsedSeconds() >= opt.timeLimit) {
            res.stopReason = "time_limit";
            break;
        }

        cbm.bestInsertion(pi, nb1);
        cbm.submatricesMoving(pi, nb1, getP(iter), iter + 1);

        if (nb1 < bStar) {
            bStar = nb1;
            piStar = pi;
            res.timeToBest = elapsedSeconds();
            res.history.push_back({iter + 1, bStar, res.timeToBest});
        }

        if (iter % 2 == 0) pi = piStar;
        cbm.perturb(pi, opt.phi, opt.V, rng);
        nb1 = cbm.countBlocks(pi);
        res.iterationsDone = iter + 1;
    }

    if (!cbm.isPermutation(piStar) || cbm.countBlocks(piStar) != bStar) throw runtime_error("best solution failed validation");
    res.bestPi = piStar;
    res.bestBlocks = bStar;
    return res;
}

static string jsonEscape(const string& s) {
    string out;
    for (char ch : s) {
        if (ch == '"' || ch == '\\') out += '\\';
        out += ch;
    }
    return out;
}

static void writeJson(ostream& out, const Options& opt, const CBM& cbm, const TSPSolver& tsp, const Result& res) {
    out << fixed << setprecision(6);
    out << "{\n  \"method\": \"ILS\",\n  \"tsp_backend\": \"" << (opt.solver == Solver::LKH ? "lkh" : "linkern") << "\",\n";
    out << "  \"instance\": \"" << jsonEscape(opt.filePath) << "\",\n  \"rows\": " << cbm.l << ",\n  \"cols\": " << cbm.c << ",\n";
    out << "  \"seed\": " << opt.seed << ",\n";
    out << "  \"parameters\": {\"iterations\": " << opt.iterations << ", \"phi\": " << opt.phi << ", \"V\": " << opt.V
        << ", \"lkh_runs\": " << opt.runs << ", \"time_limit_s\": " << opt.timeLimit << "},\n";
    out << "  \"initial_blocks\": " << res.initialBlocks << ",\n  \"best_blocks\": " << res.bestBlocks << ",\n  \"validated\": true,\n";
    out << "  \"iterations_completed\": " << res.iterationsDone << ",\n  \"stop_reason\": \"" << res.stopReason << "\",\n";
    out << "  \"time_limit_reached\": " << (res.stopReason == "time_limit" ? "true" : "false") << ",\n";
    out << "  \"elapsed_s\": " << elapsedSeconds() << ",\n  \"time_to_best_s\": " << res.timeToBest << ",\n";
    out << "  \"perturbation_seed\": " << cbm_derive_seed(opt.seed, UINT64_MAX) << ",\n";
    out << "  \"tsp\": {\"calls\": " << tsp.calls << ", \"time_s\": " << tsp.seconds << ", \"time_limit_hits\": " << tsp.timeLimitHits
        << ", \"invalid_tours\": " << cbm.rejectedTours << ", \"seeds\": [";
    for (size_t i = 0; i < tsp.seeds.size(); i++) out << (i ? ", " : "") << tsp.seeds[i];
    out << "]},\n  \"history\": [";
    for (size_t i = 0; i < res.history.size(); i++)
        out << (i ? ", " : "") << "{\"iteration\": " << res.history[i].iteration << ", \"blocks\": " << res.history[i].blocks
            << ", \"elapsed_s\": " << res.history[i].elapsed << "}";
    out << "],\n  \"permutation\": [";
    for (size_t i = 0; i < res.bestPi.size(); i++) out << (i ? ", " : "") << res.bestPi[i];
    out << "]\n}\n";
}

// Atomic: write to a temporary file, then rename over the target.
static void writeResult(const Options& opt, const CBM& cbm, const TSPSolver& tsp, const Result& res) {
    if (opt.outputPath.empty()) {
        writeJson(cout, opt, cbm, tsp, res);
        return;
    }
    const string tmp = opt.outputPath + ".tmp";
    {
        ofstream out(tmp);
        writeJson(out, opt, cbm, tsp, res);
        out.flush();
        if (!out) throw runtime_error("cannot write " + tmp);
    }
    fs::rename(tmp, opt.outputPath);
}

static Options parseArgs(int argc, char* argv[]) {
    Options opt;
    auto value = [](const string& arg, const string& key) { return arg.substr(key.size()); };
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        try {
            if (arg.rfind("--filePath=", 0) == 0)
                opt.filePath = value(arg, "--filePath=");
            else if (arg.rfind("--outputPath=", 0) == 0)
                opt.outputPath = value(arg, "--outputPath=");
            else if (arg.rfind("--solverPath=", 0) == 0)
                opt.solverPath = value(arg, "--solverPath=");
            else if (arg.rfind("--workDir=", 0) == 0)
                opt.workDir = value(arg, "--workDir=");
            else if (arg.rfind("--seed=", 0) == 0)
                opt.seed = stoull(value(arg, "--seed="));
            else if (arg.rfind("--timeLimit=", 0) == 0)
                opt.timeLimit = stod(value(arg, "--timeLimit="));
            else if (arg.rfind("--iterations=", 0) == 0)
                opt.iterations = stoi(value(arg, "--iterations="));
            else if (arg.rfind("--phi=", 0) == 0)
                opt.phi = stod(value(arg, "--phi="));
            else if (arg.rfind("--V=", 0) == 0)
                opt.V = stoi(value(arg, "--V="));
            else if (arg.rfind("--runs=", 0) == 0)
                opt.runs = stoi(value(arg, "--runs="));
            else if (arg.rfind("--algorithm=", 0) == 0) {
                string alg = value(arg, "--algorithm=");
                if (alg == "lkh")
                    opt.solver = Solver::LKH;
                else if (alg == "linkern")
                    opt.solver = Solver::LINKERN;
                else
                    throw runtime_error("unknown algorithm '" + alg + "' (use lkh or linkern)");
            } else
                throw runtime_error("unknown option " + arg);
        } catch (const invalid_argument&) {
            throw runtime_error("invalid value in " + arg);
        } catch (const out_of_range&) {
            throw runtime_error("value out of range in " + arg);
        }
    }
    if (opt.filePath.empty()) throw runtime_error("missing --filePath=<instance>");
    if (opt.solverPath.empty()) {
        const char* env = getenv(opt.solver == Solver::LKH ? "LKH_PATH" : "LINKERN_PATH");
        if (env) opt.solverPath = env;
    }
    if (opt.solverPath.empty() || access(opt.solverPath.c_str(), X_OK) != 0)
        throw runtime_error("TSP solver not executable: pass --solverPath=<binary> (or set LKH_PATH / LINKERN_PATH)");
    if (opt.timeLimit <= 0 || opt.iterations < 0) throw runtime_error("--timeLimit must be > 0 and --iterations >= 0");
    return opt;
}

int main(int argc, char* argv[]) {
    g_start = std::chrono::steady_clock::now();
    bool ownWorkDir = false;
    Options opt;
    try {
        opt = parseArgs(argc, argv);
        if (opt.workDir.empty()) {
            string tmpl = (fs::temp_directory_path() / "ils-XXXXXX").string();
            if (!mkdtemp(tmpl.data())) throw runtime_error("cannot create a work directory");
            opt.workDir = tmpl;
            ownWorkDir = true;
        } else {
            fs::create_directories(opt.workDir);
        }

        ifstream in(opt.filePath);
        if (!in) throw runtime_error("cannot open instance file " + opt.filePath);

        TSPSolver tsp(opt);
        CBM cbm(tsp);
        cbm.read(in);
        Result res = solve(cbm, opt);
        writeResult(opt, cbm, tsp, res);
    } catch (const exception& e) {
        cerr << "ILS error: " << e.what() << "\n";
        return 1;
    }
    if (ownWorkDir) fs::remove_all(opt.workDir);
    return 0;
}
