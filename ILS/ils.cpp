#include <bits/stdc++.h>
#include <chrono>

using namespace std;

// ─────────────────────────────────────────────────────────────────────────────
//  ILS for Consecutive Block Minimization
//  Supports two TSP backends, selected at runtime via --algorithm=:
//
//    --algorithm=lkh      LKH   (http://webhotel4.ruc.dk/~keld/research/LKH/)
//    --algorithm=linkern  Linkern (http://www.math.uwaterloo.ca/tsp/concorde.html)
//
//  --filePath=<path>    path to the CBM instance file
//  --outputPath=<path>  path to write the result (default: stdout)
//
//  Input format (1-indexed columns):
//    First line : l (rows)  c (cols)
//    Next l lines: n  e1 e2 ... en
// ─────────────────────────────────────────────────────────────────────────────

static const string LKH_BINARY     = "/home/pedro/cbm/src/LKH3/LKH";
static const string LINKERN_BINARY  = "/home/pedro/cbm/linkern";

// ── Solver enum ───────────────────────────────────────────────────────────────
enum class Solver { LKH, LINKERN };
static Solver g_solver = Solver::LKH;   // set by --algorithm=

static std::chrono::steady_clock::time_point g_start;
static double elapsed_seconds() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now() - g_start).count();
}
static constexpr double MAX_SECONDS = 7200.0;

// ── Output path (empty = stdout) ─────────────────────────────────────────────
static string g_outputPath;

// ── Unique filename prefix (PID + monotone counter) ──────────────────────────
static string tmpPrefix() {
    static int counter = 0;
    return "/tmp/cbm_tsp_" + to_string(getpid()) + "_" + to_string(counter++);
}

// ─────────────────────────────────────────────────────────────────────────────
//  writeTSP  –  TSPLIB EXPLICIT / FULL_MATRIX instance (used by both solvers).
//  dist is a square n×n matrix (0-based).
// ─────────────────────────────────────────────────────────────────────────────
static void writeTSP(const string& tspFile,
                     const string& name,
                     const vector<vector<int>>& dist)
{
    int n = (int)dist.size();
    ofstream f(tspFile);
    f << "NAME : "           << name << "\n"
      << "TYPE : TSP\n"
      << "DIMENSION : "      << n    << "\n"
      << "EDGE_WEIGHT_TYPE : EXPLICIT\n"
      << "EDGE_WEIGHT_FORMAT : FULL_MATRIX\n"
      << "EDGE_WEIGHT_SECTION\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            f << dist[i][j] << (j + 1 < n ? " " : "\n");
    }
    f << "EOF\n";
}

// ─────────────────────────────────────────────────────────────────────────────
//  LKH helpers
// ─────────────────────────────────────────────────────────────────────────────

static void writePAR(const string& parFile,
                     const string& tspFile,
                     const string& tourFile,
                     int runs = 1)
{
    ofstream f(parFile);
    f << "PROBLEM_FILE = "      << tspFile  << "\n"
      << "OUTPUT_TOUR_FILE = "  << tourFile << "\n"
      << "RUNS = "              << runs     << "\n"
      << "TRACE_LEVEL = 0\n";
}

// Parse TSPLIB tour format; returns 0-based node indices.
static vector<int> parseTourLKH(const string& tourFile) {
    ifstream f(tourFile);
    string line;
    bool inSection = false;
    vector<int> tour;
    while (getline(f, line)) {
        if (line.find("TOUR_SECTION") != string::npos) { inSection = true; continue; }
        if (!inSection) continue;
        istringstream iss(line);
        int v;
        while (iss >> v) {
            if (v == -1) goto done;
            tour.push_back(v - 1);   // LKH is 1-based → 0-based
        }
    }
    done:
    return tour;
}

static vector<int> runLKH(const string& tspFile,
                           const string& tourFile,
                           int runs)
{
    string parFile = tspFile + ".par";
    writePAR(parFile, tspFile, tourFile, runs);

    string cmd = LKH_BINARY + " " + parFile + " > /dev/null 2>&1";
    int rc = system(cmd.c_str());

    vector<int> tour;
    if (rc == 0) tour = parseTourLKH(tourFile);

    remove(parFile.c_str());
    return tour;
}

// ─────────────────────────────────────────────────────────────────────────────
//  Linkern helpers
//
//  Linkern output format (one line per tour edge):
//    <n_nodes>  <tour_length>
//    <node_a>  <node_b>  <edge_weight>
//    …
//  Nodes are 0-based.
// ─────────────────────────────────────────────────────────────────────────────

static vector<int> parseTourLinkern(const string& tourFile, int expectedNodes) {
    ifstream f(tourFile);
    if (!f) {
        cerr << "[Linkern] Cannot open tour file: " << tourFile << "\n";
        return {};
    }

    int nNodes, tourLen;
    if (!(f >> nNodes >> tourLen)) {
        cerr << "[Linkern] Failed to read header from " << tourFile << "\n";
        return {};
    }
    if (nNodes != expectedNodes) {
        cerr << "[Linkern] Tour has " << nNodes
             << " nodes but expected " << expectedNodes << "\n";
        return {};
    }

    vector<int> tour;
    tour.reserve(nNodes);
    for (int i = 0; i < nNodes; i++) {
        int a, b, w;
        if (!(f >> a >> b >> w)) {
            cerr << "[Linkern] Unexpected EOF at edge " << i << "\n";
            return {};
        }
        tour.push_back(a);
    }
    return tour;
}

static vector<int> runLinkern(const string& tspFile,
                               const string& tourFile,
                               int expectedNodes)
{
    string cmd = LINKERN_BINARY + " -Q -t 7200 -o " + tourFile + " " + tspFile
                 + " > /dev/null 2>&1";
    int rc = system(cmd.c_str());

    if (rc != 0) {
        cerr << "[Linkern] Process exited with code " << rc << "\n";
        return {};
    }
    return parseTourLinkern(tourFile, expectedNodes);
}

// ─────────────────────────────────────────────────────────────────────────────
//  solveTSP  –  dispatcher
// ─────────────────────────────────────────────────────────────────────────────
static vector<int> solveTSP(const vector<vector<int>>& dist, int runs = 1) {
    int n = (int)dist.size();

    string prefix   = tmpPrefix();
    string tspFile  = prefix + ".tsp";
    string tourFile = prefix + ".tour";

    writeTSP(tspFile, prefix, dist);

    vector<int> tour;
    if (g_solver == Solver::LKH)
        tour = runLKH(tspFile, tourFile, runs);
    else
        tour = runLinkern(tspFile, tourFile, n);

    remove(tspFile.c_str());
    remove(tourFile.c_str());

    // ── Fallback: greedy nearest-neighbour ───────────────────────────────────
    if (tour.empty() || (int)tour.size() != n) {
        cerr << "[Warning] Solver returned wrong tour ("
             << tour.size() << " vs " << n << "). Using greedy NN fallback.\n";
        vector<bool> visited(n, false);
        tour.clear();
        int cur = 0; visited[0] = true; tour.push_back(0);
        for (int step = 1; step < n; step++) {
            int best = -1, bestD = INT_MAX;
            for (int i = 0; i < n; i++)
                if (!visited[i] && dist[cur][i] < bestD)
                    { bestD = dist[cur][i]; best = i; }
            visited[best] = true;
            tour.push_back(best);
            cur = best;
        }
    }

    auto it = find(tour.begin(), tour.end(), 0);
    if (it != tour.end() && it != tour.begin())
        rotate(tour.begin(), it, tour.end());

    return tour;
}

// =============================================================================
struct CBM {
    int l, c;
    vector<vector<bool>> A;
    vector<vector<int>>  W;

    void read(istream& in) {
        in >> l >> c;
        A.assign(l, vector<bool>(c, false));
        for (int i = 0; i < l; i++) {
            int n; in >> n;
            for (int j = 0; j < n; j++) {
                int e; in >> e;
                A[i][e - 1] = true;
            }
        }
    }

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
                    d += bi*(1-bj) + (1-bi)*bj;
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

    void bestInsertion(vector<int>& pi, int& nb1) {
        int n = (int)pi.size();
        bool improved = true;
        while (improved) {
            improved = false;

            for (int j = 1; j < n && !improved; j++) {
                int pj_1 = (j == 0)   ? 0   : pi[j-1];
                int pj   = pi[j];
                int pj1  = (j == n-1) ? c+1 : pi[j+1];
                for (int k = 0; k <= j-2 && !improved; k++) {
                    int pk_1 = (k == 0) ? 0 : pi[k-1];
                    int pk   = pi[k];
                    int d = W[pj_1][pj] + W[pj][pj1] - W[pj_1][pj1]
                          - W[pk_1][pj] - W[pj][pk]  + W[pk_1][pk];
                    if (d > 0) {
                        nb1 -= d/2;
                        for (int i = j; i > k; i--) pi[i] = pi[i-1];
                        pi[k] = pj;
                        improved = true;
                    }
                }
            }

            for (int j = 0; j < n-1 && !improved; j++) {
                int pj_1 = (j == 0) ? 0 : pi[j-1];
                int pj   = pi[j];
                int pj1  = pi[j+1];
                for (int k = j+2; k <= n && !improved; k++) {
                    int pk_1 = pi[k-1];
                    int pk   = (k == n) ? c+1 : pi[k];
                    int d = W[pj_1][pj] + W[pj][pj1] - W[pj_1][pj1]
                          - W[pk_1][pj] - W[pj][pk]  + W[pk_1][pk];
                    if (d > 0) {
                        nb1 -= d/2;
                        for (int i = j; i < k-1; i++) pi[i] = pi[i+1];
                        pi[k-1] = pj;
                        improved = true;
                    }
                }
            }
        }
    }

    void submatricesMoving(vector<int>& pi, int& nb1, int p, int runs = 1) {
        int n = (int)pi.size();
        if (p <= 0 || n < 4) return;

        vector<pair<int,int>> adjDist;
        for (int j = 0; j < n-1; j++)
            adjDist.push_back({W[pi[j]][pi[j+1]], j});
        sort(adjDist.rbegin(), adjDist.rend());

        vector<bool> cut(n-1, false);
        int chosen = 0;
        for (auto& [d, j] : adjDist) {
            if (chosen >= p) break;
            if (j == 0 || j == n-2)     continue;
            if (j > 0   && cut[j-1])    continue;
            if (j < n-2 && cut[j+1])    continue;
            cut[j] = true;
            chosen++;
        }
        if (chosen == 0) return;

        vector<vector<int>> subs;
        int start = 0;
        for (int j = 0; j < n-1; j++) {
            if (cut[j]) {
                subs.push_back({pi.begin()+start, pi.begin()+j+1});
                start = j+1;
            }
        }
        subs.push_back({pi.begin()+start, pi.end()});
        int nsubs = (int)subs.size();

        int maxW = 0;
        for (auto& row : W) for (int v : row) maxW = max(maxW, v);
        int BIGM = maxW * (2*nsubs + 2) + 1;

        int nCities = 1 + 2*nsubs;
        vector<int> cityCol(nCities);
        cityCol[0] = 0;
        for (int s = 0; s < nsubs; s++) {
            cityCol[1 + 2*s]     = subs[s].front();
            cityCol[1 + 2*s + 1] = subs[s].back();
        }

        vector<vector<int>> D(nCities, vector<int>(nCities, 0));
        for (int i = 0; i < nCities; i++)
            for (int j = 0; j < nCities; j++) {
                if (i == j) continue;
                bool forced = false;
                for (int s = 0; s < nsubs; s++)
                    if ((i == 1+2*s && j == 1+2*s+1) ||
                        (j == 1+2*s && i == 1+2*s+1))
                    { forced = true; break; }
                D[i][j] = forced ? 0 : BIGM + W[cityCol[i]][cityCol[j]];
            }

        vector<int> tour = solveTSP(D, runs);

        vector<int> newPi;
        newPi.reserve(n);
        for (int step = 1; step < nCities; step++) {
            int idx = tour[step % nCities];
            for (int s = 0; s < nsubs; s++) {
                if (idx == 1 + 2*s) {
                    for (int col : subs[s]) newPi.push_back(col);
                    step++;
                    break;
                } else if (idx == 1 + 2*s + 1) {
                    for (int i = (int)subs[s].size()-1; i >= 0; i--)
                        newPi.push_back(subs[s][i]);
                    step++;
                    break;
                }
            }
        }

        if ((int)newPi.size() == n) {
            int nb = countBlocks(newPi);
            if (nb < nb1) { nb1 = nb; pi = newPi; }
        }
    }

    void twoOpt(vector<int>& pi, const vector<vector<int>>& Wp) {
        int n = (int)pi.size();
        bool improved = true;
        while (improved) {
            improved = false;
            for (int i = 0; i < n-1 && !improved; i++)
                for (int j = i+2; j < n && !improved; j++) {
                    int a  = (i == 0)   ? 0   : pi[i-1];
                    int b  = pi[i];
                    int cc = pi[j];
                    int d  = (j == n-1) ? c+1 : pi[j+1];
                    if (Wp[a][cc] + Wp[b][d] < Wp[a][b] + Wp[cc][d]) {
                        reverse(pi.begin()+i, pi.begin()+j+1);
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
            for (int j = i+1; j < N; j++) {
                double x = uni(rng);
                int e = Wp[i][j];
                if      (x <= phi && e >= V) e -= V;
                else if (x > 1.0 - phi)      e += V;
                Wp[i][j] = Wp[j][i] = e;
            }
        twoOpt(pi, Wp);
    }

    vector<int> initialSolution(int runs = 1) {
        int N = c + 2;
        vector<int> tour = solveTSP(W, runs);

        if ((int)tour.size() == N) {
            vector<int> pi;
            pi.reserve(c);
            for (int node : tour)
                if (node != 0 && node != c + 1)
                    pi.push_back(node);
            if ((int)pi.size() == c) return pi;
        }

        cerr << "[Warning] Solver failed for initial solution. Using identity.\n";
        vector<int> fallback(c);
        iota(fallback.begin(), fallback.end(), 1);
        return fallback;
    }

    pair<vector<int>, int> solve(int    NBITER  = 24,
                                  double phi    = 0.2,
                                  int    V      = 1,
                                  int    seed   = 42,
                                  int    runs   = 1)
    {
        mt19937 rng(seed);
        buildW();

        vector<int> piStar = initialSolution(runs);
        int bStar = countBlocks(piStar);

        vector<int> pi = piStar;
        int nb1 = bStar;

        auto getP = [&](int iter) -> int {
            if (iter < 6)       return max(1, c/5);
            else if (iter < 12) return max(1, c/10);
            else                return max(1, c/20);
        };

        for (int iter = 0; iter < NBITER; iter++) {
            double elapsed = elapsed_seconds();
            if (elapsed >= MAX_SECONDS) {
                cerr << "[TIMEOUT] " << elapsed << "s elapsed — stopping after "
                    << iter << "/" << NBITER << " iterations.\n";
                break;
            }

            bestInsertion(pi, nb1);
            submatricesMoving(pi, nb1, getP(iter), runs);

            if (nb1 < bStar) { bStar = nb1; piStar = pi; }

            if (iter % 2 == 0) pi = piStar;
            perturb(pi, phi, V, rng);
            nb1 = countBlocks(pi);
        }

        return {piStar, bStar};
    }
};

// ── Write result to file or stdout ───────────────────────────────────────────
static void writeResult(const vector<int>& bestPi, int bestBlocks) {
    ostream* out = &cout;
    ofstream fout;

    if (!g_outputPath.empty()) {
        fout.open(g_outputPath);
        if (!fout) {
            cerr << "Cannot open output file: " << g_outputPath
                 << " — falling back to stdout.\n";
        } else {
            out = &fout;
        }
    }

    *out << "Best number of bco's " << bestBlocks << "\n";
    *out << "Best column permutation (1-indexed):";
    for (int col : bestPi) *out << " " << col;
    *out << "\n";
}

// =============================================================================
static void usage(const char* prog) {
    cerr << "Usage: " << prog
         << " --filePath=<path> --algorithm=<lkh|linkern>"
            " [--outputPath=<path>]"
            " [NBITER] [phi] [V] [seed] [runs]\n";
    exit(1);
}

int main(int argc, char* argv[]) {
    g_start = std::chrono::steady_clock::now();

    string filePath;
    bool   gotFile = false, gotAlgo = false;

    int    NBITER = 24;
    double phi    = 0.2;
    int    V      = 1;
    int    seed   = 42;
    int    runs   = 1;
    int    posIdx = 0;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];

        if (arg.rfind("--filePath=", 0) == 0) {
            filePath = arg.substr(11);
            gotFile  = true;
        } else if (arg.rfind("--outputPath=", 0) == 0) {
            g_outputPath = arg.substr(13);
        } else if (arg.rfind("--algorithm=", 0) == 0) {
            string alg = arg.substr(12);
            if      (alg == "lkh")     { g_solver = Solver::LKH;     gotAlgo = true; }
            else if (alg == "linkern") { g_solver = Solver::LINKERN;  gotAlgo = true; }
            else {
                cerr << "Unknown algorithm '" << alg
                     << "'. Use 'lkh' or 'linkern'.\n";
                usage(argv[0]);
            }
        } else if (arg.rfind("--", 0) == 0) {
            cerr << "Unknown option: " << arg << "\n";
            usage(argv[0]);
        } else {
            switch (posIdx++) {
                case 0: NBITER = stoi(arg); break;
                case 1: phi    = stod(arg); break;
                case 2: V      = stoi(arg); break;
                case 3: seed   = stoi(arg); break;
                case 4: runs   = stoi(arg); break;
                default:
                    cerr << "Unexpected positional argument: " << arg << "\n";
                    usage(argv[0]);
            }
        }
    }

    if (!gotFile) { cerr << "Missing required argument --filePath=\n"; usage(argv[0]); }
    if (!gotAlgo) { cerr << "Missing required argument --algorithm=\n"; usage(argv[0]); }

    ifstream in(filePath);
    if (!in) {
        cerr << "Cannot open instance file: " << filePath << "\n";
        return 1;
    }

    CBM cbm;
    cbm.read(in);

    auto [bestPi, bestBlocks] = cbm.solve(NBITER, phi, V, seed, runs);

    writeResult(bestPi, bestBlocks);

    return 0;
}