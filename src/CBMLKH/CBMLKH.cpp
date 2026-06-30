#include "CBMLKH.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <functional>
#include <numeric>
#include <random>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <utility>

using namespace std;

namespace {
// Per-thread PRNG. A single shared std::mt19937 is not safe under the OpenMP
// parallel region in run(): concurrent access is a data race that corrupts the
// engine state. Each thread lazily constructs its own engine, seeded distinctly
// from random_device combined with the thread id.
std::mt19937& threadEngine() {
    static thread_local std::mt19937 engine(std::random_device{}() ^ static_cast<unsigned>(std::hash<std::thread::id>{}(std::this_thread::get_id())));
    return engine;
}
}  // namespace

CBMLKH::CBMLKH(const Config& cfg, shared_ptr<LKHCache> cache)
    : cfg(cfg),
      validator(cfg.instancePath),
      columnStore(cfg.instancePath),
      rows(columnStore.rows()),
      cols(columnStore.cols()),
      lkhWrapper(columnStore),
      lkhCache(cache) {
    instanceName = filesystem::path(cfg.instancePath).filename().string();
    this->cfg.resolveSegmentSizes(cols);
}

vector<Trajectory> CBMLKH::run() {
    vector<Trajectory> trajectories;
    trajectories.reserve(cfg.threads);
    vector<Solution> initialSolutions(cfg.threads);

#pragma omp parallel for num_threads(cfg.threads)
    for (int i = 0; i < cfg.threads; i++) {
        initialSolutions[i] = greedyConstruction();
        completeEval(initialSolutions[i]);
        Trajectory trajectory = LKHILS(initialSolutions[i]);
        // Concurrent push_back on a shared vector is a data race; serialize the
        // hand-off (the expensive LKHILS work above stays parallel).
#pragma omp critical
        trajectories.push_back(std::move(trajectory));
    }
    return trajectories;
}

Trajectory CBMLKH::LKHILS(Solution& initial) {
    Trajectory trajectory(initial);
    AdaptiveParameters adaptive(cfg);

    Metrics& metrics = trajectory.metrics;
    metrics.initialCost = initial.cost;
    metrics.bestCost = trajectory.bestSolution.cost;

    int iterationsWithoutImprovement = 0;

    auto start = chrono::steady_clock::now();
    auto deadline = start + chrono::seconds(cfg.maxTime);

    int i = 0;
    for (; i < cfg.maxIterations; i++) {
        auto now = chrono::steady_clock::now();
        if (now >= deadline) break;

        metrics.neighborBiasHistory.push_back(adaptive.neighborBias);

        int bestBefore = trajectory.bestSolution.cost;
        Solution neighbor = ILSNeighbor(trajectory.currentSolution, adaptive, metrics);

        if (neighbor.cost < trajectory.bestSolution.cost) {
            validator.validate(neighbor.sol, neighbor.cost);
            trajectory.record(neighbor, i, chrono::duration_cast<chrono::milliseconds>(now - start).count());
            // Always continue the search from the best solution found so far.
            trajectory.currentSolution = trajectory.bestSolution;
            metrics.acceptedMoves++;
            OperatorStats& op = metrics.opFor(neighbor.blockMovement);
            op.improvements++;
            op.totalImprovement += bestBefore - neighbor.cost;
            iterationsWithoutImprovement = 0;
            adaptive.reset(cfg);
        } else {
            metrics.rejectedMoves++;
            iterationsWithoutImprovement++;
            if (iterationsWithoutImprovement % cfg.adaptationInterval == 0) {
                adaptive.diversify(cfg);
                metrics.diversifications++;
            }
        }
    }

    metrics.iterations = i;
    metrics.bestCost = trajectory.bestSolution.cost;
    metrics.finalCost = trajectory.currentSolution.cost;
    metrics.elapsedMs = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start).count();

    return trajectory;
}

Solution CBMLKH::ILSNeighbor(const Solution& s, const AdaptiveParameters& adaptive, Metrics& metrics) {
    Solution neighbor = s;

    countBlocksPerColumn(neighbor);

    BlockMovement movement = cfg.blockMovement;
    if (movement == RANDOM) {
        uniform_int_distribution<int> coin(0, 2);
        switch (coin(threadEngine())) {
            case 0:
                movement = PEAK;
                break;
            case 1:
                movement = INTERVAL;
                break;
            case 2:
                movement = MERGE;
                break;
        }
    }

    CandidateRegion cr;
    switch (movement) {
        case PEAK:
            cr = choosePeakRegion(neighbor, adaptive);
            break;
        case INTERVAL:
            cr = chooseIntervalRegion(neighbor, adaptive);
            break;
        case MERGE:
            cr = chooseMergeRegion(neighbor, adaptive);
            break;
        default:
            throw runtime_error("Unexpected block movement in ILSNeighbor.");
    }
    neighbor.blockMovement = movement;
    metrics.opFor(movement).applications++;

    applyLKH(neighbor, cr, metrics);

    return neighbor;
}

void CBMLKH::applyLKH(Solution& s, const CandidateRegion& cr, Metrics& metrics) {
    vector<int> subsegment(s.sol.begin() + cr.start, s.sol.begin() + cr.end + 1);

    metrics.recordBlockSize(static_cast<int>(subsegment.size()));
    metrics.lkhCalls++;

    vector<int> lkhSolution;
    if (!lkhCache->get(subsegment, lkhSolution)) {
        lkhSolution = lkhWrapper.run(subsegment, instanceName, cfg.lkhMaxTime);
        lkhCache->put(subsegment, lkhSolution);
        metrics.lkhCacheMisses++;
    }

    // Best reinsertion: instead of writing the optimized block back into its
    // original slot, search for the position whose neighbors integrate best.
    reinsertBlock(s, cr, lkhSolution);
}

// Reinsert the LKH-optimized `block` into the best position of the solution.
//
// Rationale / similarity metric:
//   The CBM cost is the number of 1-blocks across rows, accumulated as a sum of
//   seam costs zerosToOnes(a, b) = #rows where column a holds 0 and column b
//   holds 1 (i.e. how many *new* 1-blocks open when b directly follows a). A low
//   seam cost therefore means b's 1-pattern is "absorbed" by a's — the two
//   columns are highly compatible/similar at that boundary. We use this seam
//   cost as the (well-defined, objective-aligned) similarity measure: the best
//   insertion gap is the one whose surrounding columns are most compatible with
//   the block's boundary columns (its first column `bf` and last column `bb`).
//
//   This is strictly stronger than a generic Hamming-similarity heuristic
//   because it is exactly the marginal objective contribution of the insertion,
//   so the chosen position can never be worse than leaving the block in place
//   (the original gap is among the candidates). It runs in O(cols): the block's
//   internal cost and the rest's internal cost are invariant across gaps, so we
//   only minimize the boundary delta, evaluated in O(1) per candidate gap.
void CBMLKH::reinsertBlock(Solution& s, const CandidateRegion& cr, const vector<int>& block) {
    int blockLen = static_cast<int>(block.size());

    // Sequence with the block carved out.
    vector<int> rest;
    rest.reserve(s.sol.size() - blockLen);
    rest.insert(rest.end(), s.sol.begin(), s.sol.begin() + cr.start);
    rest.insert(rest.end(), s.sol.begin() + cr.end + 1, s.sol.end());

    if (rest.empty()) {  // Block spanned the whole permutation: nothing to place against.
        s.sol = block;
        completeEval(s);
        return;
    }

    int bf = block.front();  // block's leading column
    int bb = block.back();   // block's trailing column
    int restSize = static_cast<int>(rest.size());

    // deltaFromRest(g): boundary cost of inserting the block at gap g, relative
    // to the invariant (restCost + block-internal cost). Minimize it.
    auto delta = [&](int g) -> int {
        if (g == 0) {
            // Block becomes the prefix; rest[0] stops being the leading column.
            return columnStore.onesCount(bf) + columnStore.zerosToOnes(bb, rest[0]) - columnStore.onesCount(rest[0]);
        }
        if (g == restSize) {
            // Block appended after the last column.
            return columnStore.zerosToOnes(rest[restSize - 1], bf);
        }
        int a = rest[g - 1];
        int b = rest[g];
        return columnStore.zerosToOnes(a, bf) + columnStore.zerosToOnes(bb, b) - columnStore.zerosToOnes(a, b);
    };

    int bestGap = 0;
    int bestDelta = delta(0);
    for (int g = 1; g <= restSize; g++) {
        int d = delta(g);
        if (d < bestDelta) {
            bestDelta = d;
            bestGap = g;
        }
    }

    s.sol.clear();
    s.sol.insert(s.sol.end(), rest.begin(), rest.begin() + bestGap);
    s.sol.insert(s.sol.end(), block.begin(), block.end());
    s.sol.insert(s.sol.end(), rest.begin() + bestGap, rest.end());

    completeEval(s);
}

size_t CBMLKH::sampleRankIndex(size_t count, double neighborBias) {
    vector<double> weights(count);
    for (size_t rank = 0; rank < count; ++rank) weights[rank] = 1.0 / pow(rank + 1, neighborBias);

    double totalWeight = accumulate(weights.begin(), weights.end(), 0.0);
    uniform_real_distribution<double> dist(0.0, totalWeight);

    double roll = dist(threadEngine());
    double cumulative = 0.0;
    size_t chosen = 0;
    for (; chosen < count; ++chosen) {
        cumulative += weights[chosen];
        if (roll < cumulative) break;
    }

    return min(chosen, count - 1);
}

CandidateRegion CBMLKH::chooseIntervalRegion(Solution& s, const AdaptiveParameters& adaptive) {
    vector<CandidateRegion> segments = findDenseSegments(s, adaptive);
    if (segments.empty()) throw runtime_error("No candidate segments found for INTERVAL neighbor.");

    return segments[sampleRankIndex(segments.size(), adaptive.neighborBias)];
}

CandidateRegion CBMLKH::choosePeakRegion(Solution& s, const AdaptiveParameters& adaptive) {
    vector<CandidateRegion> peaks = findPeakColumns(s, adaptive);
    if (peaks.empty()) throw runtime_error("No peaks found for PEAK neighbor generation.");

    return peaks[sampleRankIndex(peaks.size(), adaptive.neighborBias)];
}

CandidateRegion CBMLKH::chooseMergeRegion(Solution& s, const AdaptiveParameters& adaptive) {
    auto pickRegion = [&]() -> CandidateRegion {
        uniform_int_distribution<int> coin(0, 1);
        vector<CandidateRegion> pool = coin(threadEngine()) ? findDenseSegments(s, adaptive) : findPeakColumns(s, adaptive);
        if (pool.empty()) throw runtime_error("chooseMergeRegion: no candidates available.");
        return pool[sampleRankIndex(pool.size(), adaptive.neighborBias)];
    };

    CandidateRegion first = pickRegion();
    CandidateRegion second = pickRegion();

    // Ensure `first` precedes `second`, then merge into their spanning range
    // (overlapping regions are handled naturally by the union).
    if (first.start > second.start) swap(first, second);

    int mergeStart = first.start;
    int mergeEnd = max(first.end, second.end);
    return {mergeStart, mergeEnd, first.score + second.score};
}

int CBMLKH::completeEval(Solution& s) {
    s.cost = columnStore.onesCount(s.sol[0]);
    for (int col = 1; col < cols; col++) {
        s.cost += columnStore.zerosToOnes(s.sol[col - 1], s.sol[col]);
    }
    return s.cost;
}

Solution CBMLKH::greedyConstruction() {
    uniform_int_distribution<> colDist(0, cols - 1);
    unordered_set<int> remaining;
    Solution s;

    s.cost = 0;
    s.sol.resize(cols);
    s.blocksCount.resize(cols);

    for (int i = 0; i < cols; i++) remaining.insert(i);

    int current = colDist(threadEngine());
    int pos = 0;
    remaining.erase(current);
    s.sol[pos++] = current;

    while (!remaining.empty()) {
        current = nextInsertion(current, remaining);
        s.sol[pos++] = current;
        remaining.erase(current);
    }

    completeEval(s);
    return s;
}

int CBMLKH::nextInsertion(int current, unordered_set<int>& remaining) {
    if (remaining.empty()) throw runtime_error("No remaining columns to insert.");

    // Rank candidates by similarity to `current` (more shared rows first).
    vector<tuple<int, int>> candidates;
    candidates.reserve(remaining.size());
    for (int candidate : remaining) candidates.push_back({rows - columnStore.hamming(current, candidate), candidate});

    sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) { return get<0>(a) > get<0>(b); });

    size_t chosen = sampleRankIndex(candidates.size(), cfg.constructionBias);
    return get<1>(candidates[chosen]);
}

int CBMLKH::countBlocksPerColumn(Solution& s, int start, int end) {
    if (end == -1) end = cols - 1;

    if (start == 0) {
        s.blocksCount[0] = columnStore.onesCount(s.sol[0]);
        start = 1;
    }

    for (int i = start; i <= end; i++) {
        s.blocksCount[i] = columnStore.zerosToOnes(s.sol[i - 1], s.sol[i]);
    }

    return accumulate(s.blocksCount.begin() + start, s.blocksCount.begin() + end + 1, 0);
}

vector<CandidateRegion> CBMLKH::findDenseSegments(Solution& s, const AdaptiveParameters& adaptive) {
    vector<CandidateRegion> segments;

    int columnCount = static_cast<int>(s.blocksCount.size());
    vector<int> prefix(columnCount + 1, 0);
    for (int i = 0; i < columnCount; i++) {
        prefix[i + 1] = prefix[i] + s.blocksCount[i];
    }

    // Blocks are discovered right-to-left: the outer window anchor walks from
    // the last column down to the first, and for each anchor the segment's far
    // end is scanned from its largest admissible value down to the smallest.
    // This enumerates exactly the same set of [left, right] windows (same size
    // and score bounds) as a forward scan, only in reverse order; the final
    // sort makes the candidate set fed to selection identical.
    for (int left = columnCount - 1; left >= 0; left--) {
        int maxRight = min(columnCount - 1, left + adaptive.maxSegmentSize - 1);
        for (int right = maxRight; right >= left + cfg.minSegmentSize - 1; right--) {
            int size = right - left + 1;
            int blockSum = prefix[right + 1] - prefix[left];

            double averageDensity = static_cast<double>(blockSum) / size;
            double score = averageDensity * size;

            if (score >= adaptive.minSegmentScore) {
                segments.push_back({left, right, score});
            }
        }
    }

    sort(segments.begin(), segments.end(), [](const auto& a, const auto& b) { return a.score > b.score; });

    return segments;
}

vector<CandidateRegion> CBMLKH::findPeakColumns(Solution& s, const AdaptiveParameters& adaptive) {
    vector<CandidateRegion> peaks;
    peaks.reserve(cols);

    int halfSize = adaptive.maxSegmentSize / 2;

    // Peaks are discovered right-to-left: scan candidate columns from the last
    // to the first. Only the discovery order changes; each peak's window
    // [i - halfSize, i + halfSize] and score are computed exactly as before, and
    // the trailing sort yields the same ranked candidate set.
    for (int i = cols - 1; i >= 0; i--) {
        if (s.blocksCount[i] > 0) {
            int start = max(0, i - halfSize);
            int end = min(cols - 1, i + halfSize);
            peaks.push_back({start, end, static_cast<double>(s.blocksCount[i])});
        }
    }

    sort(peaks.begin(), peaks.end(), [](const auto& a, const auto& b) { return a.score > b.score; });

    return peaks;
}
