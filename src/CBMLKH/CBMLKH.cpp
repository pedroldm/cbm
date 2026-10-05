#include "CBMLKH.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <filesystem>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "../common/cbm_seed.h"

using namespace std;

namespace {
// Per-thread PRNG (a shared engine would be a data race under OpenMP). run()
// reseeds it at the start of every trajectory, so its stream depends only on
// the trajectory index, never on which OS thread executes the trajectory.
std::mt19937& threadEngine() {
    static thread_local std::mt19937 engine;
    return engine;
}

// Order-sensitive 64-bit hash of a canonical (sorted) LKH sub-problem.
uint64_t hashColumns(const vector<int>& columns) {
    uint64_t h = cbm_splitmix64(columns.size());
    for (int c : columns) h = cbm_splitmix64(h ^ static_cast<uint64_t>(c));
    return h;
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
    // seed = 0 asks for a nondeterministic run: draw the base once and derive
    // everything from it as usual, so the report still names every seed.
    seedBase = cfg.seed != 0 ? cfg.seed : (static_cast<unsigned long>(random_device{}()) | 1UL);
}

uint32_t CBMLKH::trajectorySeed(int globalIndex) const { return cbm_derive_seed(seedBase, static_cast<uint64_t>(globalIndex)); }

vector<Trajectory> CBMLKH::run() {
    // Indexed slots: the output order is the trajectory index, not the order in
    // which threads happen to finish.
    vector<optional<Trajectory>> slots(cfg.threads);
    exception_ptr failure;

#pragma omp parallel for num_threads(cfg.threads) schedule(static, 1)
    for (int i = 0; i < cfg.threads; i++) {
        // An exception must not escape an OpenMP region (it would terminate the
        // process without a message), so it is carried out and rethrown below.
        try {
            const int globalIndex = cfg.trajectoryOffset + i;
            const uint32_t seed = trajectorySeed(globalIndex);
            threadEngine().seed(seed);

            Solution initial = greedyConstruction();
            completeEval(initial);
            Trajectory trajectory = LKHILS(initial);
            trajectory.index = globalIndex;
            trajectory.seed = seed;
            slots[i] = std::move(trajectory);
        } catch (...) {
#pragma omp critical
            if (!failure) failure = current_exception();
        }
    }
    if (failure) rethrow_exception(failure);

    vector<Trajectory> trajectories;
    trajectories.reserve(cfg.threads);
    for (auto& slot : slots) trajectories.push_back(std::move(*slot));
    return trajectories;
}

Trajectory CBMLKH::LKHILS(Solution& initial) {
    Trajectory trajectory(initial);
    AdaptiveParameters adaptive(cfg);

    Metrics& metrics = trajectory.metrics;
    metrics.initialCost = initial.cost;
    metrics.bestCost = trajectory.bestSolution.cost;

    int iterationsWithoutImprovement = 0;

    // Histogram stop criterion: count the consecutive iterations that fail to
    // make the incumbent's block histogram any flatter (under the configured
    // FlatnessMeasure) and abandon the trajectory once cfg.histogramStopInterval
    // of them pile up. Only the incumbent is measured, since a rejected neighbor
    // is discarded and never becomes the search state; that makes this strictly
    // stronger than plain cost stagnation, because an improving move that lowers
    // the cost without flattening the histogram does not reset the counter.
    const bool histogramStopEnabled = cfg.histogramStopInterval > 0;
    int iterationsWithoutFlattening = 0;
    if (histogramStopEnabled) {
        // greedyConstruction only sizes blocksCount; fill it before scoring.
        countBlocksPerColumn(trajectory.bestSolution);
        metrics.bestFlatness = flatnessScore(cfg.histogramStopMeasure, trajectory.bestSolution.blocksCount);
        metrics.initialFlatness = metrics.bestFlatness;
    }

    auto start = chrono::steady_clock::now();
    auto deadline = start + chrono::seconds(cfg.maxTime);

    metrics.stopReason = "maxIterations";
    int i = 0;
    for (; i < cfg.maxIterations; i++) {
        auto now = chrono::steady_clock::now();
        if (now >= deadline) {
            metrics.stopReason = "maxTime";
            break;
        }

        if (histogramStopEnabled && iterationsWithoutFlattening >= cfg.histogramStopInterval) {
            metrics.histogramStop = true;
            metrics.stopReason = "histogram";
            break;
        }

        metrics.neighborBiasHistory.push_back(adaptive.neighborBias);

        int bestBefore = trajectory.bestSolution.cost;
        Solution neighbor = ILSNeighbor(trajectory.currentSolution, adaptive, metrics);

        if (neighbor.cost < trajectory.bestSolution.cost) {
            validator.validate(neighbor.sol, neighbor.cost);
            // ILSNeighbor filled neighbor.blocksCount for the pre-move solution;
            // applyLKH/reinsertBlock then reordered neighbor.sol, so recompute the
            // per-column histogram against the final permutation before recording.
            countBlocksPerColumn(neighbor);
            trajectory.record(neighbor, i, chrono::duration_cast<chrono::milliseconds>(now - start).count());
            // Always continue the search from the best solution found so far.
            trajectory.currentSolution = trajectory.bestSolution;
            metrics.acceptedMoves++;
            OperatorStats& op = metrics.opFor(neighbor.blockMovement);
            op.improvements++;
            op.totalImprovement += bestBefore - neighbor.cost;
            iterationsWithoutImprovement = 0;
            adaptive.reset(cfg);

            // countBlocksPerColumn above already refreshed the histogram of the
            // new incumbent, so scoring it costs only the measure itself.
            if (histogramStopEnabled) {
                double flatness = flatnessScore(cfg.histogramStopMeasure, neighbor.blocksCount);
                if (flatness < metrics.bestFlatness) {
                    metrics.bestFlatness = flatness;
                    metrics.flattenings++;
                    iterationsWithoutFlattening = 0;
                } else {
                    iterationsWithoutFlattening++;
                }
            }
        } else {
            metrics.rejectedMoves++;
            iterationsWithoutImprovement++;
            // The incumbent did not move, so it did not flatten either.
            if (histogramStopEnabled) iterationsWithoutFlattening++;
            if (iterationsWithoutImprovement % cfg.adaptationInterval == 0) {
                adaptive.diversify(cfg);
                metrics.diversifications++;
            }
        }
    }

    metrics.iterations = i;
    metrics.iterationsWithoutFlattening = iterationsWithoutFlattening;
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
    auto lkhStart = chrono::steady_clock::now();

    // LKH sees the sub-problem in canonical (sorted) column order and with a
    // seed derived from its contents, so its answer is a function of the column
    // set alone. That is what the cache key already assumes, and it makes the
    // result independent of which trajectory computed (and cached) it first.
    vector<int> subsegment(s.sol.begin() + cr.start, s.sol.begin() + cr.end + 1);
    sort(subsegment.begin(), subsegment.end());

    metrics.recordBlockSize(static_cast<int>(subsegment.size()));
    metrics.lkhCalls++;

    vector<int> lkhSolution;
    if (!lkhCache->get(subsegment, lkhSolution)) {
        LKHResult result = lkhWrapper.run(subsegment, cfg.lkhMaxTime, cbm_derive_seed(seedBase, hashColumns(subsegment)));
        lkhSolution = std::move(result.tour);
        lkhCache->put(subsegment, lkhSolution);
        metrics.lkhCacheMisses++;
        if (result.timeLimitHit) metrics.lkhTimeLimitHits++;
    }

    metrics.lkhTimeMs += chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - lkhStart).count();

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

    // --- Delta evaluation setup (must read s.sol before it is rebuilt below) ---
    // Internal seam cost of the block in its ORIGINAL order (as it currently sits
    // in s.sol). Needed to derive restFullCost from the old cost.
    int blockOrigInternal = 0;
    for (int i = cr.start + 1; i <= cr.end; i++) blockOrigInternal += columnStore.zerosToOnes(s.sol[i - 1], s.sol[i]);

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

    // restFullCost: cost of `rest` as a standalone solution, derived in O(1) (plus
    // the O(blockLen) blockOrigInternal above) from the old cost s.cost by removing
    // the carved block's contributions and rejoining the A|C seam. This is the
    // delta-eval counterpart of a full O(cols) rescan of `rest`.
    bool hasLeft = cr.start > 0;        // a column precedes the block (prefix A nonempty)
    bool hasRight = cr.end < cols - 1;  // a column follows the block (suffix C nonempty)
    int leftContribution = hasLeft ? columnStore.zerosToOnes(s.sol[cr.start - 1], s.sol[cr.start])
                                   : columnStore.onesCount(s.sol[cr.start]);  // block was the prefix -> head term
    int rightContribution = hasRight ? columnStore.zerosToOnes(s.sol[cr.end], s.sol[cr.end + 1]) : 0;
    int joinContribution;
    if (hasLeft && hasRight)
        joinContribution = columnStore.zerosToOnes(s.sol[cr.start - 1], s.sol[cr.end + 1]);  // A now meets C directly
    else if (!hasLeft && hasRight)
        joinContribution = columnStore.onesCount(s.sol[cr.end + 1]);  // rest[0] becomes the new head
    else
        joinContribution = 0;  // block was a suffix (rest == A): nothing rejoined
    int restFullCost = s.cost - blockOrigInternal - leftContribution - rightContribution + joinContribution;

    // Internal seam cost of the LKH-reordered block.
    int blockInternal = 0;
    for (int i = 1; i < blockLen; i++) blockInternal += columnStore.zerosToOnes(block[i - 1], block[i]);

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

    // Delta evaluation: new cost = rest cost + block-internal cost + best boundary
    // delta. Replaces a full O(cols) completeEval with O(blockLen) work.
    int deltaCost = restFullCost + blockInternal + bestDelta;
#ifdef DELTA_EVAL_VERIFY
    completeEval(s);  // authoritative recompute; the delta result must match it exactly
    if (s.cost != deltaCost) throw runtime_error("Delta eval mismatch: delta=" + to_string(deltaCost) + " full=" + to_string(s.cost));
#else
    s.cost = deltaCost;
#endif
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

    // The two regions are each bounded by maxSegmentSize, but the *gap* between
    // them is not: without this clamp the span can cover nearly the whole
    // permutation, and the sub-problem handed to LKH (an explicit full distance
    // matrix, so quadratic in the segment length) blows up. Bound the merged
    // segment at two maximal regions; when the span exceeds that, keep the
    // higher-scoring region whole and spend the remaining budget extending
    // toward the other one.
    const int maxMergeSize = min(cols, 2 * adaptive.maxSegmentSize);
    if (mergeEnd - mergeStart + 1 > maxMergeSize) {
        if (first.score >= second.score) {
            mergeEnd = mergeStart + maxMergeSize - 1;
        } else {
            mergeStart = mergeEnd - maxMergeSize + 1;
        }
    }

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

    // Ties broken by column id: `remaining` is unordered, so without it the
    // ranking would depend on the hash table's iteration order.
    sort(candidates.begin(), candidates.end(),
         [](const auto& a, const auto& b) { return get<0>(a) != get<0>(b) ? get<0>(a) > get<0>(b) : get<1>(a) < get<1>(b); });

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
    // minSegmentScore is a *relative* threshold: a window qualifies when its
    // average block density is at least minSegmentScore times the density of the
    // permutation as a whole. That keeps the knob portable across instances --
    // the raw score is a block count, so it scales with the row count and with
    // how good the current solution is, and any absolute threshold would mean
    // something different on a 100x200 instance than on a 5000x40000 one, or
    // even at different points of the same run. Read it as "how much denser than
    // average a window must be": 1.0 = at least average, 2.0 = twice as dense.
    const double meanDensity = static_cast<double>(prefix[columnCount]) / columnCount;
    const double requiredDensity = adaptive.minSegmentScore * meanDensity;

    // Best window seen regardless of the score threshold, used as the fallback
    // below.
    CandidateRegion bestBelowThreshold{-1, -1, -1.0};

    for (int left = columnCount - 1; left >= 0; left--) {
        int maxRight = min(columnCount - 1, left + adaptive.maxSegmentSize - 1);
        for (int right = maxRight; right >= left + cfg.minSegmentSize - 1; right--) {
            int size = right - left + 1;
            int blockSum = prefix[right + 1] - prefix[left];

            double averageDensity = static_cast<double>(blockSum) / size;
            // Ranking still uses the raw block sum, so wider/denser windows keep
            // outranking narrow ones exactly as before; only the admission test
            // below changed from absolute to relative.
            double score = averageDensity * size;

            if (averageDensity >= requiredDensity) {
                segments.push_back({left, right, score});
            } else if (score > bestBelowThreshold.score) {
                bestBelowThreshold = {left, right, score};
            }
        }
    }

    // No window cleared minSegmentScore: fall back to the densest one that
    // exists rather than leaving the pool empty. The threshold is a preference,
    // not a feasibility constraint, and an empty pool used to throw from inside
    // the OpenMP region in run() — which is not catchable there, so it aborted
    // the process. That is reachable in normal operation, not just from a
    // misconfiguration: `score` is the window's raw block sum, so every score
    // falls as the solution improves, and a threshold that was easily met at
    // iteration 0 can exclude every window later on.
    if (segments.empty() && bestBelowThreshold.start >= 0) {
        segments.push_back(bestBelowThreshold);
        return segments;
    }

    stable_sort(segments.begin(), segments.end(), [](const auto& a, const auto& b) { return a.score > b.score; });

    return segments;
}

vector<CandidateRegion> CBMLKH::findPeakColumns(Solution& s, const AdaptiveParameters& adaptive) {
    vector<CandidateRegion> peaks;
    peaks.reserve(cols);

    // The window spans [i - halfSize, i + halfSize], i.e. 2 * halfSize + 1
    // columns; halving maxSegmentSize - 1 keeps it at maxSegmentSize at most.
    int halfSize = (adaptive.maxSegmentSize - 1) / 2;

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

    // Every column contributes zero blocks, i.e. the solution is already optimal.
    // Hand back one window anyway: the callers would otherwise throw from inside
    // run()'s OpenMP region, which aborts the process instead of propagating.
    if (peaks.empty() && cols > 0) {
        peaks.push_back({0, min(cols - 1, max(0, adaptive.maxSegmentSize - 1)), 0.0});
        return peaks;
    }

    stable_sort(peaks.begin(), peaks.end(), [](const auto& a, const auto& b) { return a.score > b.score; });

    return peaks;
}
