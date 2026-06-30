#ifndef METRICS_HPP
#define METRICS_HPP

#include <climits>
#include <vector>

#include "Solution.hpp"

// Per-operator behavioural counters used for algorithm analysis / benchmarking.
struct OperatorStats {
    long applications = 0;      // times this operator produced a neighbor
    long improvements = 0;      // improving moves attributed to this operator
    long totalImprovement = 0;  // summed cost reduction (best_before - new) over improving moves
};

// Per-trajectory metrics gathered during a single LKHILS run. Aggregated and
// serialized to JSON at the end of main().
struct Metrics {
    int initialCost = 0;  // cost of the constructed solution
    int finalCost = 0;    // cost of currentSolution when the loop ends
    int bestCost = 0;     // best cost found
    int iterations = 0;   // ILS iterations actually executed
    long elapsedMs = 0;   // wall-clock time spent in LKHILS

    long acceptedMoves = 0;  // neighbors that improved on the best
    long rejectedMoves = 0;  // neighbors that did not

    long lkhCalls = 0;        // applyLKH invocations (sub-segment optimizations)
    long lkhCacheMisses = 0;  // applyLKH calls that actually shelled out to LKH

    long diversifications = 0;  // adaptive.diversify() calls

    OperatorStats peak;
    OperatorStats interval;
    OperatorStats merge;

    // Optimized-block size statistics (segment length handed to LKH).
    long blockSizeSum = 0;
    long blockCount = 0;
    int blockSizeMin = INT_MAX;
    int blockSizeMax = 0;

    // neighborBias sampled once per iteration (lets the decreasing schedule be
    // inspected after the fact).
    std::vector<double> neighborBiasHistory;

    OperatorStats& opFor(BlockMovement m) {
        switch (m) {
            case PEAK:
                return peak;
            case INTERVAL:
                return interval;
            default:
                return merge;
        }
    }

    void recordBlockSize(int size) {
        blockSizeSum += size;
        blockCount++;
        if (size < blockSizeMin) blockSizeMin = size;
        if (size > blockSizeMax) blockSizeMax = size;
    }

    double averageBlockSize() const { return blockCount > 0 ? static_cast<double>(blockSizeSum) / blockCount : 0.0; }
    int minBlockSize() const { return blockCount > 0 ? blockSizeMin : 0; }
};

#endif  // METRICS_HPP
