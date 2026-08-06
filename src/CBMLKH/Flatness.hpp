#ifndef FLATNESS_HPP
#define FLATNESS_HPP

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// How the "flatness" of a per-column block histogram is quantified.
//
// The histogram in question is Solution::blocksCount: entry i holds the number
// of 1-blocks opened at position i of the permutation (onesCount for position 0,
// zerosToOnes(sol[i-1], sol[i]) otherwise), so its entries sum exactly to the
// solution's cost. Flattening it therefore means either shaving down the tall
// columns (VARIANCE, PEAK) or spreading the blocks more evenly regardless of the
// total (ENTROPY, GINI — both scale-invariant, so they measure the histogram's
// *shape* and are not merely a proxy for the cost).
//
// Scoped enum on purpose: the unscoped BlockMovement in Solution.hpp already
// claims the name PEAK at global scope.
enum class FlatnessMeasure { VARIANCE, PEAK, ENTROPY, GINI };

inline const char* toString(FlatnessMeasure measure) {
    switch (measure) {
        case FlatnessMeasure::VARIANCE:
            return "VARIANCE";
        case FlatnessMeasure::PEAK:
            return "PEAK";
        case FlatnessMeasure::ENTROPY:
            return "ENTROPY";
        case FlatnessMeasure::GINI:
            return "GINI";
        default:
            return "UNKNOWN";
    }
}

inline FlatnessMeasure parseFlatnessMeasure(const std::string& value) {
    if (value == "VARIANCE") return FlatnessMeasure::VARIANCE;
    if (value == "PEAK") return FlatnessMeasure::PEAK;
    if (value == "ENTROPY") return FlatnessMeasure::ENTROPY;
    if (value == "GINI") return FlatnessMeasure::GINI;
    throw std::runtime_error("Invalid flatness measure: " + value + ". Expected VARIANCE, PEAK, ENTROPY or GINI.");
}

// Flatness score of a block histogram.
//
// CONVENTION: **lower is always flatter**, whatever the measure, so the stop
// criterion can drive every measure through the same `score < bestScore` test.
// ENTROPY is naturally the other way round (a uniform distribution maximizes it)
// and is therefore returned negated, i.e. in [-1, 0].
//
// The all-zero histogram (a zero-cost solution) scores 0 under every measure,
// which is the minimum each of them attains: a solution with no blocks at all is
// as flat as it gets.
inline double flatnessScore(FlatnessMeasure measure, const std::vector<int>& histogram) {
    const std::size_t n = histogram.size();
    if (n == 0) return 0.0;

    switch (measure) {
        case FlatnessMeasure::VARIANCE: {
            // Spread of the block counts around their mean. Lower = the blocks
            // sit at a more uniform height across the permutation.
            const double total = std::accumulate(histogram.begin(), histogram.end(), 0.0);
            const double mean = total / n;
            double squaredError = 0.0;
            for (int value : histogram) {
                const double deviation = value - mean;
                squaredError += deviation * deviation;
            }
            return squaredError / n;
        }

        case FlatnessMeasure::PEAK:
            // Height of the tallest column: the histogram only counts as flatter
            // once the worst seam in the permutation has been shaved down.
            return static_cast<double>(*std::max_element(histogram.begin(), histogram.end()));

        case FlatnessMeasure::ENTROPY: {
            // Normalized Shannon entropy of the histogram read as a distribution
            // over positions. Scale-invariant: it reacts to blocks being
            // redistributed, not to their total being reduced. Negated so that
            // lower stays flatter.
            const double total = std::accumulate(histogram.begin(), histogram.end(), 0.0);
            if (total <= 0.0 || n == 1) return 0.0;  // n == 1: log(n) == 0, nothing to normalize by
            double entropy = 0.0;
            for (int value : histogram) {
                if (value <= 0) continue;
                const double p = value / total;
                entropy -= p * std::log(p);
            }
            return -(entropy / std::log(static_cast<double>(n)));
        }

        case FlatnessMeasure::GINI: {
            // Gini coefficient of the block counts: 0 when every position carries
            // the same number of blocks, approaching 1 when they all pile onto a
            // single position. Scale-invariant like ENTROPY, but O(c log c).
            const double total = std::accumulate(histogram.begin(), histogram.end(), 0.0);
            if (total <= 0.0) return 0.0;
            std::vector<int> sorted(histogram);
            std::sort(sorted.begin(), sorted.end());
            double weighted = 0.0;
            for (std::size_t i = 0; i < n; i++) weighted += (i + 1.0) * sorted[i];
            return (2.0 * weighted) / (n * total) - (n + 1.0) / n;
        }
    }

    throw std::runtime_error("Unexpected flatness measure in flatnessScore.");
}

#endif  // FLATNESS_HPP
