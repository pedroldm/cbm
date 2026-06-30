#ifndef CBMLKH_CONFIG_HPP
#define CBMLKH_CONFIG_HPP

#include <algorithm>
#include <cmath>
#include <string>

#include "Solution.hpp"

struct Config {
    std::string instancePath;

    int threads = 1;
    int maxIterations = 1000;
    int maxTime = 3600;
    int lkhMaxTime = 5;

    // Base values
    double constructionBias = 1.0;
    double neighborBias = 1.0;
    int minSegmentSize = 5;
    double minSegmentScore = 1.0;

    // Segment sizes are configured as fractions of the column count and resolved
    // to the absolute values below by resolveSegmentSizes() once the instance is
    // loaded (e.g. maxSegmentSizeFraction = 0.1 on a 1000-column instance -> 100).
    double maxSegmentSizeFraction = 0.1;
    double maxSegmentSizeUpperBoundFraction = 0.2;

    // Resolved absolute values (filled by resolveSegmentSizes; do not set directly).
    int maxSegmentSize = 0;
    int maxSegmentSizeUpperBound = 0;

    // Adaptation control
    int adaptationInterval = 20;

    // Bounds

    double minSegmentScoreLowerBound = 0.1;

    // Lower bound for neighborBias: adaptation now *decreases* the bias over
    // stagnation (see AdaptiveParameters::diversify), so the relevant clamp is
    // a floor rather than a ceiling.
    double minNeighborBias = 0.1;

    // Multipliers
    double segmentSizeGrowthFactor = 1.2;
    double segmentScoreDecayFactor = 0.9;
    // < 1: each diversification step shrinks neighborBias toward minNeighborBias.
    double neighborBiasDecayFactor = 0.9;

    BlockMovement blockMovement = RANDOM;

    // Turn the fractional segment-size knobs into absolute column counts. Called
    // once the instance is parsed and `columnCount` is known. The upper bound is
    // clamped to be at least the base size, and both are floored at 1.
    void resolveSegmentSizes(int columnCount) {
        maxSegmentSize = std::max(1, static_cast<int>(std::lround(maxSegmentSizeFraction * columnCount)));
        maxSegmentSizeUpperBound = std::max(maxSegmentSize, static_cast<int>(std::lround(maxSegmentSizeUpperBoundFraction * columnCount)));
    }
};

struct AdaptiveParameters {
    int diversificationLevel = 0;

    int maxSegmentSize;
    double minSegmentScore;

    double neighborBias;

    explicit AdaptiveParameters(const Config& cfg)
        : maxSegmentSize(cfg.maxSegmentSize), minSegmentScore(cfg.minSegmentScore), neighborBias(cfg.neighborBias) {}

    void reset(const Config& cfg) {
        diversificationLevel = 0;

        maxSegmentSize = cfg.maxSegmentSize;
        minSegmentScore = cfg.minSegmentScore;

        neighborBias = cfg.neighborBias;
    }

    void diversify(const Config& cfg) {
        diversificationLevel++;

        maxSegmentSize = std::min(cfg.maxSegmentSizeUpperBound, static_cast<int>(maxSegmentSize * cfg.segmentSizeGrowthFactor));

        minSegmentScore = std::max(cfg.minSegmentScoreLowerBound, minSegmentScore * cfg.segmentScoreDecayFactor);

        // neighborBias schedule (decreasing): starting from cfg.neighborBias,
        // each diversification step multiplies by neighborBiasDecayFactor (< 1),
        // i.e. neighborBias_k = max(minNeighborBias, neighborBias_0 * decay^k).
        // A smaller bias flattens the rank-based roulette (weight 1/rank^bias),
        // so candidate selection becomes progressively more exploratory as
        // stagnation grows; reset() restores it to cfg.neighborBias.
        neighborBias = std::max(cfg.minNeighborBias, neighborBias * cfg.neighborBiasDecayFactor);
    }
};

#endif