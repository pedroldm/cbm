#ifndef CBMLKH_CONFIG_HPP
#define CBMLKH_CONFIG_HPP

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "Flatness.hpp"
#include "Solution.hpp"

struct Config {
    std::string instancePath;

    // iRace mode: print only the best cost found (a bare integer on stdout)
    // instead of the full JSON report, which is all a tuner's target-runner
    // consumes. Diagnostics keep going to stderr, so stdout holds just the number.
    bool iRace = false;

    // Where to write the JSON report (atomically). Empty = stdout.
    std::string outputPath;

    // LKH executable and the parent directory for scratch files. Each process
    // creates (and removes) its own private subdirectory under lkhTmpDir, so
    // concurrent processes may share it.
    std::string lkhPath;
    std::string lkhTmpDir = "/tmp";

    // RNG seed. 0 (the default) draws a nondeterministic base seed, i.e. runs
    // are not reproducible. Any other value gives trajectory i the seed
    // cbm_derive_seed(seed, trajectoryOffset + i), and every LKH sub-problem a
    // seed derived from (seed, sub-problem contents); see CBMLKH::applyLKH.
    //
    // Caveat: this pins every random choice, but a run bounded by maxTime (or
    // whose LKH calls hit lkhMaxTime) still executes a load-dependent amount of
    // work, so exact replay additionally requires neither bound to bind.
    unsigned long seed = 0;

    // Global index of this process's first trajectory. Lets one logical set of
    // repetitions be split over several executions (e.g. 10 trajectories as
    // two runs of 5 threads) with exactly the seeds a single run would use.
    int trajectoryOffset = 0;

    int threads = 1;
    int maxIterations = 1000;
    int maxTime = 3600;
    int lkhMaxTime = 5;

    // Base values
    double constructionBias = 1.0;
    double neighborBias = 1.0;
    int minSegmentSize = 5;
    // Relative density threshold for INTERVAL candidate windows, in multiples of
    // the permutation's own mean block density (see findDenseSegments): 1.0
    // admits any window at least as dense as average, 2.0 only windows twice as
    // dense. Being a ratio rather than a block count, it carries across
    // instances of different sizes and densities.
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

    // Histogram stop criterion: halt a trajectory once this many consecutive ILS
    // iterations have passed without the incumbent's block histogram getting
    // flatter (as judged by histogramStopMeasure). 0 disables the criterion, so
    // the run is bounded by maxIterations / maxTime alone.
    int histogramStopInterval = 0;
    FlatnessMeasure histogramStopMeasure = FlatnessMeasure::VARIANCE;

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
    // once the instance is parsed and `columnCount` is known. Both are floored at
    // minSegmentSize (a window narrower than that is never enumerated, so a
    // smaller cap would leave the candidate pools empty) and capped at the column
    // count; the upper bound is additionally clamped to be at least the base size.
    void resolveSegmentSizes(int columnCount) {
        const int floorSize = std::min(minSegmentSize, columnCount);
        auto resolve = [&](double fraction) {
            int size = static_cast<int>(std::lround(fraction * columnCount));
            return std::clamp(size, floorSize, columnCount);
        };
        maxSegmentSize = resolve(maxSegmentSizeFraction);
        maxSegmentSizeUpperBound = std::max(maxSegmentSize, resolve(maxSegmentSizeUpperBoundFraction));
    }

    // Reject configurations that would otherwise fail deep inside the search (or,
    // worse, silently degrade it): a non-positive segment fraction empties every
    // candidate pool, a growth factor below 1 shrinks the segment on
    // diversification, a decay factor above 1 grows the score threshold, etc.
    void validate() const {
        auto require = [](bool ok, const std::string& message) {
            if (!ok) throw std::runtime_error("Invalid configuration: " + message);
        };

        require(!instancePath.empty(), "instancePath must be set.");
        require(threads >= 1, "threads must be >= 1.");
        require(trajectoryOffset >= 0, "trajectoryOffset must be >= 0.");
        require(!lkhPath.empty(), "lkhPath must be set (config key or LKH_PATH environment variable).");
        require(maxIterations >= 1, "maxIterations must be >= 1.");
        require(maxTime >= 1, "maxTime must be >= 1 (seconds).");
        require(lkhMaxTime >= 1, "lkhMaxTime must be >= 1 (seconds).");

        require(constructionBias > 0.0, "constructionBias must be > 0.");
        require(neighborBias > 0.0, "neighborBias must be > 0.");
        require(minNeighborBias > 0.0 && minNeighborBias <= neighborBias, "minNeighborBias must be in (0, neighborBias].");

        require(minSegmentSize >= 2, "minSegmentSize must be >= 2.");
        require(maxSegmentSizeFraction > 0.0 && maxSegmentSizeFraction <= 1.0, "maxSegmentSize must be a fraction in (0, 1].");
        require(maxSegmentSizeUpperBoundFraction >= maxSegmentSizeFraction && maxSegmentSizeUpperBoundFraction <= 1.0,
                "maxSegmentSizeUpperBound must be a fraction in [maxSegmentSize, 1].");

        require(minSegmentScore > 0.0, "minSegmentScore must be > 0.");
        require(minSegmentScoreLowerBound > 0.0 && minSegmentScoreLowerBound <= minSegmentScore,
                "minSegmentScoreLowerBound must be in (0, minSegmentScore].");

        require(adaptationInterval >= 1, "adaptationInterval must be >= 1.");
        require(histogramStopInterval >= 0, "histogramStopInterval must be >= 0 (0 disables the histogram stop criterion).");
        require(segmentSizeGrowthFactor >= 1.0, "segmentSizeGrowthFactor must be >= 1.");
        require(segmentScoreDecayFactor > 0.0 && segmentScoreDecayFactor <= 1.0, "segmentScoreDecayFactor must be in (0, 1].");
        require(neighborBiasDecayFactor > 0.0 && neighborBiasDecayFactor <= 1.0, "neighborBiasDecayFactor must be in (0, 1].");
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