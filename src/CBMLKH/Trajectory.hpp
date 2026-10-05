#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include <chrono>
#include <cstdint>
#include <iostream>  // Added for std::ostream
#include <vector>

#include "Metrics.hpp"
#include "Solution.hpp"

struct TrajectoryEntry {
    int cost;
    int iteration;
    int elapsedMs;
    BlockMovement movement;
    // Per-column new-block contribution of the best solution at this check-in:
    // histogram[i] = number of 1-blocks opened at column i (zerosToOnes seam,
    // onesCount for column 0). Snapshots how dense regions "flatten" across the
    // trajectory so the improvement history can be plotted column-by-column.
    std::vector<int> histogram;
};

inline std::ostream& operator<<(std::ostream& os, const TrajectoryEntry& entry) {
    os << "[Iter: " << entry.iteration << " | Cost: " << entry.cost << " | Time: " << entry.elapsedMs << "ms"
       << " | Move: " << toString(entry.movement) << "]";
    return os;
}

class Trajectory {
   public:
    Solution currentSolution;
    Solution bestSolution;
    std::vector<TrajectoryEntry> history;
    Metrics metrics;
    int index = 0;      // global trajectory index (Config::trajectoryOffset + thread slot)
    uint32_t seed = 0;  // seed of this trajectory's RNG stream

    explicit Trajectory(const Solution& initial) : currentSolution(initial), bestSolution(initial) {}
    void record(Solution s, int iteration, int elapsedMs);
};

inline std::ostream& operator<<(std::ostream& os, const Trajectory& trajectory) {
    os << "================ Trajectory Summary ================\n";

    os << "Current Solution: " << trajectory.currentSolution << "\n";
    os << "Best Solution:    " << trajectory.bestSolution << "\n";

    os << "History Log (" << trajectory.history.size() << " check-ins):\n";

    for (const auto& entry : trajectory.history) {
        os << "  " << entry << "\n";
    }

    os << "====================================================";

    return os;
}

#endif