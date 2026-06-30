#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include <chrono>
#include <iostream>  // Added for std::ostream
#include <vector>

#include "Metrics.hpp"
#include "Solution.hpp"

struct TrajectoryEntry {
    int cost;
    int iteration;
    int elapsedMs;
    BlockMovement movement;
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