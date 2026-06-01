#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include <chrono>
#include <vector>

#include "Solution.hpp"

struct TrajectoryEntry {
    int cost;
    int iteration;
    int elapsedMs;
};

class Trajectory {
   public:
    Solution currentSolution;
    Solution bestSolution;
    std::vector<TrajectoryEntry> history;

    explicit Trajectory(const Solution& initial) : currentSolution(initial), bestSolution(initial) {}
    void record(Solution s, int iteration, int elapsedMs);
};

#endif