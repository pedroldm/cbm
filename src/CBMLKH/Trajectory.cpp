#include "Trajectory.hpp"

#include <utility>

void Trajectory::record(Solution s, int iteration, int elapsedMs) {
    if (s.cost < bestSolution.cost) {
        bestSolution = s;
    }
    // s.blocksCount is the per-column block histogram of the recorded (best)
    // solution; the caller refreshes it before recording. Moved in since s is a
    // by-value copy that is not used afterwards.
    history.push_back({bestSolution.cost, iteration, elapsedMs, s.blockMovement, std::move(s.blocksCount)});
}