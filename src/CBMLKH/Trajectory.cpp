#include "Trajectory.hpp"

void Trajectory::record(Solution s, int iteration, int elapsedMs) {
    if (s.cost < bestSolution.cost) {
        bestSolution = s;
    }
    history.push_back({bestSolution.cost, iteration, elapsedMs, s.blockMovement});
}