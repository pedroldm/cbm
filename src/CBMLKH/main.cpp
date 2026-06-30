#include <chrono>
#include <exception>
#include <iostream>
#include <memory>

#include "../IO/json.hpp"
#include "ArgsUtil.hpp"
#include "CBMLKH.hpp"
#include "LKHCache.hpp"
#include "LKHWrapper.hpp"
#include "Metrics.hpp"

using namespace std;
using json = nlohmann::json;

static json operatorToJson(const OperatorStats& op) {
    return json{
        {"applications", op.applications},
        {"improvements", op.improvements},
        {"totalImprovement", op.totalImprovement},
    };
}

static json metricsToJson(const Metrics& m) {
    return json{
        {"initialCost", m.initialCost},
        {"bestCost", m.bestCost},
        {"finalCost", m.finalCost},
        {"iterations", m.iterations},
        {"elapsedMs", m.elapsedMs},
        {"acceptedMoves", m.acceptedMoves},
        {"rejectedMoves", m.rejectedMoves},
        {"lkhCalls", m.lkhCalls},
        {"lkhCacheMisses", m.lkhCacheMisses},
        {"diversifications", m.diversifications},
        {"blockSize", {{"average", m.averageBlockSize()}, {"min", m.minBlockSize()}, {"max", m.blockSizeMax}, {"count", m.blockCount}}},
        {"operators", {{"PEAK", operatorToJson(m.peak)}, {"INTERVAL", operatorToJson(m.interval)}, {"MERGE", operatorToJson(m.merge)}}},
        {"neighborBiasHistory", m.neighborBiasHistory},
    };
}

int main(int argc, char* argv[]) {
    if (argc != 2) throw runtime_error("Usage: ./cbmlkh <config_file>");

    Config cfg = ArgsUtil::parseConfigFile(argv[1]);

    LKHWrapper::clearTmpDir();
    auto cache = make_shared<LKHCache>();
    CBMLKH cbmlkh(cfg, cache);

    auto runStart = chrono::steady_clock::now();
    vector<Trajectory> trajectories = cbmlkh.run();
    auto runtimeMs = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - runStart).count();

    // Aggregate per-trajectory metrics and locate the global best.
    int bestIndex = -1;
    OperatorStats aggPeak, aggInterval, aggMerge;
    long totalLkhCalls = 0, totalAccepted = 0, totalRejected = 0;

    json trajectoriesJson = json::array();
    for (size_t t = 0; t < trajectories.size(); t++) {
        const Trajectory& traj = trajectories[t];
        if (bestIndex == -1 || traj.bestSolution.cost < trajectories[bestIndex].bestSolution.cost) {
            bestIndex = static_cast<int>(t);
        }

        const Metrics& m = traj.metrics;
        aggPeak.applications += m.peak.applications;
        aggPeak.improvements += m.peak.improvements;
        aggPeak.totalImprovement += m.peak.totalImprovement;
        aggInterval.applications += m.interval.applications;
        aggInterval.improvements += m.interval.improvements;
        aggInterval.totalImprovement += m.interval.totalImprovement;
        aggMerge.applications += m.merge.applications;
        aggMerge.improvements += m.merge.improvements;
        aggMerge.totalImprovement += m.merge.totalImprovement;
        totalLkhCalls += m.lkhCalls;
        totalAccepted += m.acceptedMoves;
        totalRejected += m.rejectedMoves;

        json history = json::array();
        for (const auto& e : traj.history) {
            history.push_back({{"iteration", e.iteration}, {"cost", e.cost}, {"elapsedMs", e.elapsedMs}, {"move", toString(e.movement)}});
        }

        json tj = metricsToJson(m);
        tj["index"] = t;
        tj["history"] = history;
        trajectoriesJson.push_back(tj);
    }

    long long hits = cache->getHits();
    long long misses = cache->getMisses();
    long long requests = hits + misses;
    double hitRate = requests > 0 ? (100.0 * hits / requests) : 0.0;

    json output = {
        {"instance", {{"name", cbmlkh.instanceName}, {"rows", cbmlkh.rows}, {"cols", cbmlkh.cols}}},
        {"config",
         {{"threads", cfg.threads},
          {"blockMovement", toString(cfg.blockMovement)},
          {"maxIterations", cfg.maxIterations},
          {"maxTime", cfg.maxTime},
          {"lkhMaxTime", cfg.lkhMaxTime},
          {"constructionBias", cfg.constructionBias},
          {"neighborBias", cfg.neighborBias},
          {"minNeighborBias", cfg.minNeighborBias},
          {"neighborBiasDecayFactor", cfg.neighborBiasDecayFactor},
          {"minSegmentSize", cfg.minSegmentSize},
          {"maxSegmentSize", cfg.maxSegmentSize},
          {"maxSegmentSizeUpperBound", cfg.maxSegmentSizeUpperBound},
          {"minSegmentScore", cfg.minSegmentScore},
          {"minSegmentScoreLowerBound", cfg.minSegmentScoreLowerBound},
          {"segmentSizeGrowthFactor", cfg.segmentSizeGrowthFactor},
          {"segmentScoreDecayFactor", cfg.segmentScoreDecayFactor},
          {"adaptationInterval", cfg.adaptationInterval}}},
        {"global",
         {{"bestCost", bestIndex >= 0 ? trajectories[bestIndex].bestSolution.cost : -1},
          {"bestBlockMovement", bestIndex >= 0 ? toString(trajectories[bestIndex].bestSolution.blockMovement) : "NONE"},
          {"bestTrajectory", bestIndex},
          {"runtimeMs", runtimeMs},
          {"totalLkhCalls", totalLkhCalls},
          {"acceptedMoves", totalAccepted},
          {"rejectedMoves", totalRejected},
          {"lkhCache", {{"hits", hits}, {"misses", misses}, {"requests", requests}, {"hitRate", hitRate}}},
          {"operators", {{"PEAK", operatorToJson(aggPeak)}, {"INTERVAL", operatorToJson(aggInterval)}, {"MERGE", operatorToJson(aggMerge)}}}}},
        {"trajectories", trajectoriesJson},
    };

    cout << output.dump(2) << endl;

    return 0;
}
