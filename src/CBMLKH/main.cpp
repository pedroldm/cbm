#include <algorithm>
#include <chrono>
#include <exception>
#include <filesystem>
#include <fstream>
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
        {"lkhTimeMs", m.lkhTimeMs},
        {"algorithmTimeMs", m.algorithmTimeMs()},
        {"acceptedMoves", m.acceptedMoves},
        {"rejectedMoves", m.rejectedMoves},
        {"lkhCalls", m.lkhCalls},
        {"lkhCacheMisses", m.lkhCacheMisses},
        {"lkhTimeLimitHits", m.lkhTimeLimitHits},
        {"stopReason", m.stopReason},
        {"diversifications", m.diversifications},
        // Histogram stop criterion. The flatness fields are null when the
        // criterion is disabled (they stay at infinity, which JSON cannot hold).
        {"histogramStop",
         {{"triggered", m.histogramStop},
          {"flattenings", m.flattenings},
          {"initialFlatness", m.initialFlatness},
          {"bestFlatness", m.bestFlatness},
          {"iterationsWithoutFlattening", m.iterationsWithoutFlattening}}},
        {"blockSize", {{"average", m.averageBlockSize()}, {"min", m.minBlockSize()}, {"max", m.blockSizeMax}, {"count", m.blockCount}}},
        {"operators", {{"PEAK", operatorToJson(m.peak)}, {"INTERVAL", operatorToJson(m.interval)}, {"MERGE", operatorToJson(m.merge)}}},
        {"neighborBiasHistory", m.neighborBiasHistory},
    };
}

// Write-then-rename, so a reader never sees a truncated report.
static void writeAtomically(const string& path, const string& content) {
    const string tmp = path + ".tmp";
    {
        ofstream out(tmp);
        out << content;
        out.flush();
        if (!out) throw runtime_error("Cannot write " + tmp);
    }
    filesystem::rename(tmp, path);
}

static int run(int argc, char* argv[]) {
    if (argc != 2) throw runtime_error("Usage: ./cbmlkh <config_file>");

    Config cfg = ArgsUtil::parseConfigFile(argv[1]);

    LKHWrapper::configure(cfg.lkhPath, cfg.lkhTmpDir);
    auto cache = make_shared<LKHCache>();
    CBMLKH cbmlkh(cfg, cache);

    auto runStart = chrono::steady_clock::now();
    vector<Trajectory> trajectories = cbmlkh.run();
    auto runtimeMs = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - runStart).count();

    // iRace mode: emit nothing but the objective value, so the tuner's
    // target-runner can read it straight off stdout without parsing the report.
    if (cfg.iRace) {
        int bestCost = -1;
        for (const Trajectory& traj : trajectories) {
            if (bestCost < 0 || traj.bestSolution.cost < bestCost) bestCost = traj.bestSolution.cost;
        }
        cout << bestCost << endl;
        return 0;
    }

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
            history.push_back({{"iteration", e.iteration},
                               {"cost", e.cost},
                               {"elapsedMs", e.elapsedMs},
                               {"move", toString(e.movement)},
                               {"histogram", e.histogram}});
        }

        // 1-indexed column ids, as in the instance file.
        vector<int> permutation(traj.bestSolution.sol.size());
        for (size_t k = 0; k < permutation.size(); k++) permutation[k] = traj.bestSolution.sol[k] + 1;

        json tj = metricsToJson(m);
        tj["index"] = traj.index;
        tj["seed"] = traj.seed;
        tj["timeToBestMs"] = traj.history.empty() ? 0 : traj.history.back().elapsedMs;
        tj["history"] = history;
        tj["bestPermutation"] = permutation;
        trajectoriesJson.push_back(tj);
    }

    long long hits = cache->getHits();
    long long misses = cache->getMisses();
    long long requests = hits + misses;
    double hitRate = requests > 0 ? (100.0 * hits / requests) : 0.0;

    // Report the config CBMLKH actually ran with: its copy is the one whose
    // fractional segment sizes were resolved against the instance's column count.
    const Config& resolved = cbmlkh.cfg;

    json output = {
        {"instance", {{"name", cbmlkh.instanceName}, {"rows", cbmlkh.rows}, {"cols", cbmlkh.cols}}},
        {"config",
         {{"seed", resolved.seed},
          {"seedBase", cbmlkh.seedBase},
          {"reproducible", resolved.seed != 0},
          {"trajectoryOffset", resolved.trajectoryOffset},
          {"threads", resolved.threads},
          {"blockMovement", toString(resolved.blockMovement)},
          {"maxIterations", resolved.maxIterations},
          {"maxTime", resolved.maxTime},
          {"lkhMaxTime", resolved.lkhMaxTime},
          {"constructionBias", resolved.constructionBias},
          {"neighborBias", resolved.neighborBias},
          {"minNeighborBias", resolved.minNeighborBias},
          {"neighborBiasDecayFactor", resolved.neighborBiasDecayFactor},
          {"minSegmentSize", resolved.minSegmentSize},
          {"maxSegmentSizeFraction", resolved.maxSegmentSizeFraction},
          {"maxSegmentSizeUpperBoundFraction", resolved.maxSegmentSizeUpperBoundFraction},
          {"maxSegmentSize", resolved.maxSegmentSize},
          {"maxSegmentSizeUpperBound", resolved.maxSegmentSizeUpperBound},
          {"maxMergeSegmentSize", std::min(cbmlkh.cols, 2 * resolved.maxSegmentSizeUpperBound)},
          {"minSegmentScore", resolved.minSegmentScore},
          {"minSegmentScoreLowerBound", resolved.minSegmentScoreLowerBound},
          {"segmentSizeGrowthFactor", resolved.segmentSizeGrowthFactor},
          {"segmentScoreDecayFactor", resolved.segmentScoreDecayFactor},
          {"adaptationInterval", resolved.adaptationInterval},
          {"histogramStopInterval", resolved.histogramStopInterval},
          {"histogramStopMeasure", toString(resolved.histogramStopMeasure)}}},
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

    if (resolved.outputPath.empty()) {
        cout << output.dump(2) << endl;
    } else {
        writeAtomically(resolved.outputPath, output.dump(2) + "\n");
    }
    return 0;
}

int main(int argc, char* argv[]) {
    // Config errors are the common failure mode; surface the message instead of
    // letting the exception escape and abort with a bare "terminate called".
    int status;
    try {
        status = run(argc, argv);
    } catch (const exception& e) {
        cerr << e.what() << endl;
        status = 1;
    }
    LKHWrapper::cleanup();
    return status;
}
