#include "ArgsUtil.hpp"

#include <fstream>
#include <functional>
#include <stdexcept>
#include <unordered_map>

using namespace std;

Config ArgsUtil::parseConfigFile(const string& path) {
    Config cfg;

    ifstream file(path);

    if (!file.is_open()) {
        throw runtime_error("Cannot open config file: " + path);
    }

    unordered_map<string, function<void(const string&)>> parsers;

    // General
    parsers["instancePath"] = [&](const string& s) { cfg.instancePath = s; };
    parsers["threads"] = [&](const string& s) { cfg.threads = stoi(s); };
    parsers["maxIterations"] = [&](const string& s) { cfg.maxIterations = stoi(s); };
    parsers["maxTime"] = [&](const string& s) { cfg.maxTime = stoi(s); };
    parsers["lkhMaxTime"] = [&](const string& s) { cfg.lkhMaxTime = stoi(s); };

    // Base values
    parsers["constructionBias"] = [&](const string& s) { cfg.constructionBias = stod(s); };
    parsers["neighborBias"] = [&](const string& s) { cfg.neighborBias = stod(s); };
    parsers["maxSegmentSize"] = [&](const string& s) { cfg.maxSegmentSizeFraction = stod(s); };
    parsers["minSegmentScore"] = [&](const string& s) { cfg.minSegmentScore = stod(s); };

    // Adaptation control
    parsers["adaptationInterval"] = [&](const string& s) { cfg.adaptationInterval = stoi(s); };

    // Bounds
    parsers["maxSegmentSizeUpperBound"] = [&](const string& s) { cfg.maxSegmentSizeUpperBoundFraction = stod(s); };

    parsers["minSegmentScoreLowerBound"] = [&](const string& s) { cfg.minSegmentScoreLowerBound = stod(s); };

    parsers["minNeighborBias"] = [&](const string& s) { cfg.minNeighborBias = stod(s); };

    // Multipliers
    parsers["segmentSizeGrowthFactor"] = [&](const string& s) { cfg.segmentSizeGrowthFactor = stod(s); };
    parsers["segmentScoreDecayFactor"] = [&](const string& s) { cfg.segmentScoreDecayFactor = stod(s); };
    parsers["neighborBiasDecayFactor"] = [&](const string& s) { cfg.neighborBiasDecayFactor = stod(s); };

    string line;
    int lineNumber = 0;

    while (getline(file, line)) {
        ++lineNumber;

        // Ignore empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        auto pos = line.find('=');

        if (pos == string::npos) {
            throw runtime_error("Invalid line " + to_string(lineNumber) + " in config file: " + line);
        }

        string key = line.substr(0, pos);
        string value = line.substr(pos + 1);

        auto it = parsers.find(key);

        if (it != parsers.end()) {
            it->second(value);
        } else if (key == "blockMovement") {
            if (value == "RANDOM") {
                cfg.blockMovement = RANDOM;
            } else if (value == "PEAK") {
                cfg.blockMovement = PEAK;
            } else if (value == "INTERVAL") {
                cfg.blockMovement = INTERVAL;
            } else if (value == "MERGE") {
                cfg.blockMovement = MERGE;
            } else {
                throw runtime_error("Invalid blockMovement value at line " + to_string(lineNumber) + ": " + value +
                                    ". Expected RANDOM, PEAK, or INTERVAL.");
            }
        } else {
            throw runtime_error("Unknown configuration option at line " + to_string(lineNumber) + ": " + key);
        }
    }

    return cfg;
}