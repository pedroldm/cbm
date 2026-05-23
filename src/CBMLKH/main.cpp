#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "CBMLKH.hpp"
#include "PrintUtil.hpp"

using namespace std;

string stripPrefix(const string& arg, const string& prefix) {
    if (arg.rfind(prefix, 0) == 0) return arg.substr(prefix.size());
    return {};
}

template <typename T>
T parseValue(const string& val) {
    istringstream ss(val);
    T result;
    ss >> result;
    return result;
}

Config parseArgs(int argc, char* argv[]) {
    Config cfg;

    auto tryParseArg = [&](const string& arg, const string& flag, auto& field) -> bool {
        string val = stripPrefix(arg, "--" + flag + "=");
        if (val.empty()) return false;
        field = parseValue<decay_t<decltype(field)>>(val);
        return true;
    };

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        if (tryParseArg(arg, "instancePath", cfg.instancePath)) {
        } else if (tryParseArg(arg, "lkhMaxTime", cfg.lkhMaxTime)) {
        } else if (tryParseArg(arg, "constructionBias", cfg.constructionBias)) {
        } else {
            throw runtime_error("Unknown argument: " + arg);
        }
    }

    return cfg;
}

int main(int argc, char* argv[]) {
    const Config cfg = parseArgs(argc, argv);
    CBMLKH solver(cfg);
    Solution s = solver.greedyConstruction();
    vector<int> d = solver.count1BlocksPerColumn(s);
    cout << "1-blocks per column: " << d << "\n";
    return 0;
}