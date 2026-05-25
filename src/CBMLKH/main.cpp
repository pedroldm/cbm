#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CBMLKH.hpp"
#include "LKHWrapper.hpp"
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
    Config cfg = parseArgs(argc, argv);

    try {
        CBMLKH solver(cfg);
        Solution base = solver.greedyConstruction();
        LKHWrapper LKHWrapper(solver.c, solver.l, solver.tspMatrix);
        base.sol = LKHWrapper.run(base.sol, "a1", 10);
        solver.completeEval(base);
        cout << "Final solution cost: " << base.cost << "\n";
    } catch (const exception& ex) {
        cerr << "CBMLKH test failed: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}