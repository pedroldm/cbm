#include "ArgsUtil.hpp"

using namespace std;

string ArgsUtil::stripPrefix(const string& arg, const string& prefix) {
    if (arg.rfind(prefix, 0) == 0) {
        return arg.substr(prefix.size());
    }

    return {};
}

Config ArgsUtil::parseArgs(int argc, char* argv[]) {
    Config cfg;

    auto tryParseArg = [&](const string& arg, const string& flag, auto& field) -> bool {
        string val = stripPrefix(arg, "--" + flag + "=");

        if (val.empty()) {
            return false;
        }

        field = parseValue<decay_t<decltype(field)>>(val);
        return true;
    };

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];

        if (tryParseArg(arg, "instancePath", cfg.instancePath)) {
        } else if (tryParseArg(arg, "lkhMaxTime", cfg.lkhMaxTime)) {
        } else if (tryParseArg(arg, "constructionBias", cfg.constructionBias)) {
        } else if (tryParseArg(arg, "maxIterations", cfg.maxIterations)) {
        } else if (tryParseArg(arg, "maxTime", cfg.maxTime)) {
        } else if (tryParseArg(arg, "threads", cfg.threads)) {
        } else {
            throw runtime_error("Unknown argument: " + arg);
        }
    }

    return cfg;
}