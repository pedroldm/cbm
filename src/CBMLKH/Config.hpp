#ifndef CBMLKH_CONFIG_HPP
#define CBMLKH_CONFIG_HPP

#include <string>

using namespace std;

struct Config {
    string instancePath;
    int threads = 1;
    int maxIterations = 1000;
    int maxTime = 60;
    int lkhMaxTime = 5;
    float constructionBias = 1.0f;
};

#endif
