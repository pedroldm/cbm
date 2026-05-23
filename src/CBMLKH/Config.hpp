#ifndef CBMLKH_CONFIG_HPP
#define CBMLKH_CONFIG_HPP

#include <string>

using namespace std;

struct Config {
    string instancePath;
    int lkhMaxTime = 0;
    int threads = 1;
    float constructionBias = 1.0f;
};

#endif
