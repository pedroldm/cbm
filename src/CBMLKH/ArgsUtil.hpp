#ifndef ARGS_UTIL_HPP
#define ARGS_UTIL_HPP

#include <string>

#include "Config.hpp"

class ArgsUtil {
   public:
    static Config parseConfigFile(const std::string& path);
};

#endif