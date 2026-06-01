#ifndef ARGSUTIL_HPP
#define ARGSUTIL_HPP

#include <sstream>
#include <stdexcept>
#include <string>

#include "CBMLKH.hpp"

class ArgsUtil {
   private:
    static std::string stripPrefix(const std::string& arg, const std::string& prefix);

    template <typename T>
    static T parseValue(const std::string& val) {
        std::istringstream ss(val);
        T result;
        ss >> result;
        return result;
    }

   public:
    static Config parseArgs(int argc, char* argv[]);
};

#endif