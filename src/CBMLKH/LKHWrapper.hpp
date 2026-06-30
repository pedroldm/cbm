#ifndef LKHWrapper_HPP
#define LKHWrapper_HPP

#include <string>
#include <vector>

#include "ColumnStore.hpp"

class LKHWrapper {
    inline static std::string lkhPath = "/home/pedroldm/MSc/cbm/src/LKH3/LKH";
    inline static std::string tmpDir = "/tmp/LKH/";

    const ColumnStore& columns;

   public:
    explicit LKHWrapper(const ColumnStore& columns);
    std::vector<int> run(const std::vector<int>& slice, std::string instanceName, int maxTime);
    void writeTSP(const std::vector<int>& slice, std::string instanceName, std::string tspFile, long execId);
    static void clearTmpDir();
    void writePar(std::string parFile, std::string tspFile, std::string resultTourFile, int maxTime);
    void runLKH(std::string parFile);
    long getExecutionId();
    std::vector<int> getResultTour(std::string resultTourFile);
};

#endif
