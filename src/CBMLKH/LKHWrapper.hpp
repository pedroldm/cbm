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
    // Override the compiled-in defaults from the config file. Must be called
    // before clearTmpDir() and before any LKHWrapper is constructed.
    //
    // A per-process tmpDir is what makes concurrent runs safe: clearTmpDir()
    // empties the directory at startup, so several processes sharing one tmpDir
    // (e.g. irace with --parallel > 1) would delete each other's in-flight .tsp
    // and tour files.
    static void configure(const std::string& executablePath, const std::string& temporaryDir);

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
