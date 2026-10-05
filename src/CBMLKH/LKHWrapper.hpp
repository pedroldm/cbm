#ifndef LKHWrapper_HPP
#define LKHWrapper_HPP

#include <atomic>
#include <cstdint>
#include <string>
#include <vector>

#include "ColumnStore.hpp"

struct LKHResult {
    std::vector<int> tour;      // the input columns, reordered
    bool timeLimitHit = false;  // LKH stopped on TIME_LIMIT (its result is then load-dependent)
};

// Solves the open-path TSP over a set of columns with an external LKH process.
// Thread-safe: every call uses its own files inside a private per-process
// scratch directory, created by configure() and removed by cleanup().
class LKHWrapper {
    inline static std::string lkhPath;
    inline static std::string workDir;
    inline static std::atomic<long> nextCallId{0};

    const ColumnStore& columns;

   public:
    // `tmpParent` may be shared between processes: the scratch directory is a
    // fresh mkdtemp() child of it, so no process ever deletes another's files.
    static void configure(const std::string& executablePath, const std::string& tmpParent);
    static void cleanup();
    static const std::string& scratchDir() { return workDir; }

    explicit LKHWrapper(const ColumnStore& columns);
    LKHResult run(const std::vector<int>& slice, int maxTime, uint32_t seed) const;

   private:
    void writeTSP(const std::vector<int>& slice, const std::string& tspFile) const;
    static std::vector<int> readTour(const std::string& tourFile);
};

#endif
