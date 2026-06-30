#ifndef CBMLKH_HPP
#define CBMLKH_HPP

#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "ColumnStore.hpp"
#include "Config.hpp"
#include "LKHCache.hpp"
#include "LKHWrapper.hpp"
#include "Solution.hpp"
#include "Trajectory.hpp"
#include "Validator.hpp"

struct CandidateRegion {
    int start;
    int end;
    double score;
};

class CBMLKH {
   public:
    Config cfg;
    Validator validator;
    ColumnStore columnStore;

    int rows;  // number of matrix rows (l)
    int cols;  // number of columns / TSP cities (c)
    std::string instanceName;

    CBMLKH(const Config& cfg, std::shared_ptr<LKHCache> cache);

    std::vector<Trajectory> run();
    Trajectory LKHILS(Solution& initial);
    Solution ILSNeighbor(const Solution& s, const AdaptiveParameters& adaptive, Metrics& metrics);
    CandidateRegion choosePeakRegion(Solution& s, const AdaptiveParameters& adaptive);
    CandidateRegion chooseIntervalRegion(Solution& s, const AdaptiveParameters& adaptive);
    int completeEval(Solution& s);
    Solution greedyConstruction();
    int nextInsertion(int current, std::unordered_set<int>& remaining);
    int countBlocksPerColumn(Solution& s, int start = 0, int end = -1);
    std::vector<CandidateRegion> findPeakColumns(Solution& s, const AdaptiveParameters& adaptive);
    std::vector<CandidateRegion> findDenseSegments(Solution& s, const AdaptiveParameters& adaptive);
    CandidateRegion chooseMergeRegion(Solution& s, const AdaptiveParameters& adaptive);
    void applyLKH(Solution& s, const CandidateRegion& cr, Metrics& metrics);
    void reinsertBlock(Solution& s, const CandidateRegion& cr, const std::vector<int>& block);

   private:
    // Rank-based roulette: returns an index in [0, count) where rank r is
    // weighted by 1 / (r + 1)^neighborBias (lower-ranked candidates favored).
    std::size_t sampleRankIndex(std::size_t count, double neighborBias);

    LKHWrapper lkhWrapper;
    std::shared_ptr<LKHCache> lkhCache;
};

#endif  // CBMLKH_HPP
