#ifndef LKH_CACHE_HPP
#define LKH_CACHE_HPP

#include <algorithm>
#include <atomic>  // Required for thread-safe counters
#include <iostream>
#include <map>
#include <mutex>         // Required for std::unique_lock / std::shared_lock
#include <shared_mutex>  // Required for std::shared_mutex
#include <vector>

class LKHCache {
   private:
    std::map<std::vector<int>, std::vector<int>> cacheMap;

    // Mutex to protect cacheMap. Marked mutable so const functions can lock it.
    mutable std::shared_mutex cacheMutex;

    // Atomic counters to prevent data races during concurrent reads
    std::atomic<long long> hits{0};
    std::atomic<long long> misses{0};

   public:
    LKHCache() = default;

    // Explicitly disallow copying because mutexes cannot be copied
    LKHCache(const LKHCache&) = delete;
    LKHCache& operator=(const LKHCache&) = delete;

    /**
     * Thread-safe read operation. Multiple threads can call this concurrently.
     */
    bool get(const std::vector<int>& elements, std::vector<int>& outSolution) {
        std::vector<int> canonicalKey = elements;
        std::sort(canonicalKey.begin(), canonicalKey.end());

        // 'shared_lock' allows shared (read-only) access
        std::shared_lock<std::shared_mutex> lock(cacheMutex);

        auto it = cacheMap.find(canonicalKey);
        if (it != cacheMap.end()) {
            hits++;  // Safe to increment across threads because it's atomic
            outSolution = it->second;
            return true;
        }

        misses++;
        return false;
    }

    /**
     * Thread-safe write operation. Blocks all other readers and writers.
     */
    void put(const std::vector<int>& inputElements, const std::vector<int>& optimizedSolution) {
        std::vector<int> canonicalKey = inputElements;
        std::sort(canonicalKey.begin(), canonicalKey.end());

        // 'unique_lock' demands exclusive (write) access
        std::unique_lock<std::shared_mutex> lock(cacheMutex);
        cacheMap[canonicalKey] = optimizedSolution;
    }

    void clear() {
        std::unique_lock<std::shared_mutex> lock(cacheMutex);
        cacheMap.clear();
        hits = 0;
        misses = 0;
    }

    long long getHits() const { return hits.load(); }
    long long getMisses() const { return misses.load(); }

    void printStats() const {
        std::shared_lock<std::shared_mutex> lock(cacheMutex);
        long long localHits = hits.load();
        long long localMisses = misses.load();
        long long total = localHits + localMisses;
        double rate = (total > 0) ? (100.0 * localHits / total) : 0.0;
        std::cout << "[LKH Cache] Total Requests: " << total << " | Hits: " << localHits << " | Misses: " << localMisses << " | Hit Rate: " << rate
                  << "%\n";
    }
};

#endif  // LKH_CACHE_HPP