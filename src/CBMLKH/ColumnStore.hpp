#ifndef COLUMN_STORE_HPP
#define COLUMN_STORE_HPP

#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

/*
 * ColumnStore: bitset representation of a CBM instance.
 *
 * Each of the `cols` columns is stored as a packed bit vector over the `rows`
 * rows (column-major: bit `row` of a column is set iff A[row][col] == 1).
 * This replaces the previous dense cols x cols distance/transition matrices,
 * whose O(cols^2) memory and O(cols^2 * rows) precompute made very large
 * instances (columns up to ~45k => tens of GB) infeasible.
 *
 * All pairwise quantities are computed on demand from the bitsets:
 *   - hamming(a, b)     = #rows where columns a and b differ
 *   - zerosToOnes(a, b) = #rows where a == 0 and b == 1 (new 1-blocks at seam a->b)
 *   - onesToZeros(a, b) = #rows where a == 1 and b == 0
 *   - onesCount(col)    = #rows where column col == 1
 *
 * Memory is cols * ceil(rows / 64) * 8 bytes, with no O(cols^2) precompute.
 */
class ColumnStore {
   public:
    ColumnStore() = default;

    explicit ColumnStore(const std::string& instancePath) { load(instancePath); }

    void load(const std::string& instancePath) {
        std::ifstream input(instancePath);
        if (!input.is_open()) throw std::runtime_error("Error opening file: " + instancePath);

        input >> numRows >> numCols;
        if (!input) throw std::runtime_error("Malformed instance header: " + instancePath);

        wordsPerColumn = (numRows + 63) / 64;
        data.assign(static_cast<size_t>(numCols) * wordsPerColumn, 0ULL);

        for (int row = 0; row < numRows; ++row) {
            int onesInRow;
            input >> onesInRow;
            if (!input) throw std::runtime_error("Malformed instance body: " + instancePath);
            for (int k = 0; k < onesInRow; ++k) {
                int column;
                input >> column;
                if (column < 1 || column > numCols) throw std::runtime_error("Column index out of range in: " + instancePath);
                setOne(row, column - 1);
            }
        }
    }

    int rows() const { return numRows; }  // number of rows (l)
    int cols() const { return numCols; }  // number of columns / TSP cities (c)

    // #rows where column `column` has a 1.
    int onesCount(int column) const {
        const uint64_t* bits = columnBits(column);
        int total = 0;
        for (int word = 0; word < wordsPerColumn; ++word) total += __builtin_popcountll(bits[word]);
        return total;
    }

    // #rows where columns a and b differ (symmetric Hamming distance).
    int hamming(int a, int b) const {
        const uint64_t* bitsA = columnBits(a);
        const uint64_t* bitsB = columnBits(b);
        int total = 0;
        for (int word = 0; word < wordsPerColumn; ++word) total += __builtin_popcountll(bitsA[word] ^ bitsB[word]);
        return total;
    }

    // #rows where a == 0 and b == 1 (number of new 1-blocks opened at seam a -> b).
    int zerosToOnes(int a, int b) const {
        const uint64_t* bitsA = columnBits(a);
        const uint64_t* bitsB = columnBits(b);
        int total = 0;
        for (int word = 0; word < wordsPerColumn; ++word) total += __builtin_popcountll(bitsB[word] & ~bitsA[word]);
        return total;
    }

    // #rows where a == 1 and b == 0.
    int onesToZeros(int a, int b) const {
        const uint64_t* bitsA = columnBits(a);
        const uint64_t* bitsB = columnBits(b);
        int total = 0;
        for (int word = 0; word < wordsPerColumn; ++word) total += __builtin_popcountll(bitsA[word] & ~bitsB[word]);
        return total;
    }

   private:
    int numRows = 0;
    int numCols = 0;
    int wordsPerColumn = 0;
    std::vector<uint64_t> data;

    const uint64_t* columnBits(int column) const { return data.data() + static_cast<size_t>(column) * wordsPerColumn; }
    uint64_t* columnBits(int column) { return data.data() + static_cast<size_t>(column) * wordsPerColumn; }

    void setOne(int row, int column) { columnBits(column)[row >> 6] |= (1ULL << (row & 63)); }
};

#endif  // COLUMN_STORE_HPP
