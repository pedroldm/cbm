#include "Validator.hpp"

#include <fstream>
#include <stdexcept>
#include <unordered_set>

Validator::Validator(const std::string& instancePath) {
    std::ifstream f(instancePath);
    if (!f.is_open()) throw std::runtime_error("Error opening file: " + instancePath);

    f >> lines >> columns;
    binaryMatrix.assign(lines, std::vector<int>(columns, 0));

    for (int i = 0; i < lines; i++) {
        int count;
        f >> count;
        for (int j = 0; j < count; j++) {
            int e;
            f >> e;
            binaryMatrix[i][e - 1] = 1;
        }
    }
}

bool Validator::validate(const std::vector<int>& perm, int cost) {
    bool v1 = containsAllColumns(perm);
    int v2 = countOneBlocks(perm);
    if (!v1 || cost != v2)
        throw std::runtime_error("Invalid solution: containsAll=" + std::to_string(v1) + " countBlocks=" + std::to_string(v2) +
                                 " != cost=" + std::to_string(cost));
    return true;
}

int Validator::countOneBlocks(const std::vector<int>& perm) {
    int blocks = 0;
    for (int row = 0; row < lines; row++) {
        bool inBlock = false;
        for (int col : perm) {
            if (binaryMatrix[row][col] == 1) {
                if (!inBlock) {
                    blocks++;
                    inBlock = true;
                }
            } else {
                inBlock = false;
            }
        }
    }
    return blocks;
}

bool Validator::containsAllColumns(const std::vector<int>& perm) {
    std::unordered_set<int> permSet(perm.begin(), perm.end());
    for (int i = 0; i < columns; i++)
        if (permSet.find(i) == permSet.end()) return false;
    return true;
}