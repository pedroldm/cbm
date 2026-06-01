#ifndef VALIDATOR_HPP
#define VALIDATOR_HPP

#include <string>
#include <vector>

class Validator {
   public:
    int lines, columns;
    std::vector<std::vector<int>> binaryMatrix;

    explicit Validator(const std::string& instancePath);
    bool validate(const std::vector<int>& perm, int cost);

   private:
    int countOneBlocks(const std::vector<int>& perm);
    bool containsAllColumns(const std::vector<int>& perm);
};

#endif