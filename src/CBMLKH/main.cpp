#include <exception>
#include <iostream>

#include "ArgsUtil.hpp"
#include "CBMLKH.hpp"
#include "LKHWrapper.hpp"
#include "PrintUtil.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    Config cfg = ArgsUtil::parseArgs(argc, argv);

    try {
        CBMLKH solver(cfg);
        Solution base = solver.greedyConstruction();
        solver.countBlocksPerColumn(base);
        solver.findDenseSegments(base, 2, 10, 5.0);
    } catch (const exception& ex) {
        cerr << "CBMLKH test failed: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}