# Compiler
CXX = g++

# Compiler flags for development build
CXXFLAGS = -fopenmp -Wall -Wextra -std=c++17 -O2

# Compiler flags for debug build
CXXFLAGS_DEBUG = -fopenmp -Wall -Wextra -std=c++17 -g

# Compiler flags for production build
CXXFLAGS_PRD = -fopenmp -std=c++17 -O3 -march=native

# Compiler flags for delta-eval verification build: recomputes each move's cost
# with a full completeEval and asserts it matches the O(blockLen) delta result.
CXXFLAGS_VERIFY = -fopenmp -Wall -Wextra -std=c++17 -O2 -DDELTA_EVAL_VERIFY

# Source files
SRCS = ArgsUtil.cpp Validator.cpp Trajectory.cpp LKHWrapper.cpp CBMLKH.cpp main.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)
OBJS_DEBUG = $(SRCS:.cpp=.debug.o)
OBJS_PRD = $(SRCS:.cpp=.prd.o)
OBJS_VERIFY = $(SRCS:.cpp=.verify.o)

# Header files
HEADERS = ../common/cbm_lkh_params.h ../common/cbm_seed.h ../common/cbm_spawn.h ArgsUtil.hpp Config.hpp ColumnStore.hpp LKHCache.hpp Validator.hpp Trajectory.hpp LKHWrapper.hpp CBMLKH.hpp Solution.hpp PrintUtil.hpp Metrics.hpp Flatness.hpp ../IO/json.hpp

# Target executable
TARGET = main
TARGET_DEBUG = main_debug
TARGET_PRD = main_prd
TARGET_VERIFY = main_verify

# Default target (development build)
dev: $(TARGET)

# Standalone LKH baseline (whole instance as one TSP)
lkh_standalone: lkh_standalone.cpp ../common/cbm_lkh_params.h ColumnStore.hpp ../common/cbm_spawn.h ../IO/json.hpp
	$(CXX) $(CXXFLAGS_PRD) -o lkh_standalone lkh_standalone.cpp

# Debug target
debug: $(TARGET_DEBUG)

# Production target
prd: $(TARGET_PRD)

# Delta-eval verification target
verify: $(TARGET_VERIFY)

# Link the object files to create the development executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Link the object files to create the debug executable
$(TARGET_DEBUG): $(OBJS_DEBUG)
	$(CXX) $(CXXFLAGS_DEBUG) -o $(TARGET_DEBUG) $(OBJS_DEBUG)

# Link the object files to create the production executable
$(TARGET_PRD): $(OBJS_PRD)
	$(CXX) $(CXXFLAGS_PRD) -o $(TARGET_PRD) $(OBJS_PRD)

# Link the object files to create the verification executable
$(TARGET_VERIFY): $(OBJS_VERIFY)
	$(CXX) $(CXXFLAGS_VERIFY) -o $(TARGET_VERIFY) $(OBJS_VERIFY)

# Compile source files into development object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile source files into debug object files
%.debug.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS_DEBUG) -c $< -o $@

# Compile source files into production object files
%.prd.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS_PRD) -c $< -o $@

# Compile source files into verification object files
%.verify.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS_VERIFY) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(OBJS_DEBUG) $(OBJS_PRD) $(OBJS_VERIFY) \
	$(TARGET) $(TARGET_DEBUG) $(TARGET_PRD) $(TARGET_VERIFY) lkh_standalone

# Phony targets
.PHONY: dev debug prd verify clean lkh_standalone