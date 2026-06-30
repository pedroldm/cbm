# Compiler
CXX = g++

# Compiler flags for development build
CXXFLAGS = -fopenmp -Wall -Wextra -std=c++17 -O2

# Compiler flags for debug build
CXXFLAGS_DEBUG = -fopenmp -Wall -Wextra -std=c++17 -g

# Compiler flags for production build
CXXFLAGS_PRD = -fopenmp -std=c++17 -O3 -march=native

# Source files
SRCS = ArgsUtil.cpp Validator.cpp Trajectory.cpp LKHWrapper.cpp CBMLKH.cpp main.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)
OBJS_DEBUG = $(SRCS:.cpp=.debug.o)
OBJS_PRD = $(SRCS:.cpp=.prd.o)

# Header files
HEADERS = ArgsUtil.hpp Config.hpp ColumnStore.hpp LKHCache.hpp Validator.hpp Trajectory.hpp LKHWrapper.hpp CBMLKH.hpp Solution.hpp PrintUtil.hpp Metrics.hpp ../IO/json.hpp

# Target executable
TARGET = main
TARGET_DEBUG = main_debug
TARGET_PRD = main_prd

# Default target (development build)
dev: $(TARGET)

# Debug target
debug: $(TARGET_DEBUG)

# Production target
prd: $(TARGET_PRD)

# Link the object files to create the development executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Link the object files to create the debug executable
$(TARGET_DEBUG): $(OBJS_DEBUG)
	$(CXX) $(CXXFLAGS_DEBUG) -o $(TARGET_DEBUG) $(OBJS_DEBUG)

# Link the object files to create the production executable
$(TARGET_PRD): $(OBJS_PRD)
	$(CXX) $(CXXFLAGS_PRD) -o $(TARGET_PRD) $(OBJS_PRD)

# Compile source files into development object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile source files into debug object files
%.debug.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS_DEBUG) -c $< -o $@

# Compile source files into production object files
%.prd.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS_PRD) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(OBJS_DEBUG) $(OBJS_PRD) \
	$(TARGET) $(TARGET_DEBUG) $(TARGET_PRD)

# Phony targets
.PHONY: dev debug prd clean