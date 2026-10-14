# Makefile for IceNine C++ Test Harness Programs

CXX = c++
CXXFLAGS = -std=c++11 -O2 -I../.. -I../../Src -I../../XDM++/libXDM -I/opt/homebrew/include
LDFLAGS = -L../../XDM++/libXDM -lXDM++ -L/opt/homebrew/lib

# Targets
TARGET1 = generate_test_data
SRC1 = generate_test_data.cpp

TARGET2 = test_neighbors
SRC2 = test_neighbors.cpp

TARGET3 = generate_omega_test_data
SRC3 = generate_omega_test_data.cpp ../../Src/InitFilesIO.cpp

# Ensure output directory exists
OUTPUT_DIR = ../cpp_outputs

.PHONY: all clean run test_neighbors_run test_omega_ranges

all: $(TARGET1) $(TARGET2) $(TARGET3)

# Build generate_test_data
$(TARGET1): $(SRC1)
	@echo "Compiling C++ test data generator..."
	@mkdir -p $(OUTPUT_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Build test_neighbors
$(TARGET2): $(SRC2)
	@echo "Compiling neighbor validation test..."
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Build generate_omega_test_data
$(TARGET3): $(SRC3)
	@echo "Compiling omega range test data generator..."
	@mkdir -p $(OUTPUT_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

run: $(TARGET1)
	@echo "Running test data generator..."
	@mkdir -p $(OUTPUT_DIR)
	./$(TARGET1)

test_neighbors_run: $(TARGET2)
	@echo "Running neighbor validation test..."
	./$(TARGET2) ../../DataFiles/Au1007_small.mic neighbors_cpp.txt

test_omega_ranges: $(TARGET3)
	@echo "Running omega range test data generator..."
	@mkdir -p $(OUTPUT_DIR)
	./$(TARGET3) $(OUTPUT_DIR)/omega_range_test_data.json

clean:
	rm -f $(TARGET1) $(TARGET2) $(TARGET3) *.o
	rm -f $(OUTPUT_DIR)/*.json neighbors_*.txt

.PHONY: help
help:
	@echo "IceNine C++ Test Harness"
	@echo "========================"
	@echo ""
	@echo "Targets:"
	@echo "  all                   - Compile all test programs"
	@echo "  generate_test_data    - Compile test data generator only"
	@echo "  test_neighbors        - Compile neighbor validation test only"
	@echo "  generate_omega_test_data - Compile omega range test generator only"
	@echo "  run                   - Run test data generator (generates JSON)"
	@echo "  test_neighbors_run    - Run neighbor validation on Au1007_small.mic"
	@echo "  test_omega_ranges     - Generate omega range test data JSON"
	@echo "  clean                 - Remove all generated files"
	@echo ""
	@echo "Outputs:"
	@echo "  - JSON files in ../cpp_outputs/"
	@echo "  - neighbors_cpp.txt (neighbor validation output)"
	@echo "  - omega_range_test_data.json (omega test validation)"
