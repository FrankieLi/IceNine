#!/bin/bash
#
# Validate Python Detector implementation against C++
#
# This script:
# 1. Compiles the C++ validation test
# 2. Runs C++ and Python tests
# 3. Compares outputs
#

set -e

echo "=== Detector Validation: Python vs C++ ==="
echo

# Check if C++ files exist
if [ ! -f "../../Src/Detector.cpp" ]; then
    echo "ERROR: C++ source files not found"
    echo "Expected to find: ../../Src/Detector.cpp"
    echo "Current directory: $(pwd)"
    echo "Please run this script from icenine_py/tests/"
    exit 1
fi

# Compile C++ validation test
echo "1. Compiling C++ validation test..."
BOOST_PATH="/opt/homebrew/include"
g++ -std=c++11 \
    -I$BOOST_PATH \
    -I../../XDM++/libXDM \
    -I../../Src \
    validate_cpp_detector.cpp \
    ../../XDM++/libXDM/3dMath.cpp \
    ../../Src/Detector.cpp \
    -o validate_cpp_detector

if [ $? -ne 0 ]; then
    echo "ERROR: C++ compilation failed"
    exit 1
fi
echo "   ✓ C++ compiled successfully"
echo

# Run C++ test
echo "2. Running C++ validation..."
./validate_cpp_detector > cpp_detector_output.txt
echo "   ✓ C++ output saved to cpp_detector_output.txt"
echo

# Run Python test
echo "3. Running Python validation..."
python3 validate_py_detector.py > py_detector_output.txt
echo "   ✓ Python output saved to py_detector_output.txt"
echo

# Compare outputs
echo "4. Comparing outputs..."
echo

if diff -u cpp_detector_output.txt py_detector_output.txt > validation_diff.txt; then
    echo "   ✅ SUCCESS: Python and C++ outputs match exactly!"
    echo
    rm validation_diff.txt
else
    echo "   ⚠️  DIFFERENCES FOUND"
    echo
    echo "Differences:"
    cat validation_diff.txt
    echo
    echo "Full diff saved to: validation_diff.txt"
    echo
    echo "To investigate:"
    echo "  - Check cpp_detector_output.txt"
    echo "  - Check py_detector_output.txt"
    echo "  - Review validation_diff.txt"
fi

echo
echo "=== Validation Complete ==="
