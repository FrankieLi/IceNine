# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

IceNine is a forward model reconstruction tool for synchrotron X-ray diffraction data, specifically designed for three-dimensional orientation imaging of polycrystalline materials. It reconstructs crystal grain orientations from experimental detector images using adaptive search algorithms.

**Citation**: S. F. Li, and R. M. Suter, "Adaptive reconstruction method for three-dimensional orientation imaging", Journal of Applied Crystallography, 2013.

## Build System

### Requirements
- C++ compiler (gcc >= 4.2, recommended >= 4.5, or Clang/AppleClang)
- MPI implementation (Open-MPI or MPICH)
- CMake 3.15+
- Boost 1.43+ (tested with 1.89.0)
- Eigen 3+ (tested with 5.0.0)

### Modern Build (macOS with Homebrew)

The build system has been modernized for CMake 3.15+ with automatic dependency detection.

1. **Install dependencies** (on macOS):
   ```bash
   brew install cmake boost eigen open-mpi
   ```

2. **Configure with CMake**:
   ```bash
   export PATH="/opt/homebrew/bin:$PATH"
   cmake -DCMAKE_BUILD_TYPE=Release .
   ```

3. **Build**:
   ```bash
   make -j8
   ```

4. **CRITICAL**: Ensure CMAKE_BUILD_TYPE is set to "Release" - debug builds run ~100x slower due to debug checks.

### Build Modernization Changes (November 2024)

The following changes were made to modernize the build system and ensure compatibility with modern compilers and C++11:

#### CMake Modernization
- Upgraded from CMake 2.4 to CMake 3.15
- Migrated to modern target-based CMake approach
- Added automatic dependency detection using `find_package()` for MPI, Boost, and Eigen3
- Removed hardcoded paths in favor of package discovery
- Added explicit C++11 standard requirement

**Files modified**:
- [CMakeLists.txt](CMakeLists.txt) - Main build configuration
- [XDM++/libXDM/CMakeLists.txt](XDM++/libXDM/CMakeLists.txt) - Library build configuration

#### C++11 Migration Fixes

**1. Boost Lambda to C++11 Lambda Migration** ([Src/DiscreteSearch.h:285,439](Src/DiscreteSearch.h))
- Replaced Boost Lambda `bind()` expressions with C++11 lambdas
- Old: `bind(&CRecpVector::fMag, _1) > fQMax`
- New: `[fQMax](const CRecpVector& v) { return v.fMag > fQMax; }`

**2. make_pair Template Argument Fixes**
- Removed explicit template arguments from `make_pair()` calls (C++11 has automatic type deduction)
- [XDM++/libXDM/Sampling.cpp:284](XDM++/libXDM/Sampling.cpp)
- [Src/LocalOptimizationAdaptor.h:78](Src/LocalOptimizationAdaptor.h)

**3. Reference Lifetime Fix** ([XDM++/libXDM/MicMesh.h:209](XDM++/libXDM/MicMesh.h))
- Changed `const ShapePtr &` to `ShapePtr` to avoid temporary object binding

#### Dependency Management

**Eigen3 Path Resolution**:
The build system now explicitly excludes old 3rdParty installations and uses only Homebrew-managed dependencies:
```cmake
# Force CMake to use only Homebrew Eigen3 and ignore 3rdParty directory
list(APPEND CMAKE_IGNORE_PATH "/Users/sfli/Research/3rdParty")
set(CMAKE_PREFIX_PATH "/opt/homebrew" ${CMAKE_PREFIX_PATH})
find_package(Eigen3 REQUIRED NO_MODULE PATHS /opt/homebrew/share/eigen3/cmake NO_DEFAULT_PATH)
```

**Compiler Flags**:
- Added `-Wno-deprecated-declarations` to suppress Boost deprecation warnings
- Added `BOOST_ALLOW_DEPRECATED_HEADERS` definition for Boost 1.89.0 compatibility

#### Build Warnings (Non-Critical)

The following warnings appear during compilation but do not affect functionality:
- **MicIO.h**: Non-void functions not returning values in template stubs (lines 207, 211)
- **IteratorAdapter.h**: Const qualifier on reference types has no effect (lines 166, 252)

These are legacy code issues that do not impact the Release build.

### Running IceNine

```bash
./IceNine [mode] <ConfigFile>
```

**Modes**:
- `s` - Simulate detector images from sample orientations
- `r` - Reconstruct sample orientations from detector data
- `p` - Parallel reconstruction (distributed MPI)
- `G` - Paint grid visualization

**Example**:
```bash
./IceNine r ConfigFiles/ReconstructTest.config
```

## Architecture

### Core Components

**Layered Architecture**:
```
Application Layer (Driver.cpp, main.cpp)
    ↓
Reconstruction Layer (Reconstructor, BreadthFirstReconstructor, ReconstructionStrategies)
    ↓
Search & Optimization Layer (DiscreteAdaptive, DiscreteSearch, CostFunctions)
    ↓
Simulation Layer (Simulation, ForwardSimulation)
    ↓
Data Layer (Sample, Detector, ImageData, Peak)
    ↓
Utility/Physics Layer (XDM++, DiffractionCore)
```

### Key Workflows

**Forward Simulation**:
```
ConfigFile → ExperimentSetup → Sample/Detector →
Simulation → ForwardSimulation → DetectorImage output
```

**Reconstruction**:
```
ConfigFile → ExperimentSetup → Sample/DetectorImage →
Reconstructor → BreadthFirstReconstructor →
  ├── DiscreteAdaptive (orientation search)
  ├── CostFunctions (matching metrics)
  └── ReconstructionStrategies (boundary propagation)
→ Reconstructed Sample (.mic output)
```

**Parallel Reconstruction**:
- Master process (XDMServer) coordinates work distribution
- Worker processes (XDMClient) perform local reconstruction
- Communication via asynchronous MPI (XDMParallel, AsynchronousMPI)

### Critical Module Relationships

- **Reconstructor** orchestrates reconstruction and delegates to:
  - **DiscreteAdaptive**: Implements adaptive orientation search with hierarchical refinement
  - **CostFunctions**: Computes cost/confidence metrics comparing simulated vs. observed peaks
  - **Simulation**: Generates forward model predictions for candidate orientations
  - **ReconstructionStrategies**: Selects boundary voxels for propagation

- **DiscreteAdaptive** combines:
  - **DiscreteSearch**: Coarse sampling of orientation space
  - **ContinuousSearch**: Fine local optimization (gradient-based)
  - **SearchDetails**: Tracks search statistics and convergence

- **CostFunctions** uses:
  - **DiffractionCore**: Physics calculations (Bragg angles, scattering vectors)
  - **OverlapInfo**: Geometric overlap between predicted and observed peaks

## Important Design Patterns

### Template-Based Composition
- Extensive use of C++ templates for algorithm composition without runtime overhead
- `.tmpl.cpp` files contain template implementations separate from declarations
- Example: [SearchTraits.h](Src/SearchTraits.h) composes search strategies via templates

### Configuration-Driven Execution
All runtime behavior controlled via config files with sections for:
- Input/Output paths and file formats
- Beam parameters (energy, direction, height)
- Detector geometry and calibration
- Sample voxel grid and crystal structures
- Search parameters (resolution levels, candidate counts)
- Optimization parameters (cost thresholds, convergence criteria)

See [ConfigFiles/ReconstructTest.config](ConfigFiles/ReconstructTest.config) for a complete example.

### Voxel-Based Spatial Representation
- Samples discretized as 3D voxel grids stored in `.mic` format
- Each voxel contains crystal orientation (quaternion representation)
- Boundary voxels selected for adaptive reconstruction via [BoundarySelectionStrategies.h](Src/BoundarySelectionStrategies.h)

## Key Source Files

### Entry Points
- [Src/main.cpp](Src/main.cpp) - Main entry point, calls `XDMDriver()`
- [Src/Driver.cpp](Src/Driver.cpp) - Application driver and mode dispatcher

### Core Reconstruction
- [Src/Reconstructor.h](Src/Reconstructor.h) - Main reconstruction interface
- [Src/BreadthFirstReconstructor.h](Src/BreadthFirstReconstructor.h) - Breadth-first propagation strategy
- [Src/ReconstructionStrategies.h](Src/ReconstructionStrategies.h) - Neighbor selection algorithms
- [Src/SerialReconstruction.h](Src/SerialReconstruction.h) - Single-process reconstruction

### Search & Optimization
- [Src/DiscreteAdaptive.h](Src/DiscreteAdaptive.h) - Adaptive orientation search (main algorithm)
- [Src/DiscreteSearch.h](Src/DiscreteSearch.h) - Discrete candidate sampling
- [Src/ContinuousSearch.h](Src/ContinuousSearch.h) - Continuous optimization
- [Src/CostFunctions.h](Src/CostFunctions.h) - Cost/confidence metrics
- [Src/SearchDetails.h](Src/SearchDetails.h) - Search statistics and tracking

### Simulation
- [Src/Simulation.h](Src/Simulation.h) - Forward simulation engine
- [Src/ForwardSimulation.h](Src/ForwardSimulation.h) - Forward model implementation
- [Src/DiffractionCore.h](Src/DiffractionCore.h) - Diffraction physics (inline calculations)

### Data Structures
- [Src/Sample.h](Src/Sample.h) - Sample geometry and voxel data
- [Src/Detector.h](Src/Detector.h) - Detector geometry and properties
- [Src/ImageData.h](Src/ImageData.h) - Detector image representation
- [Src/Peak.h](Src/Peak.h) - Diffraction peak data structure
- [Src/ExperimentSetup.h](Src/ExperimentSetup.h) - Experimental parameters container

### Parallel Processing
- [Src/XDMServer.h](Src/XDMServer.h) - Master-worker server
- [Src/XDMClient.h](Src/XDMClient.h) - Worker client
- [Src/XDMParallel.h](Src/XDMParallel.h) - Parallel utilities
- [Src/ParallelDriver.h](Src/ParallelDriver.h) - Parallel execution orchestration

### Configuration & I/O
- [Src/ConfigFile.h](Src/ConfigFile.h) - Configuration file parsing
- [Src/InitFilesIO.h](Src/InitFilesIO.h) - Initial file I/O utilities

### XDM++ Library (XDM++/libXDM/)
- [3dMath.h](XDM++/libXDM/3dMath.h) - 3D vector/matrix math
- [Quaternion.h](XDM++/libXDM/Quaternion.h) - Quaternion operations
- [CrystalStructure.h](XDM++/libXDM/CrystalStructure.h) - Crystal structure definitions
- [MicIO.h](XDM++/libXDM/MicIO.h) - `.mic` file format I/O
- [MicGrid.h](XDM++/libXDM/MicGrid.h) - Voxel grid representation
- [Symmetry.h](XDM++/libXDM/Symmetry.h) - Crystal symmetry operations
- [AsynchronousMPI.h](XDM++/libXDM/AsynchronousMPI.h) - Async MPI communication

## File Formats

### Input Files
- **Config files** (`.config`): Text-based configuration (see ConfigFiles/)
- **Sample files** (`.mic`): Voxel grid with orientations
- **Detector files** (`.txt`): Detector geometry specification
- **Structure files** (`.dat`): Crystal structure parameters
- **Peak files** (`.bin`, `.txt`): Experimental diffraction data

### Output Files
- **Reconstructed samples** (`.mic`): Voxel grid with reconstructed orientations
- **Reduced data** (`.d`): Processed detector data

## Development Notes

### Code Style
- Heavy use of C++ templates for generic programming
- Template implementations in `.tmpl.cpp` files
- Composition over inheritance pattern throughout
- Boost libraries used extensively (ublas, MPI, serialization)

### Performance Considerations
- **Always use Release builds** - debug builds are ~100x slower
- Parallel reconstruction recommended for large samples
- Cost function evaluation is the computational bottleneck
- DiffractionCore uses inline functions for performance

### Common Modification Points
- **Adding new search strategies**: Extend [SearchTraits.h](Src/SearchTraits.h)
- **Modifying cost functions**: Edit [CostFunctions.h](Src/CostFunctions.h)
- **Changing boundary selection**: Modify [BoundarySelectionStrategies.h](Src/BoundarySelectionStrategies.h)
- **Adjusting reconstruction flow**: Update [ReconstructionStrategies.h](Src/ReconstructionStrategies.h)

### Testing
- Example configurations in [ConfigFiles/](ConfigFiles/)
- Test data typically in `DataFiles/` or `TestInput/` directories (not included in repo)
- Use small sample files for quick testing before scaling to full datasets
