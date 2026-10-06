# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

IceNine is a forward model reconstruction tool for synchrotron X-ray diffraction data for three-dimensional orientation imaging of polycrystalline materials. It reconstructs crystal grain orientations from experimental detector images using adaptive search algorithms.

**Citation**: S. F. Li, and R. M. Suter, "Adaptive reconstruction method for three-dimensional orientation imaging", Journal of Applied Crystallography, 2013.

## Build & Run (C++)

### Build

```bash
brew install cmake boost eigen open-mpi          # macOS deps
export PATH="/opt/homebrew/bin:$PATH"
cmake -DCMAKE_BUILD_TYPE=Release . && make -j8
```

**CRITICAL**: Always use `CMAKE_BUILD_TYPE=Release` — debug builds are ~100x slower.

### Run

```bash
./IceNine [mode] <ConfigFile>
```

Modes: `s` (simulate), `r` (reconstruct), `p` (parallel/MPI reconstruct), `G` (paint grid)

Example: `./IceNine r ConfigFiles/ReconstructTest.config`

## Python Port (`icenine_py/`)

Ongoing test-driven migration of IceNine's core physics to Python/PyTorch. Every function is validated against C++ ground truth (tolerance < 1e-6). C++ ground truth data lives in `icenine_py/cpp_outputs/` (JSON). See [icenine_py/MIGRATION_HISTORY.md](icenine_py/MIGRATION_HISTORY.md) for detailed migration record.

### Setup & Test

**CRITICAL**: All Python commands MUST be run using `uv`. Never use bare `python3`, `pip`, or `pytest`.

```bash
cd icenine_py
uv pip install -e ".[dev]"              # install with dev deps (pytest, black, mypy)
uv run pytest tests/                     # run all tests
uv run pytest tests/test_symmetry.py -v  # run a single test file
uv run pytest --cov=icenine tests/       # with coverage
uv run python script.py                  # run any Python script
```

### Python Style
- Black formatter, line length 100, target Python 3.9+
- Type annotations required (mypy strict)
- Dependencies: numpy, torch, pymatgen, scipy (see `pyproject.toml`)

### Critical Gotchas (Python Port)
- **Euler angle conventions**: Voxel orientations use active ZXZ Bunge convention; global sample orientation uses **passive** Euler convention (`passive_euler_matrix()`). These are NOT the same.
- **Degree/radian boundary**: `.mic` files store angles in degrees; conversion to radians happens at I/O boundaries. Config files have 15 parameters requiring degree-to-radian conversion.
- **Batching required**: Single-element physics operations are too slow. All diffraction calculations must use batched PyTorch operations.
- **Spatial indexing**: Always use `scipy.spatial.cKDTree` (not `KDTree`) — 10-100x faster, identical API.
- **Q-max mismatch is intentional**: Simulation data may use Q-max=16 while reconstruction uses Q-max=8. This does NOT degrade reconstruction quality. The cost function does not penalize under-observed peaks — peaks beyond max_q are filtered out at `VoxelCostFunction.__init__` time. Both C++ and Python handle this identically. Fewer reciprocal vectors means fewer peaks in the metric (faster but less discriminating).
- **Cost function angular sharpness**: The cost function (pixel overlap quality) is extremely sharp in orientation space — quality drops rapidly within 0.2–0.5 degrees of the correct orientation. This is fundamental Bragg diffraction physics (peaks are narrow in angle). Multi-level adaptive search is necessary: coarse Sukharev grid to find the basin, then fine MC refinement within it.

## Architecture

```
Application Layer (Driver.cpp, main.cpp)
    ↓
Reconstruction Layer (Reconstructor, BreadthFirstReconstructor, ReconstructionStrategies)
    ↓
Search & Optimization (DiscreteAdaptive, DiscreteSearch, CostFunctions)
    ↓
Simulation Layer (Simulation, ForwardSimulation)
    ↓
Data Layer (Sample, Detector, ImageData, Peak)
    ↓
Utility/Physics (XDM++/libXDM/, DiffractionCore)
```

### Key Workflows

**Forward Simulation**: ConfigFile → ExperimentSetup → Sample/Detector → Simulation → ForwardSimulation → DetectorImage

**Reconstruction**: ConfigFile → ExperimentSetup → Sample/DetectorImage → Reconstructor → BreadthFirstReconstructor → (DiscreteAdaptive + CostFunctions + ReconstructionStrategies) → Reconstructed .mic

**Parallel Reconstruction**: Master (XDMServer) distributes work to workers (XDMClient) via async MPI.

### Critical Module Relationships

- **Reconstructor** orchestrates: DiscreteAdaptive (orientation search), CostFunctions (matching metrics), Simulation (forward model), ReconstructionStrategies (boundary propagation)
- **DiscreteAdaptive** combines: DiscreteSearch (coarse sampling) + ContinuousSearch (fine optimization) + SearchDetails (convergence tracking)
- **CostFunctions** uses: DiffractionCore (Bragg angles, scattering vectors) + OverlapInfo (geometric overlap)

## Design Patterns

### Template-Based Composition (C++)
- Extensive C++ templates for zero-overhead algorithm composition
- `.tmpl.cpp` files contain template implementations separate from declarations
- [SearchTraits.h](Src/SearchTraits.h) composes search strategies via template parameters

### Configuration-Driven Execution
All runtime behavior controlled via `.config` files (beam parameters, detector geometry, sample grid, search/optimization parameters). See [ConfigFiles/ReconstructTest.config](ConfigFiles/ReconstructTest.config).

### Voxel-Based Spatial Representation
Samples are 3D voxel grids in `.mic` format; each voxel stores crystal orientation as a quaternion. Boundary voxels selected for adaptive reconstruction via [BoundarySelectionStrategies.h](Src/BoundarySelectionStrategies.h).

## Common Modification Points

- **Search strategies**: Extend [SearchTraits.h](Src/SearchTraits.h)
- **Cost functions**: Edit [CostFunctions.h](Src/CostFunctions.h)
- **Boundary selection**: Modify [BoundarySelectionStrategies.h](Src/BoundarySelectionStrategies.h)
- **Reconstruction flow**: Update [ReconstructionStrategies.h](Src/ReconstructionStrategies.h)

## Build System Notes

- CMake 3.15+ with `find_package()` for MPI, Boost, Eigen3 (no hardcoded paths)
- Eigen3 forced to Homebrew path, ignoring legacy `3rdParty/` directory
- C++11 required; Boost Lambda replaced with C++11 lambdas
- `-Wno-deprecated-declarations` and `BOOST_ALLOW_DEPRECATED_HEADERS` for Boost 1.89+ compat
- Known benign warnings: MicIO.h template stubs, IteratorAdapter.h const qualifiers
- See [CPP_MODERNIZATION.md](CPP_MODERNIZATION.md) for full build/dependency change history

## File Formats

- **Config** (`.config`): Text-based runtime configuration
- **Sample** (`.mic`): Voxel grid with crystal orientations (input/output)
- **Detector** (`.txt`): Detector geometry
- **Structure** (`.dat`): Crystal structure parameters
- **Peaks** (`.bin`, `.txt`): Experimental diffraction data
- **Reduced data** (`.d`): Processed detector data

## Performance

- Cost function evaluation is the computational bottleneck
- DiffractionCore uses inline functions for performance-critical physics
- Parallel reconstruction (mode `p`) recommended for large samples

## Gitflow Workflow

This project uses gitflow. All work must follow this branching model:

```
master           ← releases only (tagged, e.g. v3.last)
  └─ develop     ← integration branch, PRs merge here
       └─ feature/name           ← feature branches off develop
            └─ feature/name/task ← optional sub-task branches
```

### Branch Rules
- **master**: Never commit directly. Only merge from develop for releases.
- **develop**: Integration branch. Feature branches merge here via PR.
- **feature/***: All new work. Branch off develop, merge back via PR.
- **feature/\*/task**: Sub-tasks off a feature branch. Merge back to parent feature.

### Workflow Commands
- `/start-feature <name>` — create feature branch off develop
- `/start-task <name>` — create sub-task branch off current feature
- `/finish-task` — run tests, review, merge task → parent feature
- `/finish-feature` — run tests, review, create PR → develop, merge via GitHub

### Merge Rules
- **Feature → develop**: ALWAYS via GitHub PR. Never merge locally and push develop. Push the feature branch, create PR with `gh pr create`, merge with `gh pr merge --merge --delete-branch`.
- **Task → parent feature**: Local merge with `--no-ff` is fine (no PR needed).
- Always run the relevant test suite before merging:
  - Python changes: `cd icenine_py && uv run pytest tests/ -v`
  - C++ changes: `cmake -DCMAKE_BUILD_TYPE=Release . && make -j8`

### Documentation Rules
These rules apply to ALL work, not just when using slash commands:

1. **Plan documentation**: When entering plan mode, save the plan to `icenine_py/MIGRATION_HISTORY.md` under a new section (if it relates to the Python port) or the appropriate doc.
2. **Completion summaries**: After executing a plan, summarize what was accomplished in `icenine_py/MIGRATION_HISTORY.md`, then delete the plan file from `.claude/plans/`.
3. **Feature documentation**: When a new feature or module is added, update `icenine_py/README.md` to reflect it (new modules, changed structure, new dependencies, etc.).
