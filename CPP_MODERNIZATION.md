# C++ Build & Dependency Modernization

Consolidated record of C++ build system and dependency changes (November 2024).

## CMake Modernization (2.4 → 3.15+)

- Migrated to modern target-based CMake with `find_package()` discovery for MPI, Boost, Eigen3
- Removed all hardcoded paths; dependencies resolved via Homebrew on macOS
- Eigen3 forced to Homebrew path, explicitly ignoring legacy `3rdParty/` directory
- Added C++11 standard requirement

**Files modified**: `CMakeLists.txt`, `XDM++/libXDM/CMakeLists.txt`

## Boost → C++11 STL Migration

Replaced 5 Boost libraries with STL equivalents:

| Boost Library | Replacement | Files Changed |
|---------------|-------------|---------------|
| `boost::shared_ptr` | `std::shared_ptr` | 15 files |
| `boost::tuple` | `std::tuple` | 5 files |
| `boost::function` | `std::function` | 3 files |
| `boost::random` | `<random>` | 2 files |
| `boost::lambda` | C++11 lambdas | 2 files |

**9 Boost libraries remain** (ublas, MPI, serialization, etc.) with no direct C++11 equivalent.

## C++11 Compatibility Fixes

- `make_pair()`: Removed explicit template arguments (C++11 auto-deduction)
- `const ShapePtr&`: Changed to value to avoid temporary binding issues
- Compiler flags: `-Wno-deprecated-declarations`, `BOOST_ALLOW_DEPRECATED_HEADERS`

## Known Benign Warnings

- `MicIO.h`: Non-void template stubs not returning values (lines 207, 211)
- `IteratorAdapter.h`: Const qualifier on reference types (lines 166, 252)

## RNG Note

`GetRandomLocalGrid()` creates a local unseeded `mt19937` per call, producing deterministic output. This behavior was preserved from the Boost version — may be intentional.

## Future: C++17 Migration

See [CPP17_MIGRATION_GUIDE.md](CPP17_MIGRATION_GUIDE.md) for optional further Boost reduction (~3 more libraries removable with C++17 features like `std::optional`, `std::variant`, `if constexpr`). Estimated effort: 1-2 weeks.
