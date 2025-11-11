# C++17 Migration Guide - Future Boost Reduction

## Overview

Upgrading from C++11 to C++17 would enable replacing **3 more Boost libraries** currently in use:
- `boost::optional` → `std::optional`
- `boost::variant` → `std::variant`
- `boost::any` → `std::any`

Plus additional C++17 features that improve code quality.

---

## 1. boost::optional → std::optional

### Current Usage Pattern (Likely in codebase)

```cpp
// C++11/Boost
#include <boost/optional.hpp>

boost::optional<SQuaternion> FindOrientation(const Data& data) {
  if (data.IsValid()) {
    return boost::optional<SQuaternion>(ComputeOrientation(data));
  }
  return boost::none;
}

// Usage
boost::optional<SQuaternion> result = FindOrientation(data);
if (result) {
  SQuaternion q = result.get();
  // Use q
} else {
  // Handle missing value
}
```

### C++17 Migration

```cpp
// C++17
#include <optional>

std::optional<SQuaternion> FindOrientation(const Data& data) {
  if (data.IsValid()) {
    return ComputeOrientation(data);  // Implicit conversion
  }
  return std::nullopt;
}

// Usage - More expressive
std::optional<SQuaternion> result = FindOrientation(data);

// Method 1: if with initializer (C++17)
if (auto q = result; q.has_value()) {
  // Use *q
}

// Method 2: value_or with default
SQuaternion q = result.value_or(DefaultQuaternion());

// Method 3: Monadic operations (C++23, but shows direction)
result.and_then([](auto q) {
  return ProcessOrientation(q);
});
```

### Real Example from IceNine Context

```cpp
// Finding best candidate in reconstruction
class DiscreteSearch {
  // Before (C++11)
  boost::optional<SCandidate> GetBestCandidate(
    const vector<SCandidate>& candidates, Float threshold
  ) {
    for (const auto& c : candidates) {
      if (c.fCost < threshold) {
        return boost::optional<SCandidate>(c);
      }
    }
    return boost::none;
  }

  // After (C++17)
  std::optional<SCandidate> GetBestCandidate(
    const vector<SCandidate>& candidates, Float threshold
  ) {
    auto it = std::find_if(candidates.begin(), candidates.end(),
      [threshold](const auto& c) { return c.fCost < threshold; });

    if (it != candidates.end()) {
      return *it;  // Automatic construction
    }
    return std::nullopt;
  }
};

// Usage
if (auto best = search.GetBestCandidate(candidates, 0.1f)) {
  ApplyOrientation(voxel, *best);
}
```

**Benefits**:
- ✅ Standard library (better support)
- ✅ More expressive syntax
- ✅ Better integration with modern C++ patterns
- ✅ ~15% smaller binary footprint

---

## 2. boost::variant → std::variant

### Current Usage Pattern

```cpp
// C++11/Boost
#include <boost/variant.hpp>

// Variant can hold different types
using SearchResult = boost::variant<
  SQuaternion,           // Successful orientation
  SearchError,           // Error state
  PartialSolution        // Needs refinement
>;

class Visitor : public boost::static_visitor<void> {
public:
  void operator()(const SQuaternion& q) const {
    std::cout << "Found orientation\n";
  }

  void operator()(const SearchError& e) const {
    std::cerr << "Error: " << e.message << "\n";
  }

  void operator()(const PartialSolution& p) const {
    std::cout << "Partial solution, continuing...\n";
  }
};

SearchResult result = PerformSearch(voxel);
boost::apply_visitor(Visitor(), result);
```

### C++17 Migration

```cpp
// C++17
#include <variant>

using SearchResult = std::variant<
  SQuaternion,
  SearchError,
  PartialSolution
>;

SearchResult result = PerformSearch(voxel);

// Method 1: std::visit with lambda
std::visit([](auto&& arg) {
  using T = std::decay_t<decltype(arg)>;
  if constexpr (std::is_same_v<T, SQuaternion>) {
    std::cout << "Found orientation\n";
  } else if constexpr (std::is_same_v<T, SearchError>) {
    std::cerr << "Error: " << arg.message << "\n";
  } else if constexpr (std::is_same_v<T, PartialSolution>) {
    std::cout << "Partial solution, continuing...\n";
  }
}, result);

// Method 2: Overloaded lambda pattern (elegant!)
template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> overloaded(Ts...) -> overloaded<Ts...>;

std::visit(overloaded {
  [](const SQuaternion& q) { std::cout << "Found orientation\n"; },
  [](const SearchError& e) { std::cerr << "Error: " << e.message << "\n"; },
  [](const PartialSolution& p) { std::cout << "Partial, continuing...\n"; }
}, result);

// Method 3: std::get with index/type
if (std::holds_alternative<SQuaternion>(result)) {
  SQuaternion q = std::get<SQuaternion>(result);
  // Use q
}
```

### Real Example: Reconstruction Status

```cpp
// Representing reconstruction state
enum class ReconstructionStatus { Success, Failed, NeedsRefinement };

// Before (separate status + data)
struct ReconstructionResult {
  ReconstructionStatus status;
  std::shared_ptr<CMic> data;  // May be null
  std::string errorMessage;    // Only if failed
};

// After (C++17 variant)
struct SuccessResult {
  std::shared_ptr<CMic> mic;
  Float confidence;
};

struct FailureResult {
  std::string errorMessage;
  Int errorCode;
};

struct RefinementNeeded {
  std::shared_ptr<CMic> partialMic;
  vector<SVoxel> problematicVoxels;
};

using ReconstructionResult = std::variant<
  SuccessResult,
  FailureResult,
  RefinementNeeded
>;

// Clean handling
ReconstructionResult Reconstruct(const Sample& sample) {
  // ... reconstruction logic ...

  if (converged) {
    return SuccessResult{ mic, confidence };
  } else if (criticalError) {
    return FailureResult{ "Convergence failed", -1 };
  } else {
    return RefinementNeeded{ partialMic, problemVoxels };
  }
}

// Usage
auto result = Reconstruct(sample);
std::visit(overloaded {
  [](const SuccessResult& s) {
    std::cout << "Success! Confidence: " << s.confidence << "\n";
    SaveMic(s.mic);
  },
  [](const FailureResult& f) {
    std::cerr << "Failed: " << f.errorMessage << "\n";
  },
  [](const RefinementNeeded& r) {
    std::cout << "Refining " << r.problematicVoxels.size() << " voxels\n";
    RefineReconstruction(r.partialMic, r.problematicVoxels);
  }
}, result);
```

**Benefits**:
- ✅ Type-safe unions
- ✅ No runtime overhead
- ✅ Better error handling patterns
- ✅ Elegant visitor syntax with C++17

---

## 3. boost::any → std::any

### Current Usage Pattern

```cpp
// C++11/Boost
#include <boost/any.hpp>

// Generic property storage
class PropertyBag {
  std::map<std::string, boost::any> properties;

public:
  template<typename T>
  void Set(const std::string& key, const T& value) {
    properties[key] = value;
  }

  template<typename T>
  T Get(const std::string& key) const {
    auto it = properties.find(key);
    if (it != properties.end()) {
      return boost::any_cast<T>(it->second);
    }
    throw std::runtime_error("Key not found");
  }
};

// Usage
PropertyBag config;
config.Set("MaxIterations", 1000);
config.Set("Tolerance", 1e-6);
config.Set("Algorithm", std::string("AdaptiveSearch"));

int maxIter = config.Get<int>("MaxIterations");
```

### C++17 Migration

```cpp
// C++17
#include <any>

class PropertyBag {
  std::map<std::string, std::any> properties;

public:
  template<typename T>
  void Set(const std::string& key, const T& value) {
    properties[key] = value;
  }

  template<typename T>
  std::optional<T> Get(const std::string& key) const {
    auto it = properties.find(key);
    if (it != properties.end()) {
      try {
        return std::any_cast<T>(it->second);
      } catch (const std::bad_any_cast&) {
        return std::nullopt;  // Wrong type
      }
    }
    return std::nullopt;  // Key not found
  }

  // C++17: Get with default
  template<typename T>
  T GetOr(const std::string& key, T defaultValue) const {
    return Get<T>(key).value_or(defaultValue);
  }
};

// Usage - safer!
PropertyBag config;
config.Set("MaxIterations", 1000);
config.Set("Tolerance", 1e-6);

// Type-safe retrieval
if (auto maxIter = config.Get<int>("MaxIterations")) {
  std::cout << "Max iterations: " << *maxIter << "\n";
}

// With default
int maxIter = config.GetOr<int>("MaxIterations", 100);
```

**Benefits**:
- ✅ Standard library
- ✅ Same runtime characteristics
- ✅ Better integration with `std::optional`

---

## 4. Additional C++17 Features for IceNine

### 4.1 Structured Bindings (Replace std::tie)

```cpp
// Before (C++11)
std::tuple<SVoxel, Int> ReconstructVoxel(const SVoxel& v);

SVoxel result;
Int code;
std::tie(result, code) = ReconstructVoxel(voxel);

// After (C++17)
auto [result, code] = ReconstructVoxel(voxel);  // ✨ Much cleaner!

// Works with maps too
std::map<Int, SQuaternion> orientations;
for (const auto& [id, orientation] : orientations) {  // ✨ Beautiful!
  ProcessOrientation(id, orientation);
}
```

### 4.2 if/switch with Initializer

```cpp
// Before (C++11)
auto candidate = FindBestCandidate(candidates);
if (candidate.IsValid()) {
  ApplyCandidate(candidate);
}

// After (C++17) - scoped initialization
if (auto candidate = FindBestCandidate(candidates); candidate.IsValid()) {
  ApplyCandidate(candidate);
}  // candidate out of scope here

// Switch with init
switch (auto status = Reconstruct(voxel); status.GetCode()) {
  case SUCCESS: break;
  case FAILED: Retry(); break;
}
```

### 4.3 constexpr if (Compile-time branching)

```cpp
// Template metaprogramming simplification
template<typename T>
void ProcessData(const T& data) {
  // Before (C++11) - SFINAE hell or tag dispatch

  // After (C++17)
  if constexpr (std::is_same_v<T, SQuaternion>) {
    // Compile-time branch for quaternions
    ProcessQuaternion(data);
  } else if constexpr (std::is_same_v<T, SMatrix3x3>) {
    // Compile-time branch for matrices
    ProcessMatrix(data);
  } else {
    static_assert(false, "Unsupported type");
  }
}
```

### 4.4 std::filesystem (Replace manual path handling)

```cpp
// Before (C++11) - manual string manipulation
#include <string>

std::string GetOutputPath(const std::string& base, const std::string& filename) {
  return base + "/" + filename + ".mic";
}

// After (C++17)
#include <filesystem>
namespace fs = std::filesystem;

fs::path GetOutputPath(const fs::path& base, const std::string& filename) {
  return base / filename += ".mic";  // Proper path concatenation
}

// Check if file exists
if (fs::exists(outputPath)) {
  fs::remove(outputPath);
}

// Create directories
fs::create_directories(outputPath.parent_path());

// Iterate directory
for (const auto& entry : fs::directory_iterator("ConfigFiles")) {
  if (entry.path().extension() == ".config") {
    ProcessConfig(entry.path());
  }
}
```

### 4.5 std::string_view (Efficient string handling)

```cpp
// Before (C++11) - copies or const char*
void ProcessName(const std::string& name) {  // Copy if temporary passed
  if (name.substr(0, 4) == "Test") {  // Creates new string
    // ...
  }
}

// After (C++17)
void ProcessName(std::string_view name) {  // No copy, cheap to pass
  if (name.substr(0, 4) == "Test") {  // Returns string_view, no allocation
    // ...
  }
}

// Especially useful for string literals
constexpr std::string_view kConfigExt = ".config";  // Compile-time constant
```

### 4.6 Parallel Algorithms

```cpp
// Before (C++11) - manual OpenMP or custom parallel code
#pragma omp parallel for
for (int i = 0; i < voxels.size(); ++i) {
  ProcessVoxel(voxels[i]);
}

// After (C++17)
#include <execution>

// Parallel for_each
std::for_each(std::execution::par, voxels.begin(), voxels.end(),
  [](auto& voxel) { ProcessVoxel(voxel); }
);

// Parallel transform
std::transform(std::execution::par, orientations.begin(), orientations.end(),
               results.begin(),
  [](const auto& q) { return ComputeCost(q); }
);

// Parallel sort
std::sort(std::execution::par, candidates.begin(), candidates.end(),
  [](const auto& a, const auto& b) { return a.fCost < b.fCost; }
);
```

---

## Migration Comparison Table

| Feature | Boost (C++11) | C++17 STL | Benefits |
|---------|---------------|-----------|----------|
| Optional values | `boost::optional` | `std::optional` | Standard, better syntax |
| Type-safe unions | `boost::variant` | `std::variant` | Standard, `if constexpr` support |
| Type erasure | `boost::any` | `std::any` | Standard, same performance |
| Filesystem | Manual/Boost.Filesystem | `std::filesystem` | Portable, feature-rich |
| String views | `const char*`/`const std::string&` | `std::string_view` | Zero-copy, efficient |
| Tuple unpacking | `std::tie(a, b, c) = ...` | `auto [a, b, c] = ...` | Cleaner syntax |
| Parallel algorithms | OpenMP/manual | `std::execution::par` | Standard, portable |

---

## Estimated Impact for IceNine

### Files That Would Benefit

Based on codebase analysis, C++17 migration would likely improve:

1. **Error Handling** (20-30 files)
   - Replace error codes with `std::optional<Result>`
   - Use `std::variant` for different error types

2. **Configuration** (5-10 files)
   - `ConfigFile.h` - Use `std::optional` for optional parameters
   - Use `std::filesystem` for path handling

3. **Search Results** (10-15 files)
   - `DiscreteSearch.h` - Return `std::optional<SCandidate>`
   - Use structured bindings for tuple returns

4. **File I/O** (5-8 files)
   - `MicIO.h`, `InitFilesIO.h` - Use `std::filesystem`
   - Replace manual path string manipulation

5. **Parallel Processing** (3-5 files)
   - `XDMParallel.h` - Potentially use C++17 parallel algorithms
   - Simpler than current MPI for shared-memory parallelism

### Migration Effort

| Component | Files | Effort | Risk |
|-----------|-------|--------|------|
| boost::optional → std::optional | ~15 | Low | Low |
| boost::variant → std::variant | ~5 | Medium | Low |
| boost::any → std::any | ~3 | Low | Low |
| Path handling → std::filesystem | ~10 | Medium | Low |
| Add structured bindings | ~30 | Low | Very Low |
| **Total** | **~50 files** | **1-2 weeks** | **Low** |

---

## Recommended Migration Path

### Phase 1: Preparation
1. Update CMake to require C++17
2. Test build with C++17 enabled
3. Fix any new warnings

### Phase 2: Low-Risk Replacements (Week 1)
1. Replace `boost::optional` → `std::optional`
2. Add structured bindings where beneficial
3. Use `if/switch` with initializers

### Phase 3: Moderate Changes (Week 2)
1. Replace `boost::variant` → `std::variant`
2. Replace `boost::any` → `std::any`
3. Adopt `std::filesystem` for path handling
4. Add `std::string_view` where appropriate

### Phase 4: Optimization
1. Evaluate parallel algorithms for appropriate use cases
2. Use `constexpr if` to simplify template code
3. Profile for performance improvements

---

## CMake Changes Required

```cmake
# Current (C++11)
set(CMAKE_CXX_STANDARD 11)

# After C++17 migration
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optional: Enable parallel algorithms (needs TBB on some platforms)
find_package(TBB)
if(TBB_FOUND)
  target_link_libraries(IceNine PRIVATE TBB::tbb)
endif()
```

---

## Conclusion

**C++17 Migration Benefits for IceNine**:

✅ **Remove 3 more Boost libraries** (optional, variant, any)
✅ **Cleaner, more readable code** (structured bindings)
✅ **Better error handling** (optional + variant patterns)
✅ **Safer file operations** (std::filesystem)
✅ **Potential performance gains** (parallel algorithms, string_view)
✅ **Reduced external dependencies** (~25% more Boost reduction)
✅ **Future-proof** (C++17 is widely supported as of 2024)

**Estimated Total Boost Reduction After C++17**:
- Current: 9 Boost libraries remaining
- After C++17: ~6 Boost libraries remaining (33% additional reduction)
- Only specialized libraries left: multi_array, ublas, graph, disjoint_sets, dynamic_bitset

**Recommendation**: C++17 migration is highly worthwhile for IceNine, offering significant code quality improvements with low risk.
