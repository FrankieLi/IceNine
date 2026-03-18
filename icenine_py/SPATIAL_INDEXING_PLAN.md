# MicFile Spatial Indexing - Implementation Plan

**Date**: 2025-11-12
**Phase**: 3 - Spatial Queries and Neighbor Finding

## Objective

Add spatial indexing capabilities to MicFile for efficient neighbor queries, essential for:
- Reconstruction algorithms (boundary propagation)
- Grain boundary analysis
- Orientation gradient calculations
- Spatial filtering and smoothing

## Requirements from C++ Implementation

### From MicGrid.h
```cpp
vector<ShapePtr> GetNeighbors(const ShapePtr & pShape) const;
```

### From Reconstruction Code
- Find boundary voxels (fitted neighbors of unfitted voxels)
- Find k-nearest neighbors
- Find voxels within radius
- Determine if voxel is on boundary

## Proposed Python API

### 1. Spatial Index Creation
```python
class MicFile:
    def build_spatial_index(self) -> None:
        """Build KDTree for fast spatial queries."""
        self._kdtree = KDTree(self.positions.numpy())
        self._built_index = True

    def _ensure_spatial_index(self) -> None:
        """Ensure spatial index is built before queries."""
        if not hasattr(self, '_built_index') or not self._built_index:
            self.build_spatial_index()
```

### 2. Neighbor Queries
```python
def get_neighbors(
    self,
    voxel_idx: int,
    radius: float,
    max_neighbors: Optional[int] = None
) -> List[int]:
    """
    Find neighboring voxels within radius.

    Args:
        voxel_idx: Index of query voxel
        radius: Search radius in meters
        max_neighbors: Maximum number of neighbors to return

    Returns:
        List of voxel indices within radius
    """

def get_k_nearest_neighbors(
    self,
    voxel_idx: int,
    k: int
) -> Tuple[List[int], List[float]]:
    """
    Find k nearest neighboring voxels.

    Args:
        voxel_idx: Index of query voxel
        k: Number of neighbors to find

    Returns:
        Tuple of (neighbor_indices, distances)
    """
```

### 3. Boundary Detection
```python
def get_boundary_voxels(
    self,
    fitted_phase: Optional[int] = None,
    unfitted_phase: int = 0
) -> List[int]:
    """
    Find boundary voxels (fitted voxels with unfitted neighbors).

    Args:
        fitted_phase: Phase ID for fitted voxels (None = any phase > 0)
        unfitted_phase: Phase ID for unfitted voxels (default: 0)

    Returns:
        List of voxel indices on boundary
    """

def is_boundary_voxel(
    self,
    voxel_idx: int,
    radius: Optional[float] = None
) -> bool:
    """
    Check if voxel is on boundary (has unfitted neighbors).

    Args:
        voxel_idx: Index of voxel to check
        radius: Search radius (default: 2 * side_length)

    Returns:
        True if voxel is on boundary
    """
```

### 4. Spatial Queries (Advanced)
```python
def query_region(
    self,
    center: np.ndarray,
    radius: float
) -> List[int]:
    """
    Find all voxels within radius of a point.

    Args:
        center: 3D point (x, y, z) in meters
        radius: Search radius in meters

    Returns:
        List of voxel indices in region
    """

def get_orientation_gradient(
    self,
    voxel_idx: int,
    k_neighbors: int = 6
) -> torch.Tensor:
    """
    Calculate orientation gradient using neighboring voxels.

    Args:
        voxel_idx: Index of voxel
        k_neighbors: Number of neighbors to use

    Returns:
        Orientation gradient (misorientation with neighbors)
    """
```

## Implementation Strategy

### Phase 3.1: Basic Spatial Indexing
- [x] Add `scipy.spatial.KDTree` dependency (already in requirements.txt)
- [ ] Implement `build_spatial_index()`
- [ ] Implement `get_neighbors()`
- [ ] Implement `get_k_nearest_neighbors()`
- [ ] Add lazy index building

### Phase 3.2: Boundary Detection
- [ ] Implement `get_boundary_voxels()`
- [ ] Implement `is_boundary_voxel()`
- [ ] Add boundary caching for efficiency

### Phase 3.3: Advanced Queries
- [ ] Implement `query_region()`
- [ ] Implement `get_orientation_gradient()`
- [ ] Add differentiable neighbor operations (PyTorch)

### Phase 3.4: Testing & Validation
- [ ] Unit tests for each method
- [ ] Validation against C++ neighbor finding
- [ ] Performance benchmarks
- [ ] Integration with reconstruction algorithms

## Technical Design

### Data Structures
```python
class MicFile:
    # Existing
    voxels: List[Voxel]
    positions: torch.Tensor  # (N, 3)
    orientations: torch.Tensor  # (N, 3, 3)

    # New for spatial indexing
    _kdtree: Optional[scipy.spatial.KDTree] = None
    _built_index: bool = False
    _boundary_cache: Optional[Set[int]] = None
    _neighbor_cache: Dict[Tuple[int, float], List[int]] = {}
```

### KDTree Choice: scipy vs sklearn vs PyTorch
- **scipy.spatial.KDTree**: ✅ Simple, fast, already dependency
- **sklearn.neighbors.KDTree**: More features but heavier
- **PyTorch geometric**: Differentiable but complex

**Decision**: Use `scipy.spatial.KDTree` for simplicity and performance

### Lazy Index Building
```python
@property
def kdtree(self) -> KDTree:
    """Lazy-build KDTree on first access."""
    if self._kdtree is None:
        self.build_spatial_index()
    return self._kdtree
```

### Cache Invalidation
When voxels are modified:
```python
def _invalidate_spatial_cache(self):
    """Invalidate cached spatial structures."""
    self._kdtree = None
    self._built_index = False
    self._boundary_cache = None
    self._neighbor_cache.clear()
```

## Usage Examples

### Example 1: Find Neighbors
```python
mic = MicFile.read("sample.mic")

# Find all neighbors within 0.02 m of voxel 0
neighbors = mic.get_neighbors(voxel_idx=0, radius=0.02)
print(f"Found {len(neighbors)} neighbors")

# Find 6 nearest neighbors
neighbors, distances = mic.get_k_nearest_neighbors(voxel_idx=0, k=6)
print(f"Nearest neighbor at distance {distances[0]:.6f} m")
```

### Example 2: Boundary Voxels for Reconstruction
```python
mic = MicFile.read("partial_reconstruction.mic")

# Find boundary voxels (fitted voxels next to unfitted)
boundary_indices = mic.get_boundary_voxels()
print(f"Boundary has {len(boundary_indices)} voxels")

# Get boundary voxel with highest confidence
if boundary_indices:
    best_idx = max(boundary_indices, key=lambda i: mic.voxels[i].confidence)
    print(f"Best boundary voxel: {best_idx}, confidence: {mic.voxels[best_idx].confidence:.3f}")
```

### Example 3: Orientation Gradient
```python
# Calculate misorientation with neighbors
for idx in range(len(mic.voxels)):
    gradient = mic.get_orientation_gradient(idx, k_neighbors=6)
    print(f"Voxel {idx}: max misorientation = {gradient.max().item():.2f}°")
```

## Performance Considerations

### KDTree Build Time
- **Small files (<1000 voxels)**: <1 ms
- **Large files (100K voxels)**: ~50 ms
- **Solution**: Lazy building, cache tree

### Query Time
- **Single neighbor query**: O(log N)
- **Radius query**: O(log N + M) where M = results
- **K-nearest**: O(k log N)

### Memory Usage
- **KDTree**: ~2x position array size
- **Cache**: Configurable with LRU cache

## Testing Strategy

### Unit Tests
```python
def test_build_spatial_index():
    mic = MicFile.read("test.mic")
    mic.build_spatial_index()
    assert mic._built_index
    assert mic._kdtree is not None

def test_get_neighbors_radius():
    mic = MicFile.read("test.mic")
    neighbors = mic.get_neighbors(0, radius=0.1)
    assert isinstance(neighbors, list)
    assert all(isinstance(n, int) for n in neighbors)

def test_boundary_voxels():
    mic = MicFile.read("partial.mic")
    boundary = mic.get_boundary_voxels()
    # Check that boundary voxels have both fitted and unfitted neighbors
    for idx in boundary:
        neighbors = mic.get_neighbors(idx, radius=0.02)
        phases = [mic.voxels[n].phase for n in neighbors]
        assert 0 in phases  # Has unfitted neighbor
        assert any(p > 0 for p in phases)  # Has fitted neighbor
```

### Validation Against C++
Create test harness to compare neighbor finding:
```cpp
// C++ test harness
void write_neighbors(const string& mic_file) {
    MicFile<SVoxel> mic;
    mic.Read(mic_file);

    // For each voxel, write neighbors to file
    for (int i = 0; i < mic.GetNumVoxels(); i++) {
        auto neighbors = mic.GetNeighbors(mic.GetVoxels()[i]);
        // Write to neighbors_cpp.txt
    }
}
```

```python
# Python validation
def test_neighbors_match_cpp():
    mic = MicFile.read("test.mic")
    cpp_neighbors = load_cpp_neighbors("neighbors_cpp.txt")

    for idx, expected in cpp_neighbors.items():
        actual = set(mic.get_neighbors(idx, radius=0.02))
        expected = set(expected)
        assert actual == expected, f"Mismatch for voxel {idx}"
```

## Integration with Reconstruction

### Breadth-First Reconstruction Example
```python
def breadth_first_reconstruct(mic: MicFile, seed_idx: int):
    """Simple breadth-first reconstruction propagation."""
    queue = [seed_idx]
    visited = {seed_idx}

    while queue:
        current_idx = queue.pop(0)
        current_voxel = mic.voxels[current_idx]

        # Find unfitted neighbors
        neighbors = mic.get_neighbors(current_idx, radius=2 * current_voxel.side_length)
        for n_idx in neighbors:
            if n_idx not in visited and mic.voxels[n_idx].phase == 0:
                # Fit this neighbor using current orientation as seed
                # ... fitting logic ...

                queue.append(n_idx)
                visited.add(n_idx)
```

## Benefits

1. **Performance**: KDTree queries are O(log N) vs O(N) linear search
2. **Compatibility**: Matches C++ GetNeighbors() API
3. **Flexibility**: Multiple query types (radius, k-nearest, region)
4. **Integration**: Ready for reconstruction algorithms
5. **Testability**: Each method independently testable

## Risks & Mitigations

**Risk**: KDTree rebuilding cost when voxels change
**Mitigation**: Lazy building, only rebuild when needed

**Risk**: Memory usage for large files
**Mitigation**: Optional caching, configurable cache size

**Risk**: Compatibility with C++ neighbor finding
**Mitigation**: Validation tests against C++ implementation

## Timeline

- **Phase 3.1** (Basic): 2 hours
- **Phase 3.2** (Boundary): 1 hour
- **Phase 3.3** (Advanced): 2 hours
- **Phase 3.4** (Testing): 2 hours
- **Total**: ~7 hours

## Success Criteria

- [ ] All spatial query methods implemented
- [ ] 100% test coverage for new methods
- [ ] Performance: <100ms for 100K voxel queries
- [ ] Validated against C++ neighbor finding
- [ ] Documentation complete with examples

---

**Status**: Ready to implement Phase 3.1
