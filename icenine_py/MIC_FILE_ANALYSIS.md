# MicFile Migration Analysis - Progress Report

**Date**: 2025-11-11
**Status**: Architecture analysis completed, ready to design Python implementation

## Summary

Completed comprehensive analysis of the C++ MicFile system in preparation for Python/PyTorch migration. The system is a specialized file format and data structure for storing 3D microstructure data with crystal orientations.

## Files Analyzed

### Core Components (Fully Analyzed)

1. **XDM++/libXDM/MicIO.h** (806 lines)
   - Template-based file I/O system
   - Two main implementations: `MicFile<SVoxel>` (triangular) and `MicFile<SquareVoxel>` (square grid)
   - Handles Euler angle conversion (degrees in file ↔ radians in memory)

2. **XDM++/libXDM/Voxel.h** (181 lines)
   - `CEverythingVoxel` (SVoxel): Full voxel with orientation matrix, deformation, confidence
   - `SquareVoxel`: Simplified version with 4 vertices

3. **XDM++/libXDM/MicGrid.h** (910 lines)
   - High-performance 2D grid representation using boost::multi_array
   - `CMicGrid`: Triangular mesh grid with adaptive refinement
   - `RectMicGrid<ShapeT>`: Rectilinear (square) grid template
   - Spatial lookup and neighbor finding

4. **XDM++/libXDM/MicMesh.h** (456 lines)
   - Quadtree-based spatial indexing (`CMicMesh`)
   - 3D layered structure (`CGeneralShapeLocator`)
   - Overlap detection and neighbor queries

### Example Files Found

- 37 `.mic` files in `DataFiles.back/`
- 1 active file: `DataFiles/Au1007_small.mic`
- File formats confirmed (triangular and square grid)

## Key Findings

### 1. File Format Structure

**Triangular Mesh Format** (SVoxel):
```
Line 1: <SideLength>
Line 2+: <x> <y> <z> <gen> <phase> <phi1°> <Phi°> <phi2°> <confidence> <cost> <overlap> [<deformation 6×>]
```

Example from Au1007_small.mic:
```
  1.2000000E-02
 -1.2000000E-02  0.0000000E+00  0.0000000E+00           1           3           1   355.4292       5.186272       29.31929      0.1034483
```

**Square Grid Format** (SquareVoxel):
```
Line 1: <SampleSideLength> <VoxelSideLength> <OriginX> <OriginY> <OriginZ>
Line 2+: <center_x> <center_y> <center_z> <phi1°> <Phi°> <phi2°> <confidence> <phase>
```

### 2. Critical Design Details

**Orientation Representation**:
- **File storage**: Bunge Euler angles (φ₁, Φ, φ₂) in **degrees**
- **In-memory**: 3×3 rotation matrix (active Euler convention)
- **Conversion happens at I/O boundary**: `DEGREE_TO_RADIAN()` on read, `RadianToDegree()` on write

**Voxel Properties** (SVoxel):
- Position: (x, y, z) coordinates
- Orientation: 3×3 matrix `oOrientMatrix`
- Deformation: 3×3 matrix `oDeformation` (metric modifier for strain)
- Quality metrics: `fConfidence`, `fCost`, `fPixelOverlapRatio`
- Metadata: `nPhase`, `nGeneration`, `nID`

**Grid Organization**:
- 2D slices organized by z-height in 3D
- Triangular mesh supports adaptive refinement (generation levels)
- Square grid is uniform resolution
- Both support spatial queries (neighbors, overlap detection)

### 3. Data Structure Hierarchy

```
MicFile<VoxelType>                    (File I/O)
    ↓
RectMicGrid<VoxelType>                (2D grid, spatial indexing)
    ↓
CGeneralShapeLocator<ShapeContainer>  (3D layered structure)
    ↓
Individual Voxels (SVoxel/SquareVoxel)
```

### 4. Key Operations

**File I/O**:
- Read/Write text format with scientific notation
- Automatic angle unit conversion
- Support for serialization (Boost::serialization)

**Spatial Queries**:
- Point location: Find voxel at (x,y,z)
- Neighbor finding: Get adjacent voxels
- Overlap detection: Find voxels intersecting a bounding box
- Range search: Quadtree-based spatial indexing

**Grid Operations**:
- Insert voxel into spatial structure
- Refine/coarsen adaptive mesh
- Calculate Nye tensor field (for dislocation density)

## Simplification Opportunities for Python Port

### 1. Consolidate Formats
- **Current**: Two separate implementations (triangular vs square)
- **Proposed**: Single unified PyTorch-based representation
- **Benefit**: Reduced code complexity, easier maintenance

### 2. Use PyTorch Tensors Throughout
- **Orientations**: (N, 3, 3) rotation matrices (differentiable)
- **Positions**: (N, 3) coordinates
- **Properties**: (N,) for confidence, cost, etc.
- **Benefit**: GPU acceleration, automatic differentiation, batched operations

### 3. Simplify Spatial Indexing
- **Current**: Custom quadtree + boost::multi_array
- **Proposed**: Use scipy.spatial.KDTree or PyTorch geometric libraries
- **Benefit**: Leverage well-tested libraries, simpler code

### 4. Modern File Format
- **Keep**: Text-based .mic format for C++ compatibility
- **Add**: PyTorch `.pt` format for fast native serialization
- **Benefit**: Faster I/O, preserves gradients, easier Python integration

### 5. Remove Unnecessary Features
- **Skip**: Adaptive refinement (generation levels) - rare in practice
- **Skip**: Deformation tensor - often unused
- **Keep**: Core properties (position, orientation, confidence, phase)
- **Benefit**: Simpler implementation, focus on common use cases

## Migration Requirements (from User)

### 1. Regression Testing Against C++
**Strategy**:
- Generate ground truth .mic files from C++ code
- Read with Python, write back out, compare
- Validate: positions, orientations (angle conversion!), properties
- Test both triangular and square formats

**Test Cases Needed**:
- Empty file
- Single voxel
- Small grid (10×10)
- Multiple phases
- High/low confidence voxels
- Euler angle edge cases (0°, 90°, 180°, 360°)

### 2. Automatic Differentiation
**Differentiable Operations**:
- ✅ Orientation transformations (rotation matrices)
- ✅ Crystal structure calculations (G-vectors from orientations)
- ✅ Forward diffraction simulation
- ❌ File I/O (not differentiable, but doesn't need to be)
- ❌ Spatial indexing (discrete operation)

**Use Case**: Gradient-based microstructure optimization
- Example: Optimize voxel orientations to match experimental data
- Requires: `orientation_matrix.requires_grad = True`

### 3. Serialization Compatibility
**Requirements**:
- Must read existing C++ .mic files (backward compatibility)
- Must write .mic files readable by C++ code (forward compatibility)
- Should support PyTorch native format for performance

**Critical Details**:
- Preserve Euler angle convention (Bunge, active)
- Preserve degree/radian conversion at I/O boundary
- Maintain file format exactly (whitespace, precision, column order)

## Proposed Python Architecture

### Module Structure
```
icenine_py/icenine/
├── mic_file.py           # Core MicFile class
├── voxel.py              # Voxel data structures
├── mic_grid.py           # Spatial grid (optional, if needed)
└── mic_io.py             # File I/O utilities
```

### Core Classes

**1. Voxel (dataclass)**
```python
@dataclass
class Voxel:
    position: torch.Tensor          # (3,) - xyz coordinates
    orientation: torch.Tensor       # (3, 3) - rotation matrix
    confidence: float = 0.0
    cost: float = 0.0
    phase: int = 1
    id: int = -1

    # Optional advanced properties
    deformation: Optional[torch.Tensor] = None  # (3, 3)
    overlap_ratio: float = 0.0
```

**2. MicFile (main interface)**
```python
class MicFile:
    def __init__(self, voxels: List[Voxel], metadata: dict):
        # Store as PyTorch tensors for batched operations
        self.positions = torch.stack([v.position for v in voxels])      # (N, 3)
        self.orientations = torch.stack([v.orientation for v in voxels]) # (N, 3, 3)
        self.confidence = torch.tensor([v.confidence for v in voxels])   # (N,)
        # ... other properties

    @classmethod
    def read(cls, filename: str) -> 'MicFile':
        """Read .mic file (C++ compatible)"""
        # Parse text format
        # Convert Euler angles (degrees) → rotation matrices

    def write(self, filename: str):
        """Write .mic file (C++ compatible)"""
        # Convert rotation matrices → Euler angles (degrees)
        # Write text format with scientific notation

    def save_torch(self, filename: str):
        """Save PyTorch native format (fast, preserves gradients)"""
        torch.save({'positions': self.positions,
                    'orientations': self.orientations, ...}, filename)

    @classmethod
    def load_torch(cls, filename: str) -> 'MicFile':
        """Load PyTorch native format"""

    def get_neighbors(self, voxel_idx: int, radius: float) -> List[int]:
        """Find neighboring voxels (uses KDTree)"""
```

### File I/O Strategy

**Reading**:
1. Parse text file line by line
2. Extract Euler angles (φ₁, Φ, φ₂) in degrees
3. Convert to radians: `angles_rad = angles_deg * π/180`
4. Build rotation matrix: `R = euler_to_matrix(phi1, Phi, phi2, convention='ZXZ')`
5. Store as PyTorch tensor

**Writing**:
1. Extract rotation matrices from PyTorch tensors
2. Convert to Euler angles: `phi1, Phi, phi2 = matrix_to_euler(R, convention='ZXZ')`
3. Convert to degrees: `angles_deg = angles_rad * 180/π`
4. Format with scientific notation: `f"{value:13.7E}"`
5. Write text file with exact column spacing

## Testing Strategy

### Phase 1: Basic I/O
- Test 1: Read empty.mic → verify empty
- Test 2: Read oneTriangle.mic → verify 1 voxel, correct position/orientation
- Test 3: Read Au1007_small.mic → verify all voxels parsed
- Test 4: Write → Read → Compare (round-trip test)

### Phase 2: Euler Angle Conversion
- Test 5: Identity rotation (0, 0, 0)° → Identity matrix
- Test 6: 90° rotations → Verify matrix elements
- Test 7: Random rotations → Round-trip (matrix → Euler → matrix)
- Test 8: Gimbal lock cases (Φ = 0°, 180°)

### Phase 3: C++ Compatibility
- Test 9: Python write → C++ read → Verify in C++ code
- Test 10: C++ write → Python read → Verify exact match
- Test 11: Large file (500+ voxels) → Performance check

### Phase 4: Differentiation
- Test 12: Orientations with requires_grad=True
- Test 13: Compute G-vectors from orientations → Backward pass
- Test 14: Gradient of simulation loss w.r.t. orientations

## Next Steps

### Immediate (Ready to Implement)
1. Create `mic_file.py` with basic data structures
2. Implement Euler angle ↔ rotation matrix conversion
3. Implement .mic file reader
4. Create C++ test harness to generate ground truth files

### Short Term
5. Implement .mic file writer
6. Add round-trip tests
7. Validate against C++ on all 37 example files
8. Add PyTorch native serialization

### Medium Term
9. Add spatial indexing (KDTree-based neighbor queries)
10. Implement differentiable forward simulator integration
11. Create visualization utilities
12. Performance optimization (GPU, batching)

## Open Questions

1. **Deformation Tensor**: Is this used in practice? If not, can we skip it?
2. **Adaptive Refinement**: The "generation" field supports multi-resolution. Still needed?
3. **Triangular vs Square**: Which format is more common? Should we prioritize one?
4. **Spatial Indexing**: Is neighbor finding performance-critical? Do we need custom implementation?

## Risk Assessment

**Low Risk**:
- ✅ File format is well-documented in code
- ✅ Euler angle conversion is standard (scipy, pytorch3d available)
- ✅ Simple text format, easy to parse/validate

**Medium Risk**:
- ⚠️ Angle convention subtleties (active vs passive, ZXZ vs ZYZ)
- ⚠️ Floating point precision in angle conversion
- ⚠️ Performance of spatial queries for large datasets

**Mitigation**:
- Extensive testing with C++ ground truth
- Use established libraries (scipy.spatial.transform)
- Profile and optimize hot paths

## Estimated Effort

- **Basic I/O**: 1-2 days (read/write, Euler conversion)
- **Testing**: 1 day (C++ harness, regression tests)
- **Spatial indexing**: 1 day (KDTree integration)
- **Differentiation integration**: 0.5 days (already have diffraction_core)
- **Documentation**: 0.5 days

**Total**: ~4-5 days for complete, tested implementation

---

**Status**: Ready to proceed with implementation. Architecture is clear, risks are manageable, testing strategy is defined.
