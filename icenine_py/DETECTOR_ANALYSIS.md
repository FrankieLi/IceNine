# Detector Port Analysis - Answering Key Questions

## Question 1: Overlap with other Detector* C++ files

### Files Examined:
1. **Detector.h/cpp** - Core detector geometry and transformations (PRIMARY PORT TARGET)
2. **DetectorFile.h/cpp** - File I/O and parsing detector configuration
3. **DetectorData.h** - Abstract interface for detector data (IsDark/IsBright pixel queries)
4. **DetectorCalibration.h** - Detector calibration optimization (mostly empty stub)

### Analysis:

#### DetectorFile.h/cpp
**Purpose:** Parse detector configuration from text files and create `CDetectorInfo` objects
- Contains `CDetectorInfo` - a POD struct holding detector parameters
- Contains `CDetectorFile` - parses detector specification files
- `CDetectorInfo::GetDetector()` calls `CXDMDetectorFactory::MakeDetector()`

**Overlap with Detector.h/cpp:**
- **NONE** - DetectorFile is purely I/O
- The `CDetectorInfo` struct IS the parameter set we need to support
- DetectorFile.cpp line 52 has a TODO comment: "Reorganize. Manually calculating coordinate origin seems dangerous, especially when this has to happen at multiple places in the code."

**Python Port Strategy:**
- **DO NOT PORT DetectorFile** - Python has better config file handling
- **USE** `CDetectorInfo` as the template for our `DetectorParameters` dataclass
- Python can read detector configs using standard libraries (json, yaml, or ConfigParser)
- We can add a helper function `load_detector_from_dict()` if needed

#### DetectorData.h
**Purpose:** Abstract interface for querying detector pixel data
```cpp
virtual Bool IsDark( Int nJPixel, Int nKPixel ) const = 0;
virtual Bool IsBright( Int nJPixel, Int nKPixel  ) const = 0;
```

**Overlap with Detector.h/cpp:**
- **NONE** - DetectorData is about actual measured data, not geometry
- Detector.h is about geometry/transformations only
- These are separate concerns

**Python Port Strategy:**
- **DO NOT PORT DetectorData** - Wait for ImageData/Peak porting
- Python equivalent would be a protocol/ABC for detector images
- Could use numpy arrays or torch tensors directly
- This belongs with image processing, not detector geometry

#### DetectorCalibration.h
**Purpose:** Detector geometry optimization/calibration
- Currently mostly empty (stub implementation)
- Would optimize detector parameters to fit observed peaks
- Uses Eigen for optimization

**Overlap with Detector.h/cpp:**
- **MINIMAL** - Uses Detector objects, doesn't define them
- Has helper methods `GetDetParam()` and `SetDetParam()` (not implemented)
- This is a *user* of Detector, not part of its core

**Python Port Strategy:**
- **DO NOT PORT NOW** - Wait until we have full forward simulation
- Python equivalent would use PyTorch optimizer + differentiable Detector
- This will be MUCH easier in Python with autograd!
- Can leverage the differentiability we're building into Detector

### Summary - Question 1:
**NO SIGNIFICANT OVERLAP.** Detector.h/cpp is self-contained geometry/transformation code.
Other files are either:
- I/O utilities (DetectorFile) - use Python's better config handling
- Data interfaces (DetectorData) - port with image processing
- Calibration tools (DetectorCalibration) - port later with autograd

---

## Question 2: Factory Pattern Justification

### C++ Factory Analysis

The `CXDMDetectorFactory` in C++ has three methods:
1. `MakeDetector()` - Creates detector with proper coordinate origin calculation
2. `ModifyImageParameters()` - Updates detector while maintaining consistency
3. `GetImageParameters()` / `GetImageInfo()` - Extracts parameters

**Why does C++ use a factory?**

Looking at [Detector.cpp:597-607](Detector.cpp) `GetCoordOrigin()`:
```cpp
SVector3 CXDMDetectorFactory::GetCoordOrigin(...)
{
  SVector3 oCoordOrigin;
  oCoordOrigin  = - (Float) fBeamCenterJ * fPixelWidth  * oJUnitVector;
  oCoordOrigin += - (Float) fBeamCenterK * fPixelHeight * oKUnitVector;
  return oCoordOrigin;
}
```

The factory exists because:
1. **Complex initialization** - Computing `oImageCoordinateOrigin` from beam center is non-trivial
2. **Multiple constructors avoided** - C++ doesn't have great parameter validation
3. **Experiment-specific** - Comment says "This is specific to HEDM experiment"
4. **Encapsulation** - Hides coordinate system conventions from users

**Usage Analysis:**
```bash
$ grep -r "MakeDetector" Src/*.cpp
Src/DetectorFile.cpp:    CDetector oDetector =  CXDMDetectorFactory::MakeDetector(...)
```

**Only used in ONE place:** `CDetectorInfo::GetDetector()` when loading from file!

### Python Factory Justification - QUESTIONABLE

In Python, we have better alternatives:

#### Option A: **Factory Function** (simple, Pythonic)
```python
def make_detector(
    num_j_pixels: int,
    num_k_pixels: int,
    position: torch.Tensor,
    beam_center_j: float,
    beam_center_k: float,
    j_unit_vector: torch.Tensor,
    k_unit_vector: torch.Tensor,
    pixel_width: float,
    pixel_height: float,
    orientation: torch.Tensor,
    device: str = 'cpu'
) -> Detector:
    """Create detector with automatic coordinate origin calculation."""
    # Calculate image origin
    image_origin = (
        -beam_center_j * pixel_width * j_unit_vector
        - beam_center_k * pixel_height * k_unit_vector
    )

    return Detector(
        num_rows=num_k_pixels,
        num_cols=num_j_pixels,
        pixel_height=pixel_height,
        pixel_width=pixel_width,
        beam_center_j=beam_center_j,
        beam_center_k=beam_center_k,
        image_basis_j=j_unit_vector,
        image_basis_k=k_unit_vector,
        image_origin=image_origin,
        position=position,
        orientation=orientation,
        device=device
    )
```

**Pros:**
- Simple, single function
- Clearly shows coordinate origin calculation
- Standard Python pattern (like `torch.zeros()`)

**Cons:**
- None really

#### Option B: **Class Method** (also Pythonic)
```python
class Detector:
    @classmethod
    def from_beam_center(
        cls,
        num_j_pixels: int,
        num_k_pixels: int,
        position: torch.Tensor,
        beam_center_j: float,
        beam_center_k: float,
        j_unit_vector: torch.Tensor,
        k_unit_vector: torch.Tensor,
        pixel_width: float,
        pixel_height: float,
        orientation: torch.Tensor,
        device: str = 'cpu'
    ) -> 'Detector':
        """Create detector from beam center specification."""
        image_origin = cls._compute_image_origin(
            beam_center_j, beam_center_k,
            pixel_width, pixel_height,
            j_unit_vector, k_unit_vector
        )
        return cls(...)

    @staticmethod
    def _compute_image_origin(...) -> torch.Tensor:
        """Calculate image coordinate origin from beam center."""
        ...
```

**Pros:**
- Follows Python conventions (like `dict.fromkeys()`, `Path.from_uri()`)
- Encapsulates logic in the class
- Can have multiple constructors (`from_beam_center`, `from_config_dict`, etc.)

**Cons:**
- Slightly more verbose

#### Option C: **Dataclass + Builder** (over-engineered)
```python
@dataclass
class DetectorBuilder:
    """Builder pattern for complex detector initialization."""
    ...
```

**Pros:**
- Very flexible

**Cons:**
- **OVER-ENGINEERED** for this use case
- Not Pythonic
- More code to maintain

### Recommendation: **Option A + Option B hybrid**

Use BOTH:
1. **Class method `Detector.from_beam_center()`** - Primary constructor
2. **Module function `make_detector()`** - Convenience alias that calls the class method
3. **Direct `Detector()`** - For advanced users who compute image_origin themselves

```python
# detector.py

class Detector:
    def __init__(
        self,
        num_rows: int,
        num_cols: int,
        pixel_height: float,
        pixel_width: float,
        beam_center_j: float,
        beam_center_k: float,
        image_basis_j: torch.Tensor,
        image_basis_k: torch.Tensor,
        image_origin: torch.Tensor,  # <-- User provides this
        position: torch.Tensor = None,
        orientation: torch.Tensor = None,
        device: str = 'cpu',
        dtype: torch.dtype = torch.float32
    ):
        """Direct constructor - requires pre-computed image_origin."""
        ...

    @classmethod
    def from_beam_center(
        cls,
        num_j_pixels: int,
        num_k_pixels: int,
        beam_center_j: float,
        beam_center_k: float,
        pixel_width: float,
        pixel_height: float,
        j_unit_vector: torch.Tensor,
        k_unit_vector: torch.Tensor,
        position: Optional[torch.Tensor] = None,
        orientation: Optional[torch.Tensor] = None,
        device: str = 'cpu',
        dtype: torch.dtype = torch.float32
    ) -> 'Detector':
        """
        Create detector from beam center specification (recommended).

        This is the standard way to create a detector for HEDM experiments.
        The image coordinate origin is automatically calculated from the
        beam center and pixel sizes.
        """
        image_origin = cls._compute_image_origin(
            beam_center_j, beam_center_k,
            pixel_width, pixel_height,
            j_unit_vector, k_unit_vector
        )

        if position is None:
            position = torch.zeros(3, device=device, dtype=dtype)
        if orientation is None:
            orientation = torch.eye(3, device=device, dtype=dtype)

        return cls(
            num_rows=num_k_pixels,
            num_cols=num_j_pixels,
            pixel_height=pixel_height,
            pixel_width=pixel_width,
            beam_center_j=beam_center_j,
            beam_center_k=beam_center_k,
            image_basis_j=j_unit_vector,
            image_basis_k=k_unit_vector,
            image_origin=image_origin,
            position=position,
            orientation=orientation,
            device=device,
            dtype=dtype
        )

    @staticmethod
    def _compute_image_origin(
        beam_center_j: float,
        beam_center_k: float,
        pixel_width: float,
        pixel_height: float,
        j_unit_vector: torch.Tensor,
        k_unit_vector: torch.Tensor
    ) -> torch.Tensor:
        """
        Calculate image coordinate origin from beam center.

        The image origin is the point where (J=0, K=0) is located relative
        to the detector center (rotation center). This is calculated from
        the beam center position.
        """
        origin = (
            -beam_center_j * pixel_width * j_unit_vector
            - beam_center_k * pixel_height * k_unit_vector
        )
        return origin

# Convenience function (matches C++ API for easy porting)
def make_detector(*args, **kwargs) -> Detector:
    """Convenience function - calls Detector.from_beam_center()."""
    return Detector.from_beam_center(*args, **kwargs)
```

**Usage:**
```python
# Recommended: Class method
detector = Detector.from_beam_center(
    num_j_pixels=2048,
    num_k_pixels=2048,
    beam_center_j=1024.0,
    beam_center_k=1024.0,
    pixel_width=0.2,
    pixel_height=0.2,
    j_unit_vector=torch.tensor([0., 1., 0.]),
    k_unit_vector=torch.tensor([0., 0., 1.]),
    position=torch.tensor([100., 0., 0.])
)

# Alternative: Function (for C++ API compatibility)
detector = make_detector(...)

# Advanced: Direct constructor (user computes image_origin)
detector = Detector(
    num_rows=2048,
    num_cols=2048,
    image_origin=my_custom_origin,  # User calculated
    ...
)
```

**Why this is better than a Factory class:**
1. **Pythonic** - Class methods are the standard pattern
2. **Clear intent** - `from_beam_center` is self-documenting
3. **Flexible** - Can add `from_config_dict()`, `from_file()` later
4. **Simple** - Less code, no extra class
5. **Compatible** - `make_detector()` function matches C++ API

### Summary - Question 2:
**DO NOT create a Factory class.** Use class methods + convenience function instead.
This is more Pythonic and achieves the same goals with less code.

---

## Question 3: Geometry/Transformation Utilities Overlap

### Existing icenine_py Modules:

1. **constants.py** - Physical constants only (KEV_OVER_HBAR_C_IN_ANG)
2. **symmetry.py** - Crystal symmetry operations via pymatgen
3. **crystal_structure.py** - Reciprocal lattice, Miller indices
4. **diffraction_core.py** - Scattering calculations (our reference for patterns)
5. **mic_file.py** - MIC file I/O

**Analysis: No 3D geometry utilities exist yet!**

### What We Need for Detector:

1. **3D vectors** - Use `torch.Tensor` directly
2. **3x3 rotation matrices** - Use `torch.Tensor` directly
3. **Euler angles ↔ rotation matrix** - NEED TO IMPLEMENT
4. **Plane representation** - NEED TO IMPLEMENT
5. **Ray-plane intersection** - NEED TO IMPLEMENT

### Where Should Geometry Utilities Live?

#### Option A: **detector.py (internal helpers)**
```python
# detector.py

def _euler_to_rotation_matrix(phi, theta, psi):
    """ZYZ Euler angles to rotation matrix (internal)."""
    ...

def _rotation_matrix_to_euler(R):
    """Extract ZYZ Euler angles from rotation matrix (internal)."""
    ...

class Detector:
    ...
```

**Pros:**
- Self-contained
- No external dependencies
- Simple

**Cons:**
- Not reusable
- Duplicates code if Sample.py needs same utilities

#### Option B: **geometry.py (shared module)**
```python
# icenine/geometry.py

def euler_to_rotation_matrix(phi, theta, psi, convention='ZYZ'):
    """Convert Euler angles to rotation matrix."""
    ...

def rotation_matrix_to_euler(R, convention='ZYZ'):
    """Extract Euler angles from rotation matrix."""
    ...

class Plane:
    """Plane in 3D space: Ax + By + Cz + D = 0"""
    ...

class Ray:
    """Ray in 3D space: P(t) = origin + t * direction"""
    ...
```

**Pros:**
- Reusable across modules
- Clear separation of concerns
- Testable independently
- Natural place for future geometry needs (Sample, Voxel)

**Cons:**
- One more file

#### Option C: **Use scipy.spatial.transform.Rotation**
```python
from scipy.spatial.transform import Rotation as R

# Euler to matrix
rotation = R.from_euler('ZYZ', [phi, theta, psi])
matrix = rotation.as_matrix()

# Matrix to Euler
euler = R.from_matrix(matrix).as_euler('ZYZ')
```

**Pros:**
- Battle-tested
- Handles gimbal lock, convention conversions
- Already a dependency (via pymatgen)

**Cons:**
- **NOT DIFFERENTIABLE** - scipy doesn't integrate with PyTorch autograd
- Returns numpy arrays, not torch tensors
- Can't use on GPU
- **DEAL BREAKER for our use case**

#### Option D: **Use pytorch3d.transforms**
```python
from pytorch3d.transforms import euler_angles_to_matrix, matrix_to_euler_angles
```

**Pros:**
- PyTorch-based, differentiable
- GPU support
- Well-tested
- Handles multiple conventions

**Cons:**
- **Heavy dependency** - pytorch3d is large and has complex dependencies
- May not support ZYZ convention (check documentation)
- Overkill for our simple needs

### Recommendation: **Option B (geometry.py) with PyTorch**

Create `icenine/geometry.py` with:
- Euler angle ↔ rotation matrix (ZYZ convention, PyTorch, differentiable)
- Plane and Ray classes (minimal, differentiable)
- Vector operations (wrappers around torch for clarity)

**Rationale:**
1. **Differentiability required** - Rules out scipy
2. **Lightweight** - Rules out pytorch3d
3. **Reusability** - Sample.py will need rotations too
4. **Maintainability** - Clean separation of concerns
5. **Following diffraction_core pattern** - They implemented their own physics

**What to include in geometry.py:**
```python
# icenine/geometry.py
"""
Geometric primitives and transformations for X-ray diffraction.

PyTorch-based implementations for differentiability and GPU support.
"""

import torch
from typing import Tuple, Optional
from dataclasses import dataclass

# Euler angle conversions (ZYZ convention for IceNine)
def euler_to_rotation_matrix(
    phi: torch.Tensor,
    theta: torch.Tensor,
    psi: torch.Tensor
) -> torch.Tensor:
    """Convert ZYZ Euler angles to rotation matrix (batched, differentiable)."""
    ...

def rotation_matrix_to_euler(R: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Extract ZYZ Euler angles from rotation matrix (batched, differentiable)."""
    ...

# Plane representation
@dataclass
class Plane:
    """Plane: A*x + B*y + C*z + D = 0"""
    coeffs: torch.Tensor  # (4,) for single plane, (N, 4) for batch

    @property
    def normal(self) -> torch.Tensor:
        return self.coeffs[..., :3]

    def normalize(self) -> 'Plane':
        """Normalize plane equation."""
        ...

# Ray representation
@dataclass
class Ray:
    """Ray: P(t) = origin + t * direction"""
    origin: torch.Tensor     # (3,) or (N, 3)
    direction: torch.Tensor  # (3,) or (N, 3)

    def intersect_plane(self, plane: Plane) -> Tuple[torch.Tensor, torch.Tensor]:
        """Compute ray-plane intersection (batched, differentiable)."""
        ...

# Vector operations (for clarity, wraps torch)
def normalize_vector(v: torch.Tensor, dim: int = -1, eps: float = 1e-8) -> torch.Tensor:
    """Normalize vectors along dimension (batched)."""
    ...

def cross_product(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    """Cross product (batched)."""
    ...
```

**What NOT to include:**
- BBox2D - Not needed for detector geometry
- Matrix inverse, determinant - Use torch.linalg directly
- Quaternions - Not used in Detector (only in Sample/Orientation)
- Projections, reflections - Not needed yet

### Summary - Question 3:
**CREATE icenine/geometry.py** with minimal, differentiable PyTorch implementations.
- NO overlap with existing code
- Euler angles (ZYZ) + Plane + Ray only
- Reusable for future ports (Sample, Voxel, etc.)
- Keep it MINIMAL - add more only when needed

---

## Updated Recommendations

### Module Structure (REVISED):
```
icenine_py/icenine/
├── detector.py          # Main Detector class (NEW)
│   ├── class Detector
│   │   ├── __init__()                      # Direct constructor
│   │   ├── from_beam_center() classmethod  # Recommended constructor
│   │   ├── lab_to_pixel()
│   │   ├── pixel_to_lab_coordinate()
│   │   ├── set_orientation()
│   │   └── ... (all coordinate transforms)
│   └── make_detector() function           # Convenience wrapper
├── geometry.py          # Geometric primitives (NEW)
│   ├── euler_to_rotation_matrix()
│   ├── rotation_matrix_to_euler()
│   ├── class Plane
│   ├── class Ray
│   └── normalize_vector(), cross_product()
└── (no transforms.py - utilities go in geometry.py or detector.py)
```

### What NOT to Port:
1. ❌ **DetectorFile** - Use Python config parsing
2. ❌ **DetectorData** - Port with image processing later
3. ❌ **DetectorCalibration** - Port after full forward simulation
4. ❌ **CXDMDetectorFactory class** - Use class methods instead

### Implementation Priority (REVISED):
1. **Phase 1:** geometry.py (4-5 hours)
   - Euler angles (ZYZ)
   - Plane class
   - Ray class

2. **Phase 2:** detector.py (8-10 hours)
   - Detector class with all transformations
   - `from_beam_center()` class method
   - `make_detector()` convenience function

3. **Phase 3:** Testing (6-8 hours)
   - Unit tests for geometry primitives
   - Unit tests for detector transformations
   - Integration test comparing with C++ output
   - Differentiability tests

4. **Phase 4:** Documentation (3-4 hours)
   - Docstrings
   - Examples
   - Coordinate system conventions

**Total: 21-27 hours** (reduced from 25-36 hours by removing factory class and simplifying)

### Key Design Decisions (FINAL):
1. ✅ Class method `Detector.from_beam_center()` instead of factory class
2. ✅ Create `geometry.py` for reusable PyTorch-based primitives
3. ✅ Use torch.Tensor directly for vectors/matrices (no custom classes)
4. ✅ ZYZ Euler convention (matches C++)
5. ✅ Differentiable throughout
6. ✅ Support batching where beneficial
7. ✅ Device-agnostic (CPU/GPU)
