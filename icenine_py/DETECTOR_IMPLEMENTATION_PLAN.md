# Detector Implementation Plan - Final Version

**Status:** Ready to implement
**Approach:** Refactor geometry utilities during Detector port
**Timeline:** 24-32 hours

See also:
- [DETECTOR_PORT_PLAN.md](DETECTOR_PORT_PLAN.md) - Detailed technical spec
- [DETECTOR_ANALYSIS.md](DETECTOR_ANALYSIS.md) - Analysis of related C++ files
- [DETECTOR_GEOMETRY_REUSE.md](DETECTOR_GEOMETRY_REUSE.md) - Existing code analysis

## Implementation Order

### Phase 1: Create geometry.py (4-5 hours)

**Goal:** Centralized geometric primitives for IceNine Python

#### Step 1.1: Move Euler angle functions from mic_file.py (1 hour)

**Actions:**
1. Create `icenine_py/icenine/geometry.py`
2. Move these functions from mic_file.py:
   - `euler_to_matrix()` (NumPy/scipy version)
   - `matrix_to_euler()` (NumPy/scipy version)
   - `euler_to_matrix_torch()` (PyTorch version) ✅ ALREADY TESTED
3. Enhance docstrings with:
   - C++ reference (3dMath.cpp line numbers)
   - Convention clarification (ZXZ intrinsic = Bunge)
   - Batching behavior
   - Differentiability notes

**Success criteria:**
- [ ] geometry.py created with moved functions
- [ ] Docstrings reference C++ code
- [ ] No functionality changes (exact copy)

#### Step 1.2: Implement matrix_to_euler_torch() (2-3 hours)

**Goal:** Differentiable PyTorch version of matrix → Euler angles

**C++ Reference:** XDM++/libXDM/3dMath.cpp lines 180-209 (GetEulerAngles)

**Implementation:**
```python
def matrix_to_euler_torch(
    rotation_matrix: torch.Tensor,
    epsilon: float = 1e-6
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Convert rotation matrix to Bunge Euler angles (degrees) - PyTorch.

    Differentiable version for gradient-based optimization.

    C++ reference: XDM++/libXDM/3dMath.cpp GetEulerAngles (lines 180-209)

    Args:
        rotation_matrix: Shape (3, 3) or (N, 3, 3)
        epsilon: Threshold for gimbal lock detection

    Returns:
        (phi1_deg, Phi_deg, phi2_deg) in degrees

    Note:
        Handles gimbal lock cases when Phi ≈ 0 or Phi ≈ π
    """
    # Extract matrix elements
    # Handle batching (both scalar and batched inputs)

    # Case 1: Phi ≈ 0 (m[2][2] > 0.999999)
    #   phi1 = 0, Phi = 0, phi2 = atan2(m[1][0], m[0][0])

    # Case 2: Phi ≈ π (m[2][2] < -0.999999)
    #   phi1 = 0, Phi = π, phi2 = atan2(m[0][1], m[0][0])

    # Case 3: General case
    #   phi1 = atan2(m[0][2], -m[1][2])
    #   Phi = atan2(sqrt(m[2][0]² + m[2][1]²), m[2][2])
    #   phi2 = atan2(m[2][0], m[2][1])

    # Use torch.where for differentiability (not if/else)
    # Normalize angles to [0, 360] for phi1/phi2, [0, 180] for Phi
    # Convert radians to degrees
```

**Key challenges:**
- Gimbal lock handling with differentiable torch.where()
- Batching support (both scalar and batched)
- Numerical stability near singularities
- Match C++ behavior exactly

**Testing:**
1. Identity matrix → (0, 0, 0)
2. Known rotations (90° around Z, etc.)
3. Round-trip: euler → matrix → euler
4. Compare with scipy results
5. Gradient flow test
6. Batch processing test

**Success criteria:**
- [ ] Handles scalar and batched inputs
- [ ] Gimbal lock cases correct
- [ ] Matches scipy.Rotation output
- [ ] Round-trip error < 1e-5 degrees
- [ ] Gradients flow correctly

#### Step 1.3: Add Plane and Ray classes (1 hour)

**Goal:** Minimal geometric primitives for ray-detector intersection

```python
@dataclass
class Plane:
    """
    Plane in 3D space: Ax + By + Cz + D = 0

    Attributes:
        coeffs: Plane coefficients [A, B, C, D]
                Shape (4,) for single plane
                Shape (N, 4) for batched planes
    """
    coeffs: torch.Tensor

    @property
    def normal(self) -> torch.Tensor:
        """Return normal vector [A, B, C]."""
        return self.coeffs[..., :3]

    @property
    def d(self) -> torch.Tensor:
        """Return offset D."""
        return self.coeffs[..., 3]

    def normalize(self) -> 'Plane':
        """Normalize so |normal| = 1."""
        norm = torch.norm(self.normal, dim=-1, keepdim=True)
        return Plane(self.coeffs / norm)

    def distance_to_point(self, point: torch.Tensor) -> torch.Tensor:
        """Signed distance from point to plane."""
        # distance = (A*x + B*y + C*z + D) / sqrt(A² + B² + C²)
        ...

@dataclass
class Ray:
    """
    Ray in 3D space: P(t) = origin + t * direction

    Attributes:
        origin: Ray origin, shape (3,) or (N, 3)
        direction: Ray direction (should be normalized), shape (3,) or (N, 3)
    """
    origin: torch.Tensor
    direction: torch.Tensor

    def intersect_plane(
        self, plane: Plane
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Compute ray-plane intersection.

        Returns:
            intersects: Boolean tensor
            t: Intersection parameter (valid only if intersects=True)

        Math:
            Ray: P(t) = O + t*D
            Plane: N·P + d = 0
            Intersection: t = -(N·O + d) / (N·D)
            Valid if N·D != 0
        """
        # denominator = dot(plane.normal, ray.direction)
        # numerator = -(dot(plane.normal, ray.origin) + plane.d)
        # intersects = abs(denominator) > epsilon
        # t = numerator / denominator (where intersects)
        ...
```

**Success criteria:**
- [ ] Plane normalization works
- [ ] Ray-plane intersection correct
- [ ] Handles parallel rays (no intersection)
- [ ] Batching works
- [ ] Differentiable

#### Step 1.4: Update mic_file.py (30 min)

**Goal:** Maintain backward compatibility while using geometry.py

```python
# mic_file.py - add at top
from .geometry import (
    euler_to_matrix,
    matrix_to_euler,
    euler_to_matrix_torch,
    # matrix_to_euler_torch,  # Add when implemented
)

# Keep in __all__ for backward compatibility
__all__ = [
    'MicFile',
    'MicReader',
    'MicWriter',
    'euler_to_matrix',
    'matrix_to_euler',
    'euler_to_matrix_torch',
]

# Remove the old implementations (lines 670-814)
# Everything else stays the same
```

**Testing:**
- [ ] Run existing mic_file tests: `pytest tests/test_mic_file.py`
- [ ] Verify no regressions
- [ ] Import backward compatibility: `from icenine.mic_file import euler_to_matrix`

---

### Phase 2: Implement Detector class (8-10 hours)

**Goal:** Core detector geometry with all coordinate transformations

#### Step 2.1: Detector class structure (2 hours)

**File:** `icenine_py/icenine/detector.py`

**Imports:**
```python
import torch
from typing import Optional, Tuple
from dataclasses import dataclass

from .geometry import euler_to_matrix_torch, matrix_to_euler_torch, Plane, Ray
```

**Class skeleton:**
```python
class Detector:
    """
    X-ray detector geometry and coordinate transformations.

    Coordinate Systems:
    - Lab frame: 3D global coordinates (beam along +X by convention)
    - Detector frame: 2D (J, K) coordinates on detector surface (mm)
    - Pixel coordinates: 2D (col, row) discrete indices

    Conventions:
    - J direction: horizontal (columns)
    - K direction: vertical (rows)
    - Beam center: where direct beam hits detector
    - Rotation center: detector position in lab frame
    """

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
        image_origin: torch.Tensor,
        position: Optional[torch.Tensor] = None,
        orientation: Optional[torch.Tensor] = None,
        device: str = 'cpu',
        dtype: torch.dtype = torch.float32
    ):
        """Direct constructor - requires pre-computed image_origin."""
        self.device = device
        self.dtype = dtype

        # Detector dimensions
        self.num_rows = num_rows
        self.num_cols = num_cols
        self.pixel_height = pixel_height
        self.pixel_width = pixel_width
        self.pixel_half_height = pixel_height / 2.0
        self.pixel_half_width = pixel_width / 2.0

        # Beam center (in pixels)
        self.beam_center_j = beam_center_j
        self.beam_center_k = beam_center_k

        # Position and orientation
        self._position = position if position is not None else torch.zeros(3, device=device, dtype=dtype)
        self._orientation = orientation if orientation is not None else torch.eye(3, device=device, dtype=dtype)

        # Image coordinate system (detector frame - constant)
        self._det_frame_image_basis_j = image_basis_j.to(device=device, dtype=dtype)
        self._det_frame_image_basis_k = image_basis_k.to(device=device, dtype=dtype)
        self._det_frame_image_origin = image_origin.to(device=device, dtype=dtype)

        # Lab frame basis (updated by orientation)
        self._lab_frame_image_basis_j = None
        self._lab_frame_image_basis_k = None
        self._lab_frame_image_origin = None

        # Detector plane (updated by position/orientation)
        self._lab_frame_detector_plane = None

        # Update derived quantities
        self._update_lab_frame_geometry()
```

**Key internal methods:**
```python
def _update_lab_frame_geometry(self):
    """Update lab frame quantities after orientation/position change."""
    # Rotate basis vectors
    self._lab_frame_image_basis_j = self._orientation @ self._det_frame_image_basis_j
    self._lab_frame_image_basis_k = self._orientation @ self._det_frame_image_basis_k
    self._lab_frame_image_origin = self._orientation @ self._det_frame_image_origin

    # Calculate detector plane
    self._calculate_detector_plane()

def _calculate_detector_plane(self):
    """Calculate detector plane in lab frame."""
    # Define 3 points on detector plane
    p1 = self._orientation @ torch.tensor([1., 0., 0.], device=self.device, dtype=self.dtype)
    p2 = self._orientation @ torch.tensor([0., 1., 0.], device=self.device, dtype=self.dtype)
    p3 = self._orientation @ torch.tensor([0., 0., 0.], device=self.device, dtype=self.dtype)

    # Translate to detector position
    p1 = p1 + self._position
    p2 = p2 + self._position
    p3 = p3 + self._position

    # Compute plane normal via cross product
    edge1 = p3 - p1
    edge2 = p2 - p1
    normal = torch.cross(edge2, edge1)
    normal = normal / torch.norm(normal)

    # Plane equation: N·P + d = 0
    d = -torch.dot(normal, p1)

    coeffs = torch.cat([normal, d.unsqueeze(0)])
    self._lab_frame_detector_plane = Plane(coeffs)
```

**Success criteria:**
- [ ] Initialization works
- [ ] Lab frame geometry updated correctly
- [ ] Detector plane calculated correctly

#### Step 2.2: Coordinate transformations (3 hours)

Implement all transformation methods (see DETECTOR_PORT_PLAN.md Phase 2.3):
- `lab_to_detector_coordinate()`
- `detector_to_lab_coordinate()`
- `lab_to_pixel()`
- `pixel_to_lab_coordinate()`

**Success criteria:**
- [ ] All transformations implemented
- [ ] Round-trip transformations work (lab → pixel → lab)
- [ ] Batching supported

#### Step 2.3: Detector transformations (2 hours)

Implement positioning methods:
- `set_orientation()`, `set_orientation_euler()`
- `set_position()`, `translate()`, `rotate()`
- `intersect_ray()`

**Success criteria:**
- [ ] Transformations update internal state correctly
- [ ] Intersect ray works

#### Step 2.4: Properties and utilities (1 hour)

Implement properties and helper methods.

---

### Phase 3: Construction methods (2-3 hours)

Implement:
- `Detector.from_beam_center()` class method
- `make_detector()` convenience function
- `DetectorParameters` dataclass
- `get_detector_parameters()` function

**Success criteria:**
- [ ] from_beam_center() creates valid detector
- [ ] Image origin calculated correctly
- [ ] Parameters can be extracted and serialized

---

### Phase 4: Testing (6-8 hours)

#### Unit tests (4-5 hours)

**File:** `tests/test_detector.py`

```python
class TestDetectorGeometry:
    def test_initialization()
    def test_from_beam_center()
    def test_coordinate_transformations()
    def test_round_trip_transforms()

class TestDetectorTransformations:
    def test_translation()
    def test_rotation()
    def test_set_orientation()

class TestRayIntersection:
    def test_perpendicular_ray()
    def test_parallel_ray()
    def test_angled_ray()

class TestBatchedOperations:
    def test_batched_lab_to_pixel()
    def test_batched_pixel_to_lab()

class TestDifferentiability:
    def test_gradient_through_coordinate_transform()
    def test_gradient_through_orientation()
```

**File:** `tests/test_geometry.py`

```python
class TestEulerConversions:
    def test_euler_to_matrix_torch()
    def test_matrix_to_euler_torch()
    def test_round_trip()
    def test_gimbal_lock_cases()
    def test_batching()
    def test_gradients()

class TestPlaneRay:
    def test_plane_normalization()
    def test_ray_plane_intersection()
    def test_parallel_ray()
```

#### Integration tests (2-3 hours)

**File:** `tests/test_detector_cpp_compatibility.py`

Compare with C++ implementation:
1. Load detector config from C++ config file
2. Create Python detector
3. Test transformations match C++ output

**Success criteria:**
- [ ] All unit tests pass
- [ ] Integration test matches C++ within 1e-5
- [ ] Gradients flow correctly
- [ ] GPU tests pass (if CUDA available)

---

### Phase 5: Documentation (3-4 hours)

1. **Module docstrings** (1 hour)
   - geometry.py: Coordinate conventions, C++ references
   - detector.py: Usage guide, coordinate systems

2. **Examples** (2 hours)
   - `examples/detector_usage.py`: Basic usage
   - `examples/detector_calibration.py`: Differentiable optimization

3. **Update __init__.py** (0.5 hour)
   ```python
   # icenine/__init__.py
   from . import geometry
   from . import detector

   __all__ = [..., "geometry", "detector"]
   ```

4. **README updates** (0.5 hour)

---

## Testing Strategy

### During Development
- Write tests alongside implementation
- Run `pytest tests/test_geometry.py` after Phase 1
- Run `pytest tests/test_detector.py` after Phase 2

### Final Validation
```bash
# All tests
pytest icenine_py/tests/

# Specific modules
pytest icenine_py/tests/test_geometry.py -v
pytest icenine_py/tests/test_detector.py -v

# Integration
pytest icenine_py/tests/test_detector_cpp_compatibility.py -v

# Coverage
pytest --cov=icenine --cov-report=html
```

---

## Success Criteria (Final Checklist)

### Phase 1 - geometry.py
- [ ] Euler functions moved from mic_file.py
- [ ] matrix_to_euler_torch() implemented and tested
- [ ] Plane and Ray classes implemented
- [ ] mic_file.py updated, tests pass
- [ ] No breaking changes to existing code

### Phase 2 - detector.py
- [ ] Detector class implemented
- [ ] All coordinate transformations work
- [ ] Batching supported
- [ ] Differentiable

### Phase 3 - Construction
- [ ] from_beam_center() works correctly
- [ ] Parameters can be extracted

### Phase 4 - Testing
- [ ] All unit tests pass
- [ ] Integration test matches C++ (error < 1e-5)
- [ ] GPU tests pass
- [ ] Gradient tests pass

### Phase 5 - Documentation
- [ ] All docstrings complete
- [ ] Examples run without errors
- [ ] README updated

---

## Risk Mitigation

### Potential Issues

1. **Breaking mic_file.py during refactor**
   - Mitigation: Run tests immediately after changes
   - Fallback: Keep backup of original code

2. **matrix_to_euler_torch() gimbal lock handling**
   - Mitigation: Test edge cases extensively
   - Reference: Compare with scipy implementation

3. **Coordinate convention confusion**
   - Mitigation: Document clearly, add diagrams
   - Validation: Test against C++ output

4. **Device/dtype inconsistencies**
   - Mitigation: Strict type checking in __init__
   - Testing: Test CPU and GPU explicitly

---

## Next Steps

**Ready to begin implementation!**

1. Create feature branch: `git checkout -b feature/detector-port`
2. Start with Phase 1.1 (move Euler functions)
3. Commit frequently with descriptive messages
4. Run tests after each phase

**First commit should be:**
- Create geometry.py
- Move euler_to_matrix, matrix_to_euler, euler_to_matrix_torch from mic_file.py
- Update mic_file.py imports
- Verify mic_file tests still pass

Let's do this! 🚀
