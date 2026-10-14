# Detector.h/cpp Port to Python - Implementation Plan (REVISED)

## Overview
Port the C++ Detector class to Python with PyTorch for differentiable detector geometry and coordinate transformations. The detector is a critical component that handles the geometry of the X-ray detector including its position, orientation, and coordinate transformations between lab frame and pixel coordinates.

**IMPORTANT:** See [DETECTOR_ANALYSIS.md](DETECTOR_ANALYSIS.md) for detailed analysis of related files, factory pattern justification, and geometry utilities overlap.

## Goals
1. **Idiomatic Python**: Use Python best practices (class methods, type hints, properties)
2. **PyTorch-based**: All math operations in PyTorch for GPU acceleration and differentiability
3. **Differentiable**: Enable gradient flow through detector geometry for optimization
4. **Batched operations**: Support batch processing where beneficial
5. **Minimal dependencies**: Reuse PyTorch, avoid heavy dependencies like pytorch3d

## Architecture

### Module Structure (REVISED)
```
icenine_py/icenine/
├── detector.py          # Main Detector class (NEW)
│   ├── class Detector with class methods (NO factory class)
│   └── make_detector() convenience function
└── geometry.py          # Geometric primitives (NEW)
    ├── Euler angle conversions (ZYZ convention)
    ├── class Plane
    └── class Ray
```

**Note:** No `transforms.py` - utilities live in `geometry.py` or `detector.py` directly.

## Phase 1: Geometric Primitives (geometry.py)

### 1.1 Vector and Matrix Operations
Since we're using PyTorch, we can use `torch.Tensor` directly for vectors and matrices.

**Key operations needed:**
- 3D vector operations (cross product, dot product, normalization)
- 3x3 matrix operations (matrix multiplication, Euler angle conversions)
- Euler angle to rotation matrix (Active/Passive conventions)
- Matrix to Euler angles extraction

**Implementation approach:**
```python
# Use torch.Tensor directly for vectors (shape: (3,) or (N, 3))
# Use torch.Tensor for 3x3 matrices (shape: (3, 3) or (N, 3, 3))

def euler_to_rotation_matrix(phi, theta, psi, convention='ZYZ'):
    """Convert Euler angles to rotation matrix (batched)"""

def rotation_matrix_to_euler(R, convention='ZYZ'):
    """Extract Euler angles from rotation matrix (batched)"""

def normalize_vector(v, dim=-1):
    """Normalize vectors (batched)"""
```

### 1.2 Plane Class
**Purpose:** Represent a plane in 3D space (Ax + By + Cz + D = 0)

```python
@dataclass
class Plane:
    """
    Plane representation: A*x + B*y + C*z + D = 0

    Attributes:
        coeffs: Plane coefficients [A, B, C, D], shape (4,) or (N, 4)
    """
    coeffs: torch.Tensor  # (4,) or (N, 4) for batched planes

    @property
    def normal(self) -> torch.Tensor:
        """Return plane normal vector (first 3 components)"""

    def normalize(self) -> 'Plane':
        """Normalize plane equation so |normal| = 1"""

    def distance_to_point(self, point: torch.Tensor) -> torch.Tensor:
        """Calculate signed distance from point to plane"""
```

### 1.3 Ray Class
**Purpose:** Represent a ray for intersection testing

```python
@dataclass
class Ray:
    """
    Ray representation: P(t) = origin + t * direction

    Attributes:
        origin: Ray origin, shape (3,) or (N, 3)
        direction: Ray direction (normalized), shape (3,) or (N, 3)
    """
    origin: torch.Tensor
    direction: torch.Tensor

    def intersect_plane(self, plane: Plane) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Compute ray-plane intersection.
        Returns: (intersects: bool tensor, t: intersection parameter)
        """
```

### 1.4 BBox2D (Optional, may not be needed)
Simple bounding box - could use tuple of Points or dict if needed.

## Phase 2: Detector Class (detector.py)

### 2.1 Core Detector Class

**Key Attributes:**
```python
@dataclass
class DetectorParameters:
    """
    Detector configuration parameters.

    Attributes:
        num_rows: Number of detector rows (pixels in K direction)
        num_cols: Number of detector columns (pixels in J direction)
        pixel_height: Physical height of pixel in mm (K direction)
        pixel_width: Physical width of pixel in mm (J direction)
        beam_center_j: Beam center in J coordinate (pixels)
        beam_center_k: Beam center in K coordinate (pixels)
        orientation: Detector orientation matrix (3, 3) in lab frame
        position: Detector position (3,) in lab frame (rotation center)
        image_basis_j: Image J-axis unit vector (3,) in lab frame
        image_basis_k: Image K-axis unit vector (3,) in lab frame
        image_origin: Image coordinate origin (3,) relative to detector center
    """
    num_rows: int
    num_cols: int
    pixel_height: float
    pixel_width: float
    beam_center_j: float
    beam_center_k: float
    orientation: torch.Tensor  # (3, 3)
    position: torch.Tensor      # (3,)
    image_basis_j: torch.Tensor  # (3,)
    image_basis_k: torch.Tensor  # (3,)
    image_origin: torch.Tensor   # (3,)


class Detector:
    """
    X-ray detector geometry and coordinate transformations.

    The detector represents a 2D pixel array positioned and oriented in 3D lab space.
    It handles transformations between:
    - Lab frame (3D Cartesian coordinates)
    - Detector frame (2D J-K coordinates in mm)
    - Pixel coordinates (2D row-col indices)

    Key Concepts:
    - Lab frame: Global 3D coordinate system (beam propagates along X-axis by convention)
    - Detector frame: 2D coordinate system on detector surface (J: horizontal, K: vertical)
    - Pixel coordinates: Discrete 2D indices (col, row) on detector array
    - Beam center: Point where direct beam intersects detector (in J-K coordinates)
    - Image basis: Unit vectors defining detector plane orientation
    - Rotation center: Detector position in lab frame (typically beam intersection point)
    """
```

### 2.2 Core Methods - Initialization

```python
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
    device: str = 'cpu',
    dtype: torch.dtype = torch.float32
):
    """
    Initialize detector with geometric parameters.

    Args:
        num_rows: Number of pixel rows (K direction)
        num_cols: Number of pixel columns (J direction)
        pixel_height: Pixel height in mm
        pixel_width: Pixel width in mm
        beam_center_j: Beam center J coordinate (pixels)
        beam_center_k: Beam center K coordinate (pixels)
        image_basis_j: J-axis unit vector in detector frame
        image_basis_k: K-axis unit vector in detector frame
        image_origin: Coordinate origin in detector frame
        device: PyTorch device ('cpu' or 'cuda')
        dtype: PyTorch data type
    """
```

### 2.3 Core Methods - Coordinate Transformations

**Critical transformations (all differentiable):**

```python
def lab_to_detector_coordinate(self, lab_pos: torch.Tensor) -> torch.Tensor:
    """
    Transform lab frame position to detector frame (J, K) coordinates.

    Args:
        lab_pos: Position in lab frame, shape (3,) or (N, 3)

    Returns:
        Detector coordinates (J, K), shape (2,) or (N, 2)
    """

def detector_to_lab_coordinate(self, j: torch.Tensor, k: torch.Tensor) -> torch.Tensor:
    """
    Transform detector frame (J, K) to lab frame position.

    Args:
        j: J coordinates, shape () or (N,)
        k: K coordinates, shape () or (N,)

    Returns:
        Lab frame position, shape (3,) or (N, 3)
    """

def lab_to_pixel(self, lab_pos: torch.Tensor) -> torch.Tensor:
    """
    Transform lab frame position to pixel coordinates.

    Args:
        lab_pos: Position in lab frame, shape (3,) or (N, 3)

    Returns:
        Pixel coordinates (col, row), shape (2,) or (N, 2)
        Note: Returns float pixel coordinates (sub-pixel precision)
    """

def pixel_to_lab_coordinate(self, col: torch.Tensor, row: torch.Tensor) -> torch.Tensor:
    """
    Transform pixel coordinates to lab frame position.

    Args:
        col: Column indices (J direction), shape () or (N,)
        row: Row indices (K direction), shape () or (N,)

    Returns:
        Lab frame position, shape (3,) or (N, 3)
    """
```

### 2.4 Core Methods - Detector Transformations

```python
def set_orientation(self, orientation: torch.Tensor) -> None:
    """
    Set detector orientation matrix.

    Args:
        orientation: 3x3 rotation matrix
    """

def set_orientation_euler(self, phi: float, theta: float, psi: float) -> None:
    """
    Set detector orientation using Euler angles (ZYZ convention).

    Args:
        phi, theta, psi: Euler angles in radians
    """

def set_position(self, position: torch.Tensor) -> None:
    """
    Set detector position (rotation center) in lab frame.

    Args:
        position: 3D position vector
    """

def translate(self, translation: torch.Tensor) -> None:
    """
    Translate detector by offset vector.

    Args:
        translation: 3D translation vector
    """

def rotate(self, phi: float, theta: float, psi: float) -> None:
    """
    Rotate detector by Euler angles (applied incrementally).

    Args:
        phi, theta, psi: Euler angles in radians
    """
```

### 2.5 Core Methods - Geometric Queries

```python
def intersect_ray(self, ray: Ray) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Compute ray-detector plane intersection.

    Args:
        ray: Ray object (can be batched)

    Returns:
        (intersects, t): Boolean tensor and intersection parameter
    """

def in_range(self, col: torch.Tensor, row: torch.Tensor) -> torch.Tensor:
    """
    Check if pixel coordinates are within detector bounds.

    Args:
        col: Column indices, shape () or (N,)
        row: Row indices, shape () or (N,)

    Returns:
        Boolean tensor indicating valid pixels
    """

def get_detector_plane(self) -> Plane:
    """
    Get detector plane in lab frame.

    Returns:
        Plane object representing detector surface
    """
```

### 2.6 Properties (Pythonic accessors)

```python
@property
def position(self) -> torch.Tensor:
    """Detector position (rotation center) in lab frame."""

@property
def orientation(self) -> torch.Tensor:
    """Detector orientation matrix."""

@property
def pixel_count(self) -> Tuple[int, int]:
    """Number of pixels (rows, cols)."""

@property
def physical_size(self) -> Tuple[float, float]:
    """Physical detector size in mm (height, width)."""

@property
def beam_center_lab(self) -> torch.Tensor:
    """Beam center position in lab frame coordinates."""
```

## Phase 3: Detector Construction Methods (detector.py)

**REVISED:** No factory class. Use class methods + convenience function instead.

### 3.1 Primary Constructor (Class Method)

The `Detector.from_beam_center()` class method is the **recommended** way to create detectors:

```python
class Detector:
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
        Create detector from beam center specification (RECOMMENDED).

        This is the standard way to create a detector for HEDM experiments.
        The image coordinate origin is automatically calculated from the
        beam center and pixel sizes using the HEDM convention.

        Args:
            num_j_pixels: Number of pixels in J direction (columns)
            num_k_pixels: Number of pixels in K direction (rows)
            beam_center_j: Beam center in J direction (pixels)
            beam_center_k: Beam center in K direction (pixels)
            pixel_width: Pixel width in mm
            pixel_height: Pixel height in mm
            j_unit_vector: Unit vector for J direction (3,)
            k_unit_vector: Unit vector for K direction (3,)
            position: Detector position in lab frame (3,), defaults to origin
            orientation: Orientation matrix (3x3) in lab frame, defaults to identity
            device: PyTorch device
            dtype: PyTorch data type

        Returns:
            Configured Detector instance

        Example:
            >>> detector = Detector.from_beam_center(
            ...     num_j_pixels=2048,
            ...     num_k_pixels=2048,
            ...     beam_center_j=1024.0,
            ...     beam_center_k=1024.0,
            ...     pixel_width=0.2,
            ...     pixel_height=0.2,
            ...     j_unit_vector=torch.tensor([0., 1., 0.]),
            ...     k_unit_vector=torch.tensor([0., 0., 1.]),
            ...     position=torch.tensor([100., 0., 0.])
            ... )
        """
        # Calculate image origin from beam center
        image_origin = cls._compute_image_origin(
            beam_center_j, beam_center_k,
            pixel_width, pixel_height,
            j_unit_vector, k_unit_vector
        )

        # Set defaults
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
        Calculate image coordinate origin from beam center (HEDM convention).

        The image origin is the point where (J=0, K=0) is located relative
        to the detector center (rotation center). Calculated from beam center.

        This matches the C++ implementation in CXDMDetectorFactory::GetCoordOrigin().
        """
        origin = (
            -beam_center_j * pixel_width * j_unit_vector
            - beam_center_k * pixel_height * k_unit_vector
        )
        return origin
```

### 3.2 Convenience Function (C++ API Compatibility)

```python
def make_detector(
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
) -> Detector:
    """
    Create detector with automatic coordinate origin calculation.

    Convenience function that calls Detector.from_beam_center().
    Provided for C++ API compatibility and functional-style usage.

    See Detector.from_beam_center() for full documentation.
    """
    return Detector.from_beam_center(
        num_j_pixels=num_j_pixels,
        num_k_pixels=num_k_pixels,
        beam_center_j=beam_center_j,
        beam_center_k=beam_center_k,
        pixel_width=pixel_width,
        pixel_height=pixel_height,
        j_unit_vector=j_unit_vector,
        k_unit_vector=k_unit_vector,
        position=position,
        orientation=orientation,
        device=device,
        dtype=dtype
    )
```

### 3.3 Parameter Extraction

```python
@dataclass
class DetectorParameters:
    """
    Detector configuration parameters (matches C++ CDetectorInfo).

    This is a serializable representation of detector geometry,
    useful for saving/loading configurations and optimization.
    """
    num_rows: int
    num_cols: int
    pixel_height: float
    pixel_width: float
    beam_center_j: float
    beam_center_k: float
    orientation: torch.Tensor  # (3, 3)
    position: torch.Tensor      # (3,)
    j_unit_vector: torch.Tensor  # (3,)
    k_unit_vector: torch.Tensor  # (3,)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'num_rows': self.num_rows,
            'num_cols': self.num_cols,
            'pixel_height': self.pixel_height,
            'pixel_width': self.pixel_width,
            'beam_center_j': self.beam_center_j,
            'beam_center_k': self.beam_center_k,
            'orientation': self.orientation.cpu().numpy().tolist(),
            'position': self.position.cpu().numpy().tolist(),
            'j_unit_vector': self.j_unit_vector.cpu().numpy().tolist(),
            'k_unit_vector': self.k_unit_vector.cpu().numpy().tolist(),
        }


def get_detector_parameters(detector: Detector) -> DetectorParameters:
    """
    Extract detector parameters as a dataclass.

    Args:
        detector: Detector instance

    Returns:
        DetectorParameters dataclass with all detector settings

    Example:
        >>> params = get_detector_parameters(detector)
        >>> params_dict = params.to_dict()  # For JSON serialization
    """
    return DetectorParameters(
        num_rows=detector.num_rows,
        num_cols=detector.num_cols,
        pixel_height=detector.pixel_height,
        pixel_width=detector.pixel_width,
        beam_center_j=detector.beam_center_j,
        beam_center_k=detector.beam_center_k,
        orientation=detector.orientation.clone(),
        position=detector.position.clone(),
        j_unit_vector=detector._det_frame_image_basis_j.clone(),
        k_unit_vector=detector._det_frame_image_basis_k.clone(),
    )
```

## Phase 4: Utilities and Helpers

### 4.1 Coordinate Conversion Helpers

```python
def pixel_to_detector_j(pixel_col: torch.Tensor, pixel_width: float) -> torch.Tensor:
    """Convert pixel column to J coordinate in mm."""

def pixel_to_detector_k(pixel_row: torch.Tensor, pixel_height: float) -> torch.Tensor:
    """Convert pixel row to K coordinate in mm."""

def detector_j_to_pixel(j: torch.Tensor, pixel_width: float) -> torch.Tensor:
    """Convert J coordinate to pixel column (float)."""

def detector_k_to_pixel(k: torch.Tensor, pixel_height: float) -> torch.Tensor:
    """Convert K coordinate to pixel row (float)."""
```

### 4.2 Batch Support Utilities

```python
def create_detector_batch(
    detector: Detector,
    orientations: torch.Tensor,  # (N, 3, 3)
    positions: torch.Tensor      # (N, 3)
) -> List[Detector]:
    """
    Create multiple detector instances with different orientations/positions.
    Useful for scanning or optimization over detector geometry.
    """
```

## Phase 5: Testing Strategy

### 5.1 Unit Tests (test_detector.py)

```python
class TestDetectorGeometry:
    """Test basic detector geometry and initialization."""

    def test_initialization():
        """Test detector can be created with valid parameters."""

    def test_coordinate_transformations():
        """Test lab ↔ detector ↔ pixel conversions are consistent."""

    def test_identity_transforms():
        """Test round-trip transformations preserve values."""

    def test_rotation_matrices():
        """Test Euler angle conversions."""


class TestDetectorTransformations:
    """Test detector positioning and orientation changes."""

    def test_translation():
        """Test detector translation updates position correctly."""

    def test_rotation():
        """Test detector rotation updates orientation correctly."""

    def test_set_orientation():
        """Test setting orientation directly."""


class TestRayIntersection:
    """Test ray-detector intersection calculations."""

    def test_perpendicular_ray():
        """Test ray perpendicular to detector intersects correctly."""

    def test_parallel_ray():
        """Test parallel ray returns no intersection."""

    def test_angled_ray():
        """Test ray at angle intersects at correct location."""


class TestBatchedOperations:
    """Test batched coordinate transformations."""

    def test_batched_lab_to_pixel():
        """Test batch transformation of multiple points."""

    def test_batched_pixel_to_lab():
        """Test batch transformation of multiple pixels."""


class TestDifferentiability:
    """Test PyTorch autodiff through detector operations."""

    def test_gradient_through_coordinate_transform():
        """Test gradients flow through coordinate transformations."""

    def test_gradient_through_orientation():
        """Test gradients with respect to detector orientation."""

    def test_gradient_through_position():
        """Test gradients with respect to detector position."""
```

### 5.2 Integration Tests (test_detector_integration.py)

```python
class TestCppCompatibility:
    """Compare results with C++ implementation."""

    def test_against_cpp_detector():
        """
        Load detector configuration from C++ config file,
        perform transformations, compare results.
        """

    def test_pixel_positions_match():
        """
        For known detector geometry, verify pixel→lab
        positions match C++ output.
        """


class TestDetectorFactory:
    """Test factory functions produce correct detectors."""

    def test_make_detector():
        """Test make_detector creates properly configured detector."""

    def test_parameter_extraction():
        """Test get_detector_parameters extracts correct values."""
```

### 5.3 Performance Tests

```python
class TestPerformance:
    """Benchmark detector operations."""

    def test_batch_transform_speed():
        """Benchmark batched coordinate transformations."""

    def test_gpu_acceleration():
        """Verify GPU acceleration works correctly."""
```

## Phase 6: Documentation

### 6.1 Module Docstring
- Overview of detector geometry concepts
- Coordinate system conventions
- Usage examples
- Links to relevant papers/documentation

### 6.2 Class and Method Docstrings
- Full parameter descriptions with types and shapes
- Return value descriptions
- Usage examples for complex methods
- Notes on batching behavior
- Notes on differentiability

### 6.3 Examples (examples/detector_usage.py)

```python
"""
Example 1: Basic detector creation and coordinate transformation
Example 2: Detector positioning for HEDM experiment
Example 3: Ray-detector intersection for forward simulation
Example 4: Batch processing multiple points
Example 5: Gradient-based detector calibration (differentiable)
"""
```

## Implementation Notes

### Differentiability Considerations
1. **Use torch operations throughout:** All math operations must use PyTorch functions
2. **Avoid in-place operations where needed:** May break autograd
3. **Handle special cases carefully:** e.g., normalization (divide by zero), acos/asin (domain errors)
4. **Clamp values for numerical stability:** e.g., in Euler angle extraction

### Batching Strategy
1. **Support both single and batched inputs:** Use broadcasting and reshape
2. **Consistent shape conventions:** Document expected shapes clearly
3. **Memory efficiency:** Avoid unnecessary copies

### PyTorch Best Practices
1. **Device management:** Allow user to specify device, respect it throughout
2. **Dtype consistency:** Maintain consistent dtype throughout calculations
3. **No-grad contexts:** Use where appropriate for non-trainable operations

### API Design Principles
1. **Properties over getters:** Use `@property` for read-only attributes
2. **Immutability options:** Consider making detector parameters immutable after creation
3. **Builder pattern:** For complex initialization, consider a builder
4. **Type hints everywhere:** Full type annotations for all public APIs

### Compatibility with C++
1. **Match coordinate conventions:** Ensure same lab frame definition
2. **Match Euler angle convention:** Use ZYZ convention
3. **Match transformation order:** Rotation then translation
4. **Numerical precision:** Use float32 by default (matches C++ Float)

## Dependencies

### Required
- `torch >= 2.0`
- `numpy` (for testing and I/O)
- `pytest` (for testing)

### Optional
- `scipy` (for rotation utilities as reference)
- `matplotlib` (for visualization in examples)

## Migration from C++

### Key Differences
1. **No BBox2D needed:** Can use simple tuples or skip entirely
2. **No separate ImageData class initially:** Focus on detector geometry
3. **Factory is functions, not class:** More Pythonic
4. **Properties instead of Get methods:** More Pythonic
5. **No explicit const:** Use immutability through design
6. **Batching is first-class:** Support batch operations natively

### Not Porting (Initially)
1. **AddDirectBeam:** This is about image manipulation, not detector geometry
2. **GetPixelExtent(voxel):** This is voxel-specific, can be added later if needed
3. **State saving (oPrevOrientationMatrix):** Not needed initially

## Success Criteria

1. ✅ All unit tests pass
2. ✅ Integration test matches C++ detector within numerical precision
3. ✅ Batched operations work correctly
4. ✅ Gradients can be computed through all transformations
5. ✅ GPU execution works correctly
6. ✅ Performance is acceptable (benchmark defined)
7. ✅ Documentation is complete and clear
8. ✅ Examples run without errors

## Timeline Estimate (REVISED)

- **Phase 1 (Geometry primitives):** 4-5 hours
- **Phase 2 (Core Detector):** 8-10 hours
- **Phase 3 (Construction methods):** 2-3 hours
- **Phase 4 (Utilities):** 1-2 hours (reduced - minimal helpers)
- **Phase 5 (Testing):** 6-8 hours
- **Phase 6 (Documentation):** 3-4 hours

**Total: 24-32 hours** (~3-4 days of focused work)

**Reduced from original estimate by:**
- Removing factory class complexity
- Simplifying utilities
- Using class methods (more Pythonic)

## Open Questions for Review (UPDATED)

### Resolved (see DETECTOR_ANALYSIS.md):
1. ✅ **Factory pattern:** Use class methods, NOT factory class
2. ✅ **Geometry utilities:** Create `geometry.py` with minimal PyTorch implementations
3. ✅ **DetectorFile/Data:** DO NOT port, use Python config libraries

### Remaining Questions:

1. **Batching strategy:** Should Detector support batched internal state, or just batched operations?
   - **Recommendation:** Batched operations only (like diffraction_core.py)
   - Detector instance represents single detector
   - Operations accept batched inputs (N, 3) for positions

2. **Immutability:** Should detector parameters be immutable after creation?
   - **Recommendation:** Mutable for flexibility, but provide immutable views if needed
   - Allow `set_orientation()`, `translate()` for calibration workflows
   - Document which operations modify state

3. **Device handling:** Auto-detect device from input tensors, or require explicit specification?
   - **Recommendation:** Explicit device in constructor (follows PyTorch convention)
   - All tensors must be on same device
   - Provide `.to(device)` method for moving detector

4. **Numerical precision:** Float32 (matches C++) vs Float64 (more accurate)?
   - **Recommendation:** Float32 default (matches C++ and diffraction_core.py)
   - Allow user to specify dtype if needed
   - Document precision limits for coordinate transforms

5. **Euler convention:** Confirm ZYZ is correct, document clearly
   - **TODO:** Verify with C++ code (Detector.cpp line 171)
   - Document convention in geometry.py docstrings
   - Add tests comparing with C++ outputs

6. **Lab frame convention:** Confirm beam along X-axis is correct
   - **TODO:** Verify with C++ code and ExperimentSetup
   - Document in detector.py module docstring
   - Add diagram in documentation

7. **Image coordinate system:** Confirm J=horizontal, K=vertical convention
   - **TODO:** Verify with C++ code
   - J = columns, K = rows (appears consistent)
   - Document clearly in Detector class docstring

8. **Integration with existing code:** How will this interface with DiffractionCore and Sample?
   - Detector will be used in forward simulation
   - Need to check Sample.py porting plan for coordinate conventions
   - Ensure consistent device/dtype handling

9. **Configuration file I/O:** Should we support reading/writing detector configs from files?
   - **Recommendation:** Add helper functions later as needed
   - Use JSON/YAML for Python configs (not custom parser)
   - `DetectorParameters.to_dict()` already supports serialization
