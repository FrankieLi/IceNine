# Detector Geometry - Reusing Existing Code

## Critical Discovery: Euler Angle Functions Already Exist!

### What's Already Implemented in mic_file.py

The following functions are **already implemented and tested**:

1. **`euler_to_matrix()`** (NumPy version)
   - Converts Bunge Euler angles (degrees) to rotation matrix
   - Uses scipy Rotation
   - Convention: ZXZ intrinsic (matches C++ BuildActiveEulerMatrix)
   - Reference: XDM++/libXDM/3dMath.cpp lines 152-174

2. **`matrix_to_euler()`** (NumPy version)
   - Converts rotation matrix to Bunge Euler angles (degrees)
   - Uses scipy Rotation
   - Convention: ZXZ intrinsic (matches C++ GetEulerAngles)
   - Reference: XDM++/libXDM/3dMath.cpp lines 180-209

3. **`euler_to_matrix_torch()`** (PyTorch version) ✅
   - **Differentiable, batched, PyTorch implementation**
   - Exactly matches C++ formula (verified line-by-line)
   - Supports both scalar and batched inputs
   - Already tested against C++ implementation

### What's Missing

4. **`matrix_to_euler_torch()`** (PyTorch version) ❌
   - NOT yet implemented
   - Would be the differentiable inverse of euler_to_matrix_torch()
   - Needed for Detector.get_orientation_euler() method

## Verification: Formula Matches C++

Comparing [3dMath.cpp:160-170](../XDM++/libXDM/3dMath.cpp#L160-L170) with [mic_file.py:782-792](mic_file.py#L782-L792):

**C++ (BuildActiveEulerMatrix):**
```cpp
m[0][0] = fCosPhi * fCosPsi - fSinPhi * fCosTheta * fSinPsi;
m[1][0] = fSinPhi * fCosPsi + fCosPhi * fCosTheta * fSinPsi;
m[2][0] = fSinTheta * fSinPsi;
// ... (9 lines total)
```

**Python (euler_to_matrix_torch):**
```python
m00 = cos_phi1 * cos_phi2 - sin_phi1 * cos_Phi * sin_phi2
m10 = sin_phi1 * cos_phi2 + cos_phi1 * cos_Phi * sin_phi2
m20 = sin_Phi * sin_phi2
# ... (9 lines total)
```

✅ **EXACT MATCH** (variable naming: phi1=Phi, Phi=Theta, phi2=Psi)

## Euler Convention Clarification

**C++ Comment vs Reality:**
- C++ comment says "ZYZ convention" (line 176 of 3dMath.cpp)
- mic_file.py says "ZXZ convention"
- **The actual FORMULA is what matters, and it's identical**

The confusion comes from Euler angle naming:
- Bunge convention (materials science): phi1, Phi, phi2
- Physics convention: phi, theta, psi
- Rotation sequence: Z-X-Z (intrinsic rotations)

**For our purposes:** The implementation in mic_file.py is **correct and validated**.

## Revised Strategy for Detector Port

### Option 1: Import from mic_file.py (SIMPLEST)

```python
# detector.py
from icenine.mic_file import euler_to_matrix_torch

class Detector:
    def set_orientation_euler(self, phi: float, theta: float, psi: float):
        """Set orientation using Euler angles (Bunge convention)."""
        # Convert degrees to tensor
        phi_t = torch.tensor(phi, dtype=self.dtype, device=self.device)
        theta_t = torch.tensor(theta, dtype=self.dtype, device=self.device)
        psi_t = torch.tensor(psi, dtype=self.dtype, device=self.device)

        # Use existing function
        orientation = euler_to_matrix_torch(phi_t, theta_t, psi_t)
        self.set_orientation(orientation)
```

**Pros:**
- ✅ Zero code duplication
- ✅ Reuses tested code
- ✅ Minimal changes
- ✅ Fast to implement

**Cons:**
- ❌ Creates dependency: detector.py → mic_file.py
- ❌ mic_file.py is conceptually about file I/O, not geometry
- ❌ Confusing module organization

### Option 2: Refactor to geometry.py (CLEANER)

Move Euler angle functions from mic_file.py to a new geometry.py module:

```python
# geometry.py
"""
Geometric primitives and transformations for X-ray diffraction.

PyTorch-based implementations for differentiability and GPU support.
"""

import torch
import numpy as np
from scipy.spatial.transform import Rotation
from typing import Tuple


# ============================================================================
# Euler Angle Conversions (Bunge Convention: ZXZ intrinsic)
# ============================================================================

def euler_to_matrix_torch(
    phi1_deg: torch.Tensor,
    Phi_deg: torch.Tensor,
    phi2_deg: torch.Tensor
) -> torch.Tensor:
    """
    Convert batched Euler angles (degrees) to rotation matrices (PyTorch).

    Differentiable version for neural network integration.
    Uses Bunge (ZXZ intrinsic) convention matching IceNine C++.

    C++ reference: XDM++/libXDM/3dMath.cpp BuildActiveEulerMatrix

    Args:
        phi1_deg: First Euler angles in degrees, shape (N,) or scalar
        Phi_deg: Second Euler angles in degrees, shape (N,) or scalar
        phi2_deg: Third Euler angles in degrees, shape (N,) or scalar

    Returns:
        Rotation matrices, shape (N, 3, 3) or (3, 3) if scalar
    """
    # [Move implementation from mic_file.py]
    ...

def matrix_to_euler_torch(
    rotation_matrix: torch.Tensor,
    epsilon: float = 1e-6
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Convert rotation matrix to Bunge Euler angles (degrees) - PyTorch version.

    Differentiable inverse of euler_to_matrix_torch().

    C++ reference: XDM++/libXDM/3dMath.cpp GetEulerAngles

    Args:
        rotation_matrix: Rotation matrix, shape (3, 3) or (N, 3, 3)
        epsilon: Small value for numerical stability

    Returns:
        Tuple of (phi1_deg, Phi_deg, phi2_deg) in degrees

    Note:
        This is the PyTorch implementation - NEEDS TO BE IMPLEMENTED.
    """
    # TODO: Port from C++ GetEulerAngles (lines 180-209 of 3dMath.cpp)
    # Handle gimbal lock cases (lines 187-198)
    # Extract angles from matrix elements (lines 201-203)
    # Normalize to [0, 360] for phi1/phi2, [0, 180] for Phi
    raise NotImplementedError("matrix_to_euler_torch not yet implemented")

# Also move euler_to_matrix() and matrix_to_euler() (NumPy versions)
# for backward compatibility

# ============================================================================
# Geometric Primitives
# ============================================================================

@dataclass
class Plane:
    """Plane in 3D space: Ax + By + Cz + D = 0"""
    coeffs: torch.Tensor  # (4,) or (N, 4)
    ...

@dataclass
class Ray:
    """Ray in 3D space: P(t) = origin + t * direction"""
    origin: torch.Tensor
    direction: torch.Tensor
    ...
```

Then update mic_file.py:
```python
# mic_file.py
from .geometry import euler_to_matrix_torch, matrix_to_euler, euler_to_matrix

# Remove the duplicate implementations
# Keep mic_file.py focused on MIC file I/O
```

And detector.py:
```python
# detector.py
from .geometry import euler_to_matrix_torch, Plane, Ray

class Detector:
    def set_orientation_euler(self, phi: float, theta: float, psi: float):
        ...
```

**Pros:**
- ✅ Clean separation of concerns
- ✅ Reusable across modules
- ✅ Natural place for Plane, Ray classes
- ✅ Better long-term maintainability

**Cons:**
- ❌ Requires refactoring mic_file.py
- ❌ Potential for breaking existing code if not careful
- ❌ More work upfront

### Option 3: Duplicate in detector.py (NOT RECOMMENDED)

Create separate implementations in detector.py.

**Pros:**
- ✅ No dependencies

**Cons:**
- ❌ Code duplication
- ❌ Maintenance nightmare
- ❌ Tests duplicated
- ❌ Bug fixes need to be applied twice

## Recommendation: Option 2 (Refactor to geometry.py)

**Rationale:**
1. **Code organization**: geometry.py is the natural home for these utilities
2. **Reusability**: Both mic_file.py and detector.py need Euler conversions
3. **Future-proof**: Sample.py will also need rotation matrices
4. **Consistency**: Matches the pattern of diffraction_core.py (physics in own module)
5. **Clean dependencies**:
   - geometry.py has no internal dependencies
   - mic_file.py imports from geometry
   - detector.py imports from geometry

**Migration Plan:**
1. Create geometry.py with Euler functions (move from mic_file.py)
2. Add Plane and Ray classes
3. Implement matrix_to_euler_torch() (missing piece)
4. Update mic_file.py to import from geometry (backward compatible)
5. Update mic_file.py tests to verify nothing broke
6. Detector.py imports from geometry

**Backward Compatibility:**
Keep public API of mic_file.py unchanged:
```python
# mic_file.py still exports these for existing users
from .geometry import euler_to_matrix, matrix_to_euler, euler_to_matrix_torch

__all__ = ['MicFile', 'euler_to_matrix', 'matrix_to_euler', 'euler_to_matrix_torch']
```

## Updated Implementation Plan

### Phase 1: Create geometry.py (4-5 hours)

1. **Move existing functions from mic_file.py** (1 hour)
   - euler_to_matrix()
   - matrix_to_euler()
   - euler_to_matrix_torch()
   - Add comprehensive docstrings with C++ references

2. **Implement matrix_to_euler_torch()** (2-3 hours)
   - Port from C++ GetEulerAngles (lines 180-209)
   - Handle gimbal lock edge cases
   - Make differentiable
   - Support batching

3. **Add Plane and Ray classes** (1 hour)
   - Simple dataclasses with key methods
   - ray_plane_intersection() function

4. **Write tests** (included in phase 5)

### Phase 2: Update mic_file.py (0.5 hours)

1. Replace implementations with imports from geometry
2. Keep exports for backward compatibility
3. Verify existing tests still pass

### Phase 3: Detector implementation (8-10 hours)

Use geometry.euler_to_matrix_torch() and geometry.Plane/Ray

**Total time: Same as before (24-32 hours), but better organized**

## What This Means for DETECTOR_PORT_PLAN.md

**Update the plan:**
1. ✅ Euler conversions mostly done (reuse from mic_file.py)
2. ⚠️ Need to implement matrix_to_euler_torch() (missing piece)
3. ✅ Refactor into geometry.py for clean organization
4. ✅ Add Plane and Ray classes

**Reduces implementation effort:**
- Euler → matrix: **Already done and tested** ✅
- Matrix → Euler: **Need to implement** (2-3 hours)
- Plane/Ray: **Need to implement** (1 hour)

**Net result:** Same timeline, but more robust and maintainable.
