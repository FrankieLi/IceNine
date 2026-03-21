"""
SO(3) uniform sampling via Sukharev grids on the 3-sphere.

Direct port of XDM++/libXDM/Sampling.h/cpp (Yershova & LaValle, 2003).
Provides deterministic low-discrepancy orientation grids for discrete search
in reconstruction.

C++ Reference:
    XDM++/libXDM/Sampling.h   — CQuaternionGrid class
    XDM++/libXDM/Sampling.cpp — Implementation
    Src/DiscreteSearch.cpp     — GenerateLocalGrid wrapper
"""

import math
from typing import List, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation


# ---------------------------------------------------------------------------
# Constants matching C++ Sampling.h
# ---------------------------------------------------------------------------
_NUM_VERTICES = 8
_CUBE_DIM = 3
_NUM_HYPERFACES = 4

# Gray code vertex ordering for Sukharev grid (Sampling.h:120-130)
_VERTEX_ORDER_GRAY_CODE = np.array([
    [0, 0, 0],
    [1, 1, 1],
    [0, 0, 1],
    [1, 1, 0],
    [0, 1, 0],
    [1, 0, 1],
    [0, 1, 1],
    [1, 0, 0],
], dtype=np.float64)


# ---------------------------------------------------------------------------
# Sukharev grid sequence (Lindemann & LaValle, 2003)
# ---------------------------------------------------------------------------

def get_sukarev_grid_point(state: int, ndim: int = _CUBE_DIM) -> np.ndarray:
    """
    Map integer state to a point in [0,1)^ndim using Gray-code Sukharev sequence.

    C++ Reference: Sampling.cpp:123-149 GetSukarevGridPoint
    """
    sample = np.zeros(ndim, dtype=np.float64)
    current_factor = 0.5
    current_index = state % _NUM_VERTICES
    blah = state // _NUM_VERTICES

    while blah > 0:
        for j in range(ndim):
            sample[j] += current_factor * _VERTEX_ORDER_GRAY_CODE[current_index, j]
        current_factor *= 0.5
        current_index = blah % _NUM_VERTICES
        blah = blah // _NUM_VERTICES

    for j in range(ndim):
        sample[j] += current_factor * _VERTEX_ORDER_GRAY_CODE[current_index, j]

    return sample


def get_layered_sukarev_grid_point(state: int, ndim: int = _CUBE_DIM) -> np.ndarray:
    """
    Layered Sukharev grid: incremental, deterministic low-discrepancy sequence.

    C++ Reference: Sampling.cpp:82-107 GetLayeredSukharevGridPoint
    """
    state += 1  # C++ does nState++
    level = 0
    offset = 0.5
    samp_index = 2.0 ** (ndim * level)

    while state > samp_index:
        state -= int(samp_index)
        offset = offset / 2.0
        level += 1
        samp_index = 2.0 ** (ndim * level)

    state -= 1  # C++ does nState--

    blah = np.full(ndim, offset, dtype=np.float64)
    sample = blah + get_sukarev_grid_point(state, ndim)
    return sample


def make_sukarev_grid_points(level: int, side_width: float = 1.0) -> np.ndarray:
    """
    Generate a regular 3D grid of Sukharev points at given resolution level.

    Returns (N, 3) array where N = (2^level)^3.
    Points are centered in cells of size side_width / 2^level.

    C++ Reference: Sampling.cpp:40-67 MakeSukarevGridPoints
    """
    n_div = 2 ** level
    scale = side_width / n_div
    n_points_per_axis = n_div
    offset = 0.5 * scale

    points = []
    for nx in range(n_points_per_axis):
        for ny in range(n_points_per_axis):
            for nz in range(n_points_per_axis):
                point = np.array([
                    nx * scale + offset,
                    ny * scale + offset,
                    nz * scale + offset,
                ], dtype=np.float64)
                points.append(point)

    return np.array(points, dtype=np.float64)


# ---------------------------------------------------------------------------
# SLERP (Spherical Linear Interpolation)
# ---------------------------------------------------------------------------

def slerp(q1: np.ndarray, q2: np.ndarray, t: float) -> np.ndarray:
    """
    Spherical linear interpolation between two unit quaternions.

    Quaternion convention: [w, x, y, z] (scalar-first).

    C++ Reference: Quaternion.cpp:404-455 Interpolate
    """
    cos_theta = np.dot(q1, q2)

    # If on opposite hemisphere, flip
    flip = False
    if cos_theta < 0.0:
        cos_theta = -cos_theta
        flip = True

    if 1.0 - cos_theta < 0.01:
        # Near-parallel: linear interpolation
        beta = 1.0 - t
        alpha = t
    else:
        theta = math.acos(cos_theta)
        sin_theta = math.sin(theta)
        beta = math.sin(theta - t * theta) / sin_theta
        alpha = math.sin(t * theta) / sin_theta

    if flip:
        alpha = -alpha

    result = beta * q1 + alpha * q2
    return result


# ---------------------------------------------------------------------------
# CQuaternionGrid — uniform sampling on SO(3)
# ---------------------------------------------------------------------------

class QuaternionGrid:
    """
    Generates uniform orientation grids on SO(3) using quaternion parameterization.

    Uses the Yershova & LaValle Sukharev grid mapped to 4 hyperfaces of the
    upper hemi-hypersphere of S^3.

    C++ Reference: Sampling.h/cpp CQuaternionGrid
    """

    def __init__(self):
        self._hyperface_vertices = self._initialize_faces()

    @staticmethod
    def _initialize_faces() -> np.ndarray:
        """
        Initialize 4 hyperface vertex sets, 8 vertices each.
        Each vertex is a unit quaternion [w, x, y, z].

        C++ Reference: Sampling.cpp:170-213 InitializeFaces
        """
        faces = np.zeros((_NUM_HYPERFACES, _NUM_VERTICES, 4), dtype=np.float64)

        # HyperFace 0 (identity face: w = +0.5)
        faces[0] = [
            [0.5, -0.5, -0.5, -0.5],
            [0.5,  0.5, -0.5, -0.5],
            [0.5, -0.5, -0.5,  0.5],
            [0.5,  0.5, -0.5,  0.5],
            [0.5, -0.5,  0.5, -0.5],
            [0.5,  0.5,  0.5, -0.5],
            [0.5, -0.5,  0.5,  0.5],
            [0.5,  0.5,  0.5,  0.5],
        ]

        # HyperFace 1 (x = +0.5)
        faces[1] = [
            [-0.5,  0.5, -0.5, -0.5],
            [ 0.5,  0.5, -0.5, -0.5],
            [-0.5,  0.5, -0.5,  0.5],
            [ 0.5,  0.5, -0.5,  0.5],
            [-0.5,  0.5,  0.5, -0.5],
            [ 0.5,  0.5,  0.5, -0.5],
            [-0.5,  0.5,  0.5,  0.5],
            [ 0.5,  0.5,  0.5,  0.5],
        ]

        # HyperFace 2 (y = +0.5)
        faces[2] = [
            [-0.5, -0.5,  0.5, -0.5],
            [ 0.5, -0.5,  0.5, -0.5],
            [-0.5, -0.5,  0.5,  0.5],
            [ 0.5, -0.5,  0.5,  0.5],
            [-0.5,  0.5,  0.5, -0.5],
            [ 0.5,  0.5,  0.5, -0.5],
            [-0.5,  0.5,  0.5,  0.5],
            [ 0.5,  0.5,  0.5,  0.5],
        ]

        # HyperFace 3 (z = +0.5)
        faces[3] = [
            [-0.5, -0.5, -0.5,  0.5],
            [ 0.5, -0.5, -0.5,  0.5],
            [-0.5, -0.5,  0.5,  0.5],
            [ 0.5, -0.5,  0.5,  0.5],
            [-0.5,  0.5, -0.5,  0.5],
            [ 0.5,  0.5, -0.5,  0.5],
            [-0.5,  0.5,  0.5,  0.5],
            [ 0.5,  0.5,  0.5,  0.5],
        ]

        return faces

    def barycentric_to_quaternion(
        self, face_index: int, alpha: float, beta: float, gamma: float
    ) -> np.ndarray:
        """
        Map barycentric coordinates on a hyperface to a unit quaternion.

        Uses tri-linear SLERP across the 8 vertices of the hyperface cube.
        Result is normalized to the positive-w hemisphere.

        Args:
            face_index: Hyperface index [0, 3]
            alpha, beta, gamma: Barycentric coordinates in [0, 1]

        Returns:
            Unit quaternion [w, x, y, z]

        C++ Reference: Sampling.cpp:222-239 BarycentricToQuaternion
        """
        v = self._hyperface_vertices[face_index]

        x1 = slerp(v[0], v[1], alpha)
        x2 = slerp(v[2], v[3], alpha)
        x3 = slerp(v[4], v[5], alpha)
        x4 = slerp(v[6], v[7], alpha)

        y1 = slerp(x1, x3, beta)
        y2 = slerp(x2, x4, beta)

        result = slerp(y1, y2, gamma)

        # ToConvention: ensure positive hemisphere (w >= 0)
        if result[0] < 0:
            result = -result

        return result

    def get_near_identity_point(self, x: float, y: float, z: float) -> np.ndarray:
        """
        Get a quaternion near identity by offsetting from center of face 0.

        C++ Reference: Sampling.cpp:448-459 GetNearIdentityPoint
        """
        origin = np.array([0.5, 0.5, 0.5])
        new_pos = np.array([x, y, z]) + origin
        return self.barycentric_to_quaternion(0, new_pos[0], new_pos[1], new_pos[2])

    def get_grid_by_count(self, n_points: int) -> np.ndarray:
        """
        Generate n_points uniformly distributed on SO(3).

        Points are distributed across 4 hyperfaces using the layered
        Sukharev grid sequence.

        Args:
            n_points: Number of orientations to generate

        Returns:
            (n_points, 4) array of unit quaternions [w, x, y, z]

        C++ Reference: Sampling.cpp:392-408 GetGrid(Int)
        """
        quats = np.zeros((n_points, 4), dtype=np.float64)
        for i in range(n_points):
            face_index = i % _NUM_HYPERFACES
            index = i // _NUM_HYPERFACES
            v_interp = get_layered_sukarev_grid_point(index)
            quats[i] = self.barycentric_to_quaternion(
                face_index, v_interp[0], v_interp[1], v_interp[2]
            )
        return quats

    def get_grid_by_dispersion(self, max_distance: float) -> np.ndarray:
        """
        Generate uniform SO(3) grid with bounded dispersion.

        Uses Proposition 4.5 of Yershova & LaValle:
            n ~ 0.5 * [(2*pi / d_rho)^3 - 1]

        Args:
            max_distance: Maximum dispersion (Euclidean distance on S^3)

        Returns:
            (N, 4) array of unit quaternions

        C++ Reference: Sampling.cpp:430-441 GetGrid(Float)
        """
        n_points = math.ceil(
            0.5 * ((2.0 * math.pi / max_distance) ** 3 - 1.0)
        )
        return self.get_grid_by_count(n_points)

    def get_random_local_grid(
        self, max_distance: float, n_points: int, rng: Optional[np.random.Generator] = None
    ) -> np.ndarray:
        """
        Generate random local grid around identity.

        Note: C++ creates a local unseeded RNG per call (deterministic).
        We preserve this behavior with a default unseeded generator.

        C++ Reference: Sampling.cpp:466-479 GetRandomLocalGrid
        """
        if rng is None:
            rng = np.random.default_rng(seed=None)

        quats = np.zeros((n_points, 4), dtype=np.float64)
        half = max_distance / 2.0
        for i in range(n_points):
            x = rng.uniform(-half, half)
            y = rng.uniform(-half, half)
            z = rng.uniform(-half, half)
            quats[i] = self.get_near_identity_point(x, y, z)
        return quats

    def get_structured_local_grid(
        self, side_length: float, level: int = 2
    ) -> np.ndarray:
        """
        Generate structured local grid around identity quaternion.

        Produces an approximately Euclidean grid in SO(3) centered at identity.
        The grid covers a misorientation range of `side_length` radians.

        Args:
            side_length: Angular diameter in radians (misorientation space)
            level: Resolution level (grid spacing ~ side_length / 2^level)

        Returns:
            (N, 4) array of unit quaternions where N = (2^level)^3

        C++ Reference: Sampling.cpp:500-521 GetStructuredLocalGrid
        """
        # Project from angular space to parametric space
        side_length = 2.0 * math.sin(math.pi / 8.0) * math.tan(side_length / 2.0)

        offsets = make_sukarev_grid_points(level, side_length)
        center = np.full(3, side_length / 2.0)
        origin = np.array([0.5, 0.5, 0.5]) - center

        quats = np.zeros((len(offsets), 4), dtype=np.float64)
        for i, offset in enumerate(offsets):
            new_pos = offset + origin
            quats[i] = self.barycentric_to_quaternion(
                0, new_pos[0], new_pos[1], new_pos[2]
            )
        return quats


# ---------------------------------------------------------------------------
# High-level helpers
# ---------------------------------------------------------------------------

def generate_local_grid(
    angular_coverage: float, level: int
) -> np.ndarray:
    """
    Generate local orientation grid as rotation matrices.

    This is the Python equivalent of OrientationSearch::Utilities::GenerateLocalGrid.
    Produces small rotation matrices near identity for local search refinement.

    Args:
        angular_coverage: Diameter of the grid in radians
        level: Resolution level (0 = coarsest, higher = finer)

    Returns:
        (N, 3, 3) array of rotation matrices

    C++ Reference: DiscreteSearch.cpp:64-79 GenerateLocalGrid
    """
    grid_gen = QuaternionGrid()
    quats = grid_gen.get_structured_local_grid(angular_coverage, level)
    # Convert quaternions [w, x, y, z] to rotation matrices via scipy
    # scipy uses [x, y, z, w] convention, so reorder
    quats_scipy = quats[:, [1, 2, 3, 0]]
    rotations = Rotation.from_quat(quats_scipy)
    return rotations.as_matrix()


def generate_local_grid_multi_level(
    angular_coverage: float, min_level: int, max_level: int
) -> np.ndarray:
    """
    Generate local grid across multiple resolution levels.

    C++ Reference: DiscreteSearch.cpp:84-90 GenerateLocalGrid (multi-level)
    """
    grids = []
    for level in range(min_level, max_level + 1):
        grids.append(generate_local_grid(angular_coverage, level))
    return np.concatenate(grids, axis=0)


def load_fundamental_zone_file(filename: str) -> np.ndarray:
    """
    Load fundamental zone orientations from file.

    File format: text file with Euler angle triples in degrees (phi, theta, psi),
    one triple per line (space/comma/tab separated). Active ZXZ Bunge convention.

    Returns:
        (N, 3, 3) array of rotation matrices

    C++ Reference: InitFilesIO.cpp:406-442 ReadFundamentalZoneFile
                   ReconstructionSetup.cpp:129-140 (usage)
    """
    from icenine.geometry import active_euler_matrix

    angles_deg = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.replace(',', ' ').split()
            if len(parts) >= 3:
                phi = float(parts[0])
                theta = float(parts[1])
                psi = float(parts[2])
                angles_deg.append((phi, theta, psi))

    n = len(angles_deg)
    matrices = np.zeros((n, 3, 3), dtype=np.float64)
    for i, (phi, theta, psi) in enumerate(angles_deg):
        # active_euler_matrix takes radians, returns torch.Tensor
        mat = active_euler_matrix(
            math.radians(phi), math.radians(theta), math.radians(psi)
        )
        matrices[i] = mat.numpy()

    return matrices


def reduce_to_fundamental_zone(
    q: np.ndarray, symmetry_quats: np.ndarray
) -> np.ndarray:
    """
    Reduce a quaternion to the fundamental zone under crystal symmetry.

    Applies all symmetry operators and selects the result with largest |w|
    (closest to identity rotation).

    Args:
        q: Quaternion [w, x, y, z]
        symmetry_quats: (N_sym, 4) array of symmetry operator quaternions [w, x, y, z]

    Returns:
        Reduced quaternion [w, x, y, z]

    C++ Reference: Symmetry.tmpl.cpp:79-98 ReduceToFundamentalZone
    """
    best = q.copy()
    max_abs_w = abs(q[0])

    for sym_q in symmetry_quats:
        # Quaternion product: q * sym_q
        product = _quat_multiply(q, sym_q)
        if abs(product[0]) > max_abs_w:
            best = product
            max_abs_w = abs(best[0])

    # Ensure positive-w hemisphere (C++ ToConvention)
    if best[0] < 0:
        best = -best

    return best


def get_misorientation(
    q1: np.ndarray, q2: np.ndarray, symmetry_quats: np.ndarray
) -> float:
    """
    Compute symmetry-reduced misorientation angle between two orientations.

    Args:
        q1, q2: Quaternions [w, x, y, z]
        symmetry_quats: (N_sym, 4) symmetry operator quaternions

    Returns:
        Misorientation angle in radians

    C++ Reference: Symmetry.tmpl.cpp:61-71 GetMisorientation (quaternion version)
    """
    q_prod = _quat_multiply(_quat_inverse(q1), q2)
    q_fz = reduce_to_fundamental_zone(q_prod, symmetry_quats)
    angle = 2.0 * math.acos(min(1.0, q_fz[0]))
    return angle


def is_in_fundamental_zone(
    q: np.ndarray, symmetry_quats: np.ndarray
) -> bool:
    """
    Check if quaternion is already in the fundamental zone.

    C++ Reference: Symmetry.tmpl.cpp:105-116 IsInFundamentalZone
    """
    q_fz = reduce_to_fundamental_zone(q, symmetry_quats)
    diff = q - q_fz
    return np.linalg.norm(diff) < 0.01


# ---------------------------------------------------------------------------
# Quaternion arithmetic helpers (scalar-first [w, x, y, z] convention)
# ---------------------------------------------------------------------------

def _quat_multiply(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """Hamilton product of two quaternions [w, x, y, z]."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
        w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2,
    ], dtype=np.float64)


def _quat_inverse(q: np.ndarray) -> np.ndarray:
    """Inverse of a unit quaternion (conjugate)."""
    return np.array([q[0], -q[1], -q[2], -q[3]], dtype=np.float64)


def quaternion_to_matrix(q: np.ndarray) -> np.ndarray:
    """
    Convert quaternion [w, x, y, z] to 3x3 rotation matrix.
    Uses scipy for numerical stability.
    """
    r = Rotation.from_quat([q[1], q[2], q[3], q[0]])  # scipy uses [x, y, z, w]
    return r.as_matrix()


def matrix_to_quaternion(m: np.ndarray) -> np.ndarray:
    """
    Convert 3x3 rotation matrix to quaternion [w, x, y, z].
    Uses scipy for numerical stability.
    """
    r = Rotation.from_matrix(m)
    q_scipy = r.as_quat()  # [x, y, z, w]
    return np.array([q_scipy[3], q_scipy[0], q_scipy[1], q_scipy[2]])
