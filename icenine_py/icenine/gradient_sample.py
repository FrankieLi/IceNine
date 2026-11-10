"""A sample whose grains hold smooth orientation gradients (autoresearch AR3, step 3).

The simulated test layer has exactly one orientation per grain. A real strain-free grain can hold
micro-gradients, so every voxel may have its own orientation. ``gradient_orientations`` builds a
copy of such a layer's truth in which each grain keeps its old orientation at its centre and turns
smoothly with position inside the grain.

The model is a choice, not a measurement: the rotation vector is linear in position. Only numpy
is used; no other module changes.
"""

from typing import Dict

import numpy as np

LEVELS_DEG = (0.0, 0.1, 0.2, 0.5, 1.0)


def rotation_matrices(rotvec: np.ndarray) -> np.ndarray:
    """Rotation matrices (n, 3, 3) of rotation vectors (n, 3) in radians (Rodrigues)."""
    rv = np.asarray(rotvec, dtype=np.float64)
    th = np.linalg.norm(rv, axis=1)
    safe = np.where(th > 0, th, 1.0)
    k = rv / safe[:, None]
    K = np.zeros((len(rv), 3, 3))
    K[:, 0, 1], K[:, 0, 2] = -k[:, 2], k[:, 1]
    K[:, 1, 0], K[:, 1, 2] = k[:, 2], -k[:, 0]
    K[:, 2, 0], K[:, 2, 1] = -k[:, 1], k[:, 0]
    s, c = np.sin(th)[:, None, None], np.cos(th)[:, None, None]
    out = np.eye(3)[None] + s * K + (1.0 - c) * (K @ K)
    out[th == 0] = np.eye(3)
    return out


def gradient_orientations(
    centres: np.ndarray,
    grain_ids: np.ndarray,
    grain_orientations: Dict[int, np.ndarray],
    levels_deg: Dict[int, float],
    rng: np.random.Generator,
) -> np.ndarray:
    """Orientation (n, 3, 3) of every voxel of a layer whose grains hold gradients.

    ``centres`` (n, 2) or (n, 3): voxel centres (only the first two coordinates are used; the
    layer is flat). ``grain_ids`` (n,): grain of each voxel. ``grain_orientations``: grain id ->
    (3, 3) orientation at the grain centre. ``levels_deg``: grain id -> spread in degrees. ``rng``:
    numpy generator.

    For each grain, in increasing id order, one 3 by 2 matrix ``M`` of standard normal numbers is
    drawn (also for a grain of level 0 or of one voxel, so the draws do not depend on the levels).
    With ``c`` the mean of the grain's voxel centres, the rotation vector at position ``x`` is
    ``M (x - c)``, scaled so that its largest length over the grain's voxels equals the level
    (degrees, converted to radians). The voxel orientation is that small rotation LEFT-multiplied
    onto the grain orientation: ``R_voxel = exp([r]x) @ R_grain``. The misorientation between
    ``R_voxel`` and ``R_grain`` is then the rotation angle, which is at most the level (equal at
    the farthest voxel). A grain of one voxel, and a grain of level 0, keep their orientation
    exactly (bit for bit).
    """
    xy = np.asarray(centres, dtype=np.float64)[:, :2]
    ids = np.asarray(grain_ids)
    out = np.empty((len(ids), 3, 3), dtype=np.float64)
    for g in sorted(int(x) for x in np.unique(ids)):
        M = rng.standard_normal((3, 2))
        sel = np.nonzero(ids == g)[0]
        Rg = np.asarray(grain_orientations[g], dtype=np.float64)
        out[sel] = Rg
        level = float(levels_deg[g])
        if level == 0.0 or len(sel) < 2:
            continue
        d = xy[sel] - xy[sel].mean(axis=0)
        rv = d @ M.T  # (m, 3)
        top = np.linalg.norm(rv, axis=1).max()
        if top <= 0.0:
            continue
        rv = rv * (np.radians(level) / top)
        out[sel] = rotation_matrices(rv) @ Rg
    return out


def spread_deg(
    orientations: np.ndarray,
    grain_ids: np.ndarray,
    grain_orientations: Dict[int, np.ndarray],
    sym: np.ndarray,
) -> Dict[int, float]:
    """Measured spread per grain: the largest misorientation (with the crystal's symmetry
    operators ``sym``) between a voxel and its grain's centre orientation, in degrees."""
    from .ray_indexing import misorientation_deg

    ids = np.asarray(grain_ids)
    res: Dict[int, float] = {}
    for g in sorted(int(x) for x in np.unique(ids)):
        sel = np.nonzero(ids == g)[0]
        Rg = np.broadcast_to(np.asarray(grain_orientations[g]), (len(sel), 3, 3))
        res[g] = float(np.max(misorientation_deg(orientations[sel], Rg, sym)))
    return res
