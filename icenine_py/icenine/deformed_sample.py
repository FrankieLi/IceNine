"""Samples whose grains hold random voxel-to-voxel misorientations (AR6 grainopt, ideas I11, I12).

A plastically deformed grain (geometrically necessary dislocations) holds orientation variations
that are not a smooth low-order field. ``random_field_orientations`` builds such a layer from a
one-orientation-per-grain truth: each voxel's orientation is a small rotation left-multiplied onto
its grain's orientation, ``R_v = exp([r_v]x) R_grain``, with the rotation vectors ``r_v`` a random
field over the grain:

* ``corr_len = 0``: independent Gaussian vectors per voxel (white noise; neighbours unrelated).
* ``corr_len > 0``: white noise smoothed by a Gaussian kernel of that length (same units as the
  centres) over the grain's voxels, giving spatially correlated variations.

In both cases the field's mean over the grain is removed (the grain's mean orientation stays the
grain orientation to first order) and the field is scaled so its RMS rotation angle over the grain
equals the grain's level (degrees). The model is a choice, not a measurement. Only numpy is used.
"""

from typing import Dict, Tuple

import numpy as np

from .gradient_sample import rotation_matrices

LEVELS_DEG = (0.5, 1.0, 2.0, 5.0)


def random_rotvec_field(
    xy: np.ndarray, level_deg: float, corr_len: float, rng: np.random.Generator
) -> np.ndarray:
    """Rotation vectors (m, 3) (radians) over one grain's voxel centres ``xy`` (m, 2): zero mean,
    RMS angle ``level_deg``; white (``corr_len`` 0) or Gaussian-smoothed with length ``corr_len``."""
    m = len(xy)
    w = rng.standard_normal((m, 3))
    if m < 2 or level_deg == 0.0:
        return np.zeros((m, 3))
    if corr_len > 0:
        d2 = ((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1)
        K = np.exp(-0.5 * d2 / corr_len**2)
        w = K @ w
    w = w - w.mean(0)
    rms = float(np.sqrt((w * w).sum(1).mean()))
    if rms <= 0:
        return np.zeros((m, 3))
    return w * (np.radians(level_deg) / rms)


def random_field_orientations(
    centres: np.ndarray,
    grain_ids: np.ndarray,
    grain_orientations: Dict[int, np.ndarray],
    levels_deg: Dict[int, float],
    corr_len: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Orientation (n, 3, 3) of every voxel: ``exp([r_v]x) R_grain`` with ``r_v`` drawn by
    ``random_rotvec_field`` per grain, in increasing grain-id order (one draw per grain, also for a
    level-0 grain, so the draws do not depend on the levels)."""
    xy = np.asarray(centres, dtype=np.float64)[:, :2]
    ids = np.asarray(grain_ids)
    out = np.empty((len(ids), 3, 3), dtype=np.float64)
    for g in sorted(int(x) for x in np.unique(ids)):
        sel = np.nonzero(ids == g)[0]
        rv = random_rotvec_field(xy[sel], float(levels_deg[g]), corr_len, rng)
        Rg = np.asarray(grain_orientations[g], dtype=np.float64)
        out[sel] = rotation_matrices(rv) @ Rg if np.any(rv) else Rg
    return out


def neighbour_misorientation_deg(
    R: np.ndarray, neighbours: "list[list[int]]", sym: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """(pairs (k, 2), misorientation deg (k,)) over unique neighbour pairs."""
    from .ray_indexing import misorientation_deg

    pairs = np.array(
        [(v, u) for v in range(len(neighbours)) for u in neighbours[v] if u > v], dtype=np.int64
    ).reshape(-1, 2)
    if not len(pairs):
        return pairs, np.zeros(0)
    return pairs, np.asarray(misorientation_deg(R[pairs[:, 0]], R[pairs[:, 1]], sym), float)
