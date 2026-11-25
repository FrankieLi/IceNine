"""tex1: the pure parts of the textured sample ``tex24`` (numpy and scipy only, no file access).

``tex1_make_sample.py`` and ``tex1_sample_check.py`` use these; ``tests/test_tex1.py`` tests them.

Terms. Class of a grain: 0 sharp, 1 noise only, 2 gradient only, 3 both. Edge pair: two voxels
sharing an edge; touching pair: two voxels sharing an edge or a corner. ``R_g``: the one
orientation drawn per grain.

Gradient field of a gradient grain, ``u = (x - c) / rho`` (``c`` the centroid of the grain's voxel
centres, ``rho`` the RMS distance of the centres from it), rotation vector in radians
``r(u) = s_g (M1 u + CURVATURE (M2 (q(u) - mean q)))`` with ``q = (ux^2, ux uy, uy^2)``, ``M1``
3 by 2 and ``M2`` 3 by 3 standard normal; ``R = exp(r) R_g``. ``s_g`` is set so that the largest
edge-pair misorientation inside the grain equals the grain's strength ``t_g``.
Noise: per voxel an independent rotation vector uniform in the ball of ``NOISE_BALL_DEG``
degrees, ``R = exp(n_v) exp(r) R_g`` (``r = 0`` in a grain with no gradient): the noise parts of
any two voxels differ by at most twice the radius.
"""

from typing import Any, Dict, List, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

NOISE_BALL_DEG = 1.0
GRADIENT_CAP_DEG = 0.5
GRADIENT_MIN_DEG = 0.1
CURVATURE = 0.25
SHARP, NOISE_ONLY, GRADIENT_ONLY, BOTH = 0, 1, 2, 3
CLASS_NAMES = ("sharp", "noise", "gradient", "both")
FILE_ROUNDING_DEG = 5.0e-4  # a .mic row keeps phi1 to 1e-4 degree (``MicFile.write``)


# ----------------------------------------------------------------------------- neighbour pairs


def edge_pairs(centres: np.ndarray, side: float) -> np.ndarray:
    """Pairs (m, 2) of triangles of one side length that share an edge. Edge neighbours have
    their centres side/sqrt(3) apart (0.577 side); the next distances are 1.0 side (same
    direction) and 1.155 side (vertex only), so a cut at 0.7 side selects exactly the edge
    neighbours."""
    tree = cKDTree(np.asarray(centres, dtype=np.float64))
    p = tree.query_pairs(0.7 * side, output_type="ndarray")
    return np.asarray(p, dtype=np.int64).reshape(-1, 2)


def touching_pairs(centres: np.ndarray, side: float) -> np.ndarray:
    """Pairs (m, 2) of triangles that share an edge or a corner: centre distances 0.577, 1.0 and
    1.155 side (three of the 12 around a corner each); the next distance is 1.528 side, so a cut
    at 1.2 side selects exactly them."""
    tree = cKDTree(np.asarray(centres, dtype=np.float64))
    p = tree.query_pairs(1.2 * side, output_type="ndarray")
    return np.asarray(p, dtype=np.int64).reshape(-1, 2)


def relative_angle_deg(Ma: np.ndarray, Mb: np.ndarray) -> np.ndarray:
    """Rotation angle (degrees) of Ma Mb^T, for stacks (n, 3, 3); accurate for tiny angles
    (arctan2 of the skew part and the trace part)."""
    R = np.einsum("nij,nkj->nik", Ma, Mb)
    skew = 0.5 * np.stack(
        [R[:, 2, 1] - R[:, 1, 2], R[:, 0, 2] - R[:, 2, 0], R[:, 1, 0] - R[:, 0, 1]], axis=1
    )
    sin_t = np.linalg.norm(skew, axis=1)
    cos_t = 0.5 * (np.trace(R, axis1=1, axis2=2) - 1.0)
    return np.degrees(np.arctan2(sin_t, cos_t))


def rotation_matrices(rotvec: np.ndarray) -> np.ndarray:
    """Rotation matrices (n, 3, 3) of rotation vectors (n, 3) in radians (the library's)."""
    from icenine.gradient_sample import rotation_matrices as rm

    return rm(np.asarray(rotvec, dtype=np.float64).reshape(-1, 3))


# ----------------------------------------------------------------------------- the classes


def draw_classes(
    grain_ids: Sequence[int], seed: int, noise_fraction: float, gradient_fraction: float
) -> Tuple[List[int], List[int]]:
    """The noise grains (``round(noise_fraction * n)`` of them, ``default_rng([seed, 1])``) and the
    gradient grains (``default_rng([seed, 2])``), drawn independently from the sorted ids."""
    ids = np.array(sorted(int(g) for g in grain_ids))
    n = len(ids)
    noise = np.random.default_rng([seed, 1]).choice(
        ids, size=int(round(noise_fraction * n)), replace=False
    )
    grad = np.random.default_rng([seed, 2]).choice(
        ids, size=int(round(gradient_fraction * n)), replace=False
    )
    return sorted(int(g) for g in noise), sorted(int(g) for g in grad)


def class_of(is_noise: bool, is_gradient: bool) -> int:
    return (1 if is_noise else 0) + (2 if is_gradient else 0)


# ----------------------------------------------------------------------------- the gradient field


def field_vectors(xy: np.ndarray, M1: np.ndarray, M2: np.ndarray) -> np.ndarray:
    """Unit-strength rotation vectors (n, 3) (radians) of the field at the grain's voxel centres
    ``xy`` (n, 2): ``M1 u + CURVATURE M2 (q - mean q)``, ``u = (x - c) / rho``."""
    xy = np.asarray(xy, dtype=np.float64)[:, :2]
    d = xy - xy.mean(axis=0)
    rho = float(np.sqrt((d**2).sum(axis=1).mean()))
    u = d / rho
    q = np.stack([u[:, 0] ** 2, u[:, 0] * u[:, 1], u[:, 1] ** 2], axis=1)
    return u @ np.asarray(M1).T + CURVATURE * ((q - q.mean(axis=0)) @ np.asarray(M2).T)


def scale_to_edge_max(
    rv: np.ndarray, la: np.ndarray, lb: np.ndarray, target_deg: float
) -> Tuple[float, float]:
    """The factor ``s`` such that the largest misorientation over the pairs (la, lb) of the
    rotations ``exp(s rv)`` equals ``target_deg`` (iterated: the angle is linear in ``s`` only to
    first order), and that largest misorientation (never above the target). ``rv`` radians."""
    M = rotation_matrices(rv)
    cur = float(relative_angle_deg(M[la], M[lb]).max())
    scale = 1.0
    for _ in range(60):
        scale *= target_deg / cur
        M = rotation_matrices(rv * scale)
        cur = float(relative_angle_deg(M[la], M[lb]).max())
        if abs(cur - target_deg) <= 1e-12 * target_deg:
            break
    if cur > target_deg:
        scale *= target_deg / cur * (1.0 - 1e-12)
        M = rotation_matrices(rv * scale)
        cur = float(relative_angle_deg(M[la], M[lb]).max())
    return float(scale), cur


def ball_vectors(n: int, rng: np.random.Generator, radius_deg: float) -> np.ndarray:
    """``n`` rotation vectors (radians) uniform in the ball of ``radius_deg``."""
    g = rng.standard_normal((n, 3))
    g /= np.linalg.norm(g, axis=1, keepdims=True)
    r = radius_deg * np.cbrt(rng.random(n))
    return g * np.radians(r)[:, None]


def textured_orientations(
    centres: np.ndarray,
    grain_ids: np.ndarray,
    grain_orientations: Dict[int, np.ndarray],
    noise_grains: Sequence[int],
    gradient_grains: Sequence[int],
    pairs: np.ndarray,
    seed: int,
    field_params: Dict[int, Dict[str, Any]] = None,  # type: ignore[assignment]
) -> Tuple[np.ndarray, np.ndarray, List[Dict[str, Any]]]:
    """The textured truth. Returns (R (n, 3, 3), R_smooth (n, 3, 3) the same without noise,
    records of the gradient grains).

    Draws: gradient grains in increasing id order from ``default_rng([seed, 4])``, per grain
    ``M1`` (3 by 2), ``M2`` (3 by 3), strength ``t_g`` uniform on [GRADIENT_MIN_DEG,
    GRADIENT_CAP_DEG); noise from ``default_rng([seed, 3])``, per noise grain in increasing id
    order, one ball vector per voxel in row order. With ``field_params`` (id -> M1, M2, scale) the
    gradient is rebuilt from those instead of drawn (the check script's path).
    """
    xy = np.asarray(centres, dtype=np.float64)[:, :2]
    ids = np.asarray(grain_ids)
    n = len(ids)
    R = np.empty((n, 3, 3))
    R_smooth = np.empty((n, 3, 3))
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    same = ids[pairs[:, 0]] == ids[pairs[:, 1]] if len(pairs) else np.zeros(0, dtype=bool)
    pairs = pairs[same]
    pair_grain = ids[pairs[:, 0]] if len(pairs) else np.zeros(0, dtype=ids.dtype)
    grng = np.random.default_rng([seed, 4])
    nrng = np.random.default_rng([seed, 3])
    gset, nset = set(int(g) for g in gradient_grains), set(int(g) for g in noise_grains)
    records: List[Dict[str, Any]] = []
    for g in sorted(int(x) for x in np.unique(ids)):
        sel = np.nonzero(ids == g)[0]
        Rg = np.asarray(grain_orientations[g], dtype=np.float64)
        Rs = np.broadcast_to(Rg, (len(sel), 3, 3)).copy()
        if g in gset:
            if field_params is None:
                M1 = grng.standard_normal((3, 2))
                M2 = grng.standard_normal((3, 3))
                t_g = float(grng.uniform(GRADIENT_MIN_DEG, GRADIENT_CAP_DEG))
            else:
                M1, M2 = np.asarray(field_params[g]["M1"]), np.asarray(field_params[g]["M2"])
                t_g = float(field_params[g]["t_g_deg"])
            rv = field_vectors(xy[sel], M1, M2)
            gp = pairs[pair_grain == g]
            rec: Dict[str, Any] = dict(
                id=g,
                n_voxels=int(len(sel)),
                M1=M1.tolist(),
                M2=M2.tolist(),
                t_g_deg=t_g,
                n_edge_pairs=int(len(gp)),
            )
            if field_params is not None:
                scale = float(field_params[g]["scale"])
                edge_max = None
            elif len(gp):
                local = {int(v): i for i, v in enumerate(sel)}
                la = np.array([local[int(a)] for a in gp[:, 0]])
                lb = np.array([local[int(b)] for b in gp[:, 1]])
                scale, edge_max = scale_to_edge_max(rv, la, lb, t_g)
            else:  # no edge pair inside the grain: no gradient (reported)
                scale, edge_max = 0.0, 0.0
            rv = rv * scale
            Rs = rotation_matrices(rv) @ Rg
            rec.update(
                scale=scale,
                edge_max_deg=edge_max,
                spread_to_centre_deg=float(np.degrees(np.linalg.norm(rv, axis=1).max())),
            )
            records.append(rec)
        R_smooth[sel] = Rs
        if g in nset:
            nv = ball_vectors(len(sel), nrng, NOISE_BALL_DEG)
            R[sel] = rotation_matrices(nv) @ Rs
        else:
            R[sel] = Rs
    return R, R_smooth, records


# ----------------------------------------------------------------------------- statistics


def percentiles(values: np.ndarray, qs: Sequence[float] = (50, 90, 99, 99.9)) -> Dict[str, float]:
    """Percentiles of ``values`` as {"p50": ..}."""
    v = np.asarray(values, dtype=np.float64)
    if len(v) == 0:
        return {f"p{q:g}": float("nan") for q in qs}
    return {f"p{q:g}": float(np.percentile(v, q)) for q in qs}


def summary(values: np.ndarray) -> Dict[str, Any]:
    """n, max, median, p99 of ``values``."""
    v = np.asarray(values, dtype=np.float64)
    if len(v) == 0:
        return dict(n=0, max=None, median=None, p99=None)
    return dict(
        n=int(len(v)),
        max=float(v.max()),
        median=float(np.median(v)),
        p99=float(np.percentile(v, 99)),
    )


def nearest_centre_rows(left_xy: np.ndarray, n: int) -> np.ndarray:
    """Sorted row indices of the ``n`` voxels nearest the layer centre (0, 0), by the left
    vertex (x, y) of the .mic row: the expression of ``render_full.py --n-voxels``
    (``np.argsort`` of the distance, first ``n``)."""
    d = np.linalg.norm(np.asarray(left_xy, dtype=np.float64)[:, :2] - np.array([0.0, 0.0]), axis=1)
    return np.sort(np.argsort(d)[:n])


def gate(checks: Dict[str, bool]) -> Dict[str, Any]:
    """The hard gate: every named check passes."""
    return dict(checks=dict(checks), passed=bool(all(checks.values())))
