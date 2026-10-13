"""Pure helpers for the Phase D sample: grain grouping, orientation draws, separation checks.

Everything here is symmetry-agnostic: the symmetry operators come in as an (N_sym, 4) array of
proper-rotation quaternions [w, x, y, z] (see icenine.orientation_search.get_symmetry_quaternions),
and misorientation is the minimum over those operators.
"""

from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation

from icenine.geometry import euler_to_matrix, matrix_to_euler
from icenine.sampling import matrix_to_quaternion, quaternion_to_matrix, reduce_to_fundamental_zone

EULER_DECIMALS = 6  # decimals written to the .mic (degrees)


def group_grains(euler_deg: np.ndarray, decimals: int = 6) -> Tuple[np.ndarray, np.ndarray]:
    """Group voxels with identical orientation. Returns (grain_id per voxel, unique eulers).

    Grain ids are numbered by first appearance in the file (stable and deterministic).
    """
    key = np.round(np.asarray(euler_deg, dtype=np.float64), decimals)
    uniq, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = np.asarray(inv).reshape(-1)
    order = np.argsort(first)  # unique-row index -> rank of first appearance
    rank = np.empty_like(order)
    rank[order] = np.arange(len(order))
    return rank[inv].astype(np.int64), uniq[order]


def euler_to_quat(euler_deg: np.ndarray) -> np.ndarray:
    """(N, 3) Bunge ZXZ Euler angles in degrees -> (N, 4) quaternions [w, x, y, z]."""
    out = np.empty((len(euler_deg), 4))
    for i, e in enumerate(np.asarray(euler_deg, dtype=np.float64)):
        out[i] = matrix_to_quaternion(euler_to_matrix(*e).astype(np.float64))
    return out


def quat_to_euler(q: np.ndarray) -> np.ndarray:
    """(N, 4) quaternions [w, x, y, z] -> (N, 3) Bunge ZXZ Euler angles in degrees."""
    return np.array([matrix_to_euler(quaternion_to_matrix(x)) for x in q], dtype=np.float64)


def _qmul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Hamilton product, broadcasting over leading dims; last dim is [w, x, y, z]."""
    w1, x1, y1, z1 = np.moveaxis(a, -1, 0)
    w2, x2, y2, z2 = np.moveaxis(b, -1, 0)
    return np.stack(
        [
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
            w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2,
        ],
        axis=-1,
    )


def misorientation_matrix_deg(qa: np.ndarray, qb: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """Symmetry-reduced misorientation (deg) between every row of qa (A,4) and qb (B,4)."""
    conj = qa * np.array([1.0, -1.0, -1.0, -1.0])
    rel = _qmul(conj[:, None, :], qb[None, :, :])  # (A, B, 4)
    prod = _qmul(rel[:, :, None, :], np.asarray(sym)[None, None, :, :])  # (A, B, S, 4)
    w = np.abs(prod[..., 0]).max(axis=-1)
    return np.degrees(2.0 * np.arccos(np.clip(w, 0.0, 1.0)))


def draw_uniform_quats(n: int, rng: np.random.Generator) -> np.ndarray:
    """n quaternions uniform on SO(3) (normalised 4-D Gaussians), canonical sign w >= 0."""
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    q[q[:, 0] < 0] *= -1.0
    return q


def reduce_quats(q: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """Reduce each quaternion to the fundamental zone (largest |w| over the symmetry group)."""
    return np.array([reduce_to_fundamental_zone(x, sym) for x in q])


def draw_orientation(seed: int, grain: int, attempt: int, sym: np.ndarray) -> np.ndarray:
    """One FZ-reduced uniform orientation, deterministic in (seed, grain, attempt), as the
    quaternion of the Euler angles rounded to the .mic precision."""
    rng = np.random.default_rng([seed, grain, attempt])
    q = reduce_quats(draw_uniform_quats(1, rng), sym)
    eul = np.round(quat_to_euler(q), EULER_DECIMALS)
    return euler_to_quat(eul)[0]


def grain_adjacency(
    positions: np.ndarray, grain: np.ndarray, radius: float
) -> List[Tuple[int, int]]:
    """Unordered grain pairs (a < b) with at least one voxel pair closer than `radius`."""
    tree = cKDTree(positions)
    pairs = tree.query_pairs(radius, output_type="ndarray")
    if len(pairs) == 0:
        return []
    ga, gb = grain[pairs[:, 0]], grain[pairs[:, 1]]
    keep = ga != gb
    lo, hi = np.minimum(ga[keep], gb[keep]), np.maximum(ga[keep], gb[keep])
    return sorted({(int(a), int(b)) for a, b in zip(lo, hi)})


def draw_separated(
    n_grains: int,
    old_q: np.ndarray,
    adjacency: Sequence[Tuple[int, int]],
    sym: np.ndarray,
    seed: int = 0,
    min_sep_deg: float = 1.0,
    max_rounds: int = 200,
) -> Tuple[np.ndarray, Dict[str, float]]:
    """One new orientation per grain, with no new orientation within `min_sep_deg` of any old
    orientation or of an adjacent grain's new one (symmetry-reduced). Violators are redrawn
    (attempt counter +1, deterministic in seed/grain/attempt). Returns (new_q, stats)."""
    attempts = np.zeros(n_grains, dtype=np.int64)
    new_q = np.array([draw_orientation(seed, g, 0, sym) for g in range(n_grains)])
    adj = np.asarray(adjacency, dtype=np.int64).reshape(-1, 2)
    for _ in range(max_rounds):
        bad = set()
        d_old = misorientation_matrix_deg(new_q, old_q, sym).min(axis=1)
        bad.update(np.nonzero(d_old < min_sep_deg)[0].tolist())
        if len(adj):
            d_adj = np.array(
                [
                    misorientation_matrix_deg(new_q[a : a + 1], new_q[b : b + 1], sym)[0, 0]
                    for a, b in adj
                ]
            )
            for a, b in adj[d_adj < min_sep_deg]:
                bad.add(int(max(a, b)))  # redraw the later grain of the pair
        if not bad:
            break
        for g in sorted(bad):
            attempts[g] += 1
            new_q[g] = draw_orientation(seed, g, int(attempts[g]), sym)
    else:
        raise RuntimeError("draw_separated did not converge")
    d_old = misorientation_matrix_deg(new_q, old_q, sym).min(axis=1)
    d_adj = (
        np.array(
            [
                misorientation_matrix_deg(new_q[a : a + 1], new_q[b : b + 1], sym)[0, 0]
                for a, b in adj
            ]
        )
        if len(adj)
        else np.array([np.inf])
    )
    stats = {
        "n_redrawn": int((attempts > 0).sum()),
        "max_attempts": int(attempts.max()),
        "min_new_vs_old_deg": float(d_old.min()),
        "min_neighbour_new_deg": float(d_adj.min()),
        "n_adjacent_pairs": int(len(adj)),
    }
    return new_q, stats


def write_mic_with_euler(src: str, dst: str, grain: np.ndarray, new_euler_deg: np.ndarray) -> None:
    """Copy `src` .mic replacing the Euler columns (7-9) by new_euler_deg[grain]; every other
    token (positions, direction, generation, phase, confidence) is copied verbatim."""
    with open(src) as f:
        lines = f.read().splitlines()
    out = [lines[0]]
    body = [ln for ln in lines[1:] if ln.strip()]
    assert len(body) == len(grain)
    for ln, g in zip(body, grain):
        t = ln.split()
        e = new_euler_deg[g]
        t[6], t[7], t[8] = (
            f"{e[0]:.{EULER_DECIMALS}f}",
            f"{e[1]:.{EULER_DECIMALS}f}",
            f"{e[2]:.{EULER_DECIMALS}f}",
        )
        out.append("\t".join(t))
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


def so3_geodesic_deg(ra: np.ndarray, rb: np.ndarray) -> float:
    """Plain (no symmetry) angle between two rotation matrices, degrees."""
    return float(np.degrees(np.linalg.norm(Rotation.from_matrix(ra @ rb.T).as_rotvec())))
