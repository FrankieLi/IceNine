"""Cubic coincidence-site-lattice (CSL) tools: the Sigma table, relatives of an orientation,
Brandon classification, and the symmetry-reduced misorientation (numpy/scipy only).

Conventions: an orientation R maps crystal to sample coordinates; the cubic crystal symmetry acts
on the right (R ~ R S). The crystal-frame misorientation of two orientations is M = Ra^T Rb, and
Rb is a Sigma-n relative of Ra when S1 M S2 equals the ideal CSL rotation C0 of Sigma n for some
cubic S1, S2. Relatives of g are {g S1 C0}, deduplicated modulo the right symmetry.
"""

from typing import Dict, List, NamedTuple, Optional, Tuple

import numpy as np
from scipy.spatial.transform import Rotation

BRANDON_DEG = 15.0
FAMILIES = {"<100>": (1, 0, 0), "<110>": (1, 1, 0), "<111>": (1, 1, 1)}


def cubic_ops() -> np.ndarray:
    """The 24 proper rotations of the cubic point group, (24, 3, 3)."""
    gens = [
        Rotation.from_rotvec([0, 0, np.pi / 2]).as_matrix(),
        Rotation.from_rotvec(np.array([1, 1, 1]) * (2 * np.pi / 3) / np.sqrt(3)).as_matrix(),
    ]
    ops = [np.eye(3)]
    changed = True
    while changed:
        changed = False
        for g in gens:
            for o in list(ops):
                n = g @ o
                if not any(np.allclose(n, m, atol=1e-9) for m in ops):
                    ops.append(n)
                    changed = True
    out = np.stack(ops)
    assert out.shape == (24, 3, 3)
    return out


OPS = cubic_ops()


def reduced_misorientation_deg(Ra: np.ndarray, Rb: np.ndarray) -> np.ndarray:
    """Cubic-symmetry-reduced angle (deg) between orientations, min over S of angle(Ra S Rb^T)...
    Batched over leading axes (broadcasting)."""
    Ra = np.asarray(Ra, dtype=np.float64)
    Rb = np.asarray(Rb, dtype=np.float64)
    M = (Ra[..., None, :, :] @ OPS) @ np.swapaxes(Rb, -1, -2)[..., None, :, :]
    ang = np.degrees(Rotation.from_matrix(M.reshape(-1, 3, 3)).magnitude()).reshape(M.shape[:-2])
    return ang.min(axis=-1)


class CSLEntry(NamedTuple):
    sigma: int
    label: str  # e.g. "3", "13a"
    angle_deg: float
    axis: Tuple[int, int, int]
    C0: np.ndarray  # ideal rotation, crystal frame


_TABLE: Optional[List[CSLEntry]] = None


def _reduced_rotation_angle(C: np.ndarray) -> float:
    both = OPS[:, None] @ C @ OPS[None]
    return float(np.degrees(Rotation.from_matrix(both.reshape(-1, 3, 3)).magnitude().min()))


def csl_table(max_sigma: int = 29) -> List[CSLEntry]:
    """Cubic CSL rotations with odd Sigma <= max_sigma, one entry per distinct (Sigma, reduced
    angle): for an axis [uvw] and integer m, tan(theta/2) = n sqrt(N)/m with N = u^2+v^2+w^2, gcd(m, n) = 1, gives
    Sigma = m^2 + n^2 N (halved while even). Sorted by Sigma; variants of one Sigma get a, b, ..."""
    global _TABLE
    if _TABLE is None or max(e.sigma for e in _TABLE) < max_sigma:
        seen: Dict[Tuple[int, float], Tuple[int, float, Tuple[int, int, int]]] = {}
        for u in range(5):
            for v in range(u + 1):
                for w in range(v + 1):
                    if (u, v, w) == (0, 0, 0) or np.gcd.reduce([u, v, w]) != 1:
                        continue
                    N = u * u + v * v + w * w
                    for m, n in [(m, n) for n in range(1, 5) for m in range(1, 120)]:
                        if np.gcd(m, n) != 1:
                            continue
                        sig = m * m + n * n * N
                        while sig % 2 == 0:
                            sig //= 2
                        if sig == 1 or sig > max_sigma:
                            continue
                        th = 2 * np.degrees(np.arctan(n * np.sqrt(N) / m))
                        C = Rotation.from_rotvec(
                            np.radians(th) * np.array([u, v, w]) / np.sqrt(N)
                        ).as_matrix()
                        key = (sig, round(_reduced_rotation_angle(C), 2))
                        if key not in seen:
                            seen[key] = (sig, float(th), (u, v, w))
        rows = sorted(seen.values())
        count: Dict[int, int] = {}
        for sig, _, _ in rows:
            count[sig] = count.get(sig, 0) + 1
        idx: Dict[int, int] = {}
        out: List[CSLEntry] = []
        for sig, th, ax in rows:
            idx[sig] = idx.get(sig, 0) + 1
            label = f"{sig}" if count[sig] == 1 else f"{sig}{'abcdefgh'[idx[sig] - 1]}"
            C0 = Rotation.from_rotvec(
                np.radians(th) * np.array(ax) / np.linalg.norm(ax)
            ).as_matrix()
            out.append(CSLEntry(sig, label, th, ax, C0))
        _TABLE = out
    return [e for e in _TABLE if e.sigma <= max_sigma]


def dedup_orientations(R: np.ndarray, tol_deg: float = 0.01) -> np.ndarray:
    """Indices of the first occurrence of each distinct orientation modulo the cubic symmetry."""
    keep: List[int] = []
    for i in range(len(R)):
        if not keep or reduced_misorientation_deg(R[i], R[keep]).min() > tol_deg:
            keep.append(i)
    return np.array(keep, dtype=np.int64)


def csl_relatives(
    R: np.ndarray, sigmas: Optional[List[int]] = None, max_sigma: int = 29
) -> Tuple[np.ndarray, List[str]]:
    """All distinct Sigma relatives of orientation R (3,3): {R S C0} modulo right symmetry, for
    every table entry with Sigma in `sigmas` (default: all odd Sigma <= max_sigma). Returns
    (n, 3, 3) matrices and their Sigma labels."""
    R = np.asarray(R, dtype=np.float64)
    mats: List[np.ndarray] = []
    labels: List[str] = []
    for e in csl_table(max_sigma):
        if sigmas is not None and e.sigma not in sigmas:
            continue
        cand = R[None] @ OPS @ e.C0  # (24, 3, 3)
        keep = dedup_orientations(cand)
        mats.extend(cand[keep])
        labels.extend([e.label] * len(keep))
    return np.stack(mats), labels


def csl_classify(R_true: np.ndarray, R_est: np.ndarray, max_sigma: int = 29) -> Dict[str, object]:
    """Brandon classification of the crystal-frame misorientation M = R_true^T R_est: the angle
    (reduced), and the lowest Sigma <= max_sigma whose ideal rotation is within
    15 deg / sqrt(Sigma) of M over the cubic symmetry on both sides (0 = none, 1 = identity within
    1 deg)."""
    M = np.asarray(R_true).T @ np.asarray(R_est)
    both = OPS[:, None] @ M @ OPS[None]
    ang = float(np.degrees(Rotation.from_matrix(both.reshape(-1, 3, 3)).magnitude().min()))
    if ang < 1.0:
        return dict(sigma=1, label="1", angle=ang, deviation=ang)
    flat = both.reshape(-1, 3, 3)
    for e in csl_table(max_sigma):
        dev = float(np.degrees(Rotation.from_matrix(flat @ e.C0.T).magnitude().min()))
        if dev <= BRANDON_DEG / np.sqrt(e.sigma):
            return dict(sigma=e.sigma, label=e.label, angle=ang, deviation=dev)
    return dict(sigma=0, label="none", angle=ang, deviation=float("nan"))


def invariant_reflection_mask(
    hkls: np.ndarray, R: np.ndarray, relatives: np.ndarray, tol: float = 1e-3
) -> np.ndarray:
    """Which reflections (crystal-frame directions hkl, (n, 3), any lattice) of orientation R point
    to the same sample direction under some relative R' (R' hkl' = R hkl for a symmetry-equivalent
    hkl' of the same family): reflections the truth shares with its CSL relatives. Returns a bool
    (n_rel, n) array, one row per relative."""
    h = np.asarray(hkls, dtype=np.float64)
    n = np.linalg.norm(h, axis=1, keepdims=True)
    u = h / n
    out = np.zeros((len(relatives), len(h)), dtype=bool)
    sample_dirs = (R @ u.T).T  # (n, 3)
    fam = np.einsum("sij,nj->nsi", OPS, u)  # (n, 24, 3) symmetry-equivalent directions
    for k, Rk in enumerate(relatives):
        rel_dirs = np.einsum("ij,nsj->nsi", Rk, fam)  # (n, 24, 3)
        d = np.linalg.norm(rel_dirs - sample_dirs[:, None, :], axis=-1).min(axis=1)
        out[k] = d < tol
    return out
