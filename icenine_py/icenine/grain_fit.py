"""Full-grain orientation fitting: fit every voxel of a grain together against the images.

Plan: ``docs/grain_opt/findings.md`` (and the plan doc linked there). The pieces:

* **Grouping** (``group_grains``): flood fill over voxel neighbours whose misorientation is below
  a threshold (2 deg); grains are labelled largest first.
* **Spot model** (``SpotModel``): the forward model's peak geometry for many voxels at once, each
  with its own orientation and triangle (the math of ``SmoothCost.predict``, vectorised over
  voxels; differentiable in the orientations).
* **Coverage points** (``coverage_points``): each projected voxel triangle is cut into ``n_sub**2``
  equal sub-triangles; their centroids carry the triangle's projected area (pixels^2) split
  equally. The projection along a fixed ray onto a plane is affine, so the sub-triangle centroids
  are exact barycentric combinations of the projected vertices. Mass per pixel is then the
  geometric coverage of that pixel by the voxel's spot (in the limit of many sub-triangles),
  which is what the stock forward model paints.
* **Region** (``GrainRegion``): per reflection of the grain (key: reflection, omega branch,
  detector), a padded window in pixels and frames around the spot at the start orientation, the
  observed lit pixels in it (binary, and their clean intensity when a store is given), and the
  fixed prediction of the shell voxels (neighbours of the grain) that falls in the window.
* **Objectives**: per-reflection Sinkhorn divergence between the predicted and the observed spot
  as unit-mass distributions over (col, row, frame) (``sinkhorn_divergence``; main), and least
  squares with one scale per reflection on a soft raster of the prediction (``ls_objective``;
  baseline arm).
* **Fields** (``FieldA``): grain rotation plus a linear (optionally quadratic) rotation-vector
  field in position, ``R_v = exp([w0 + J d_v (+ Q q(d_v))]) R_ref``; ``d_v`` is the voxel's offset
  from the grain centre divided by the grain radius, so every coefficient is the rotation (rad) at
  the grain's edge. Stage B adds a free rotation vector per voxel on top of stage A.

The strict binary cost (``VoxelCostFunction``) is the referee and is not computed here; scripts
under ``scripts/grain_opt/`` call it per voxel.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch

from .ray_indexing import misorientation_deg
from .diffraction_core import BraggConstants, scattering_omegas_fixed_beam

Tensor = torch.Tensor
DT = torch.float64


# --------------------------------------------------------------------------------------------
# rotations
# --------------------------------------------------------------------------------------------
def rotvec_to_matrix(w: Tensor) -> Tensor:
    """Rotation matrices (..., 3, 3) of rotation vectors (..., 3) (Rodrigues; smooth at 0)."""
    th2 = (w * w).sum(-1, keepdim=True)[..., None]  # (..., 1, 1)
    small = th2 < 1e-12
    th2s = torch.where(small, torch.ones_like(th2), th2)
    th = torch.sqrt(th2s)
    a = torch.where(small, 1.0 - th2 / 6.0, torch.sin(th) / th)
    b = torch.where(small, 0.5 - th2 / 24.0, (1.0 - torch.cos(th)) / th2s)
    z = torch.zeros_like(w[..., 0])
    K = torch.stack(
        [
            torch.stack([z, -w[..., 2], w[..., 1]], -1),
            torch.stack([w[..., 2], z, -w[..., 0]], -1),
            torch.stack([-w[..., 1], w[..., 0], z], -1),
        ],
        -2,
    )
    eye = torch.eye(3, dtype=w.dtype).expand(K.shape)
    return eye + a * K + b * (K @ K)


def matrix_to_rotvec(R: np.ndarray) -> np.ndarray:
    """Rotation vectors (n, 3) of rotation matrices (n, 3, 3) (numpy, via scipy)."""
    from scipy.spatial.transform import Rotation

    return Rotation.from_matrix(np.asarray(R, dtype=np.float64)).as_rotvec()


def align_to(R: np.ndarray, Rref: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """The symmetry equivalents ``R S`` (n, 3, 3) of ``R`` closest to ``Rref`` (3, 3)."""
    R = np.asarray(R, dtype=np.float64).reshape(-1, 3, 3)
    M = np.einsum("ji,njk->nik", Rref, R)  # Rref^T R
    tr = M.reshape(-1, 9) @ sym.transpose(0, 2, 1).reshape(len(sym), 9).T
    s = np.argmax(tr, axis=1)
    return np.einsum("nij,njk->nik", R, sym[s])


def align_ops(R: np.ndarray, Rref: np.ndarray, sym: np.ndarray) -> np.ndarray:
    """The symmetry operators S (n, 3, 3) with ``R S`` closest to ``Rref`` (as ``align_to``)."""
    R = np.asarray(R, dtype=np.float64).reshape(-1, 3, 3)
    M = np.einsum("ji,njk->nik", Rref, R)
    tr = M.reshape(-1, 9) @ sym.transpose(0, 2, 1).reshape(len(sym), 9).T
    return sym[np.argmax(tr, axis=1)]


def rotate_strain(strain: np.ndarray, S: np.ndarray) -> np.ndarray:
    """Crystal-frame strain (n, 6) [e11 e22 e33 e23 e13 e12] re-expressed for the symmetry-
    equivalent orientation R S: eps' = S^T eps S (then g = (R S) M' (S^T g0) = R M g0)."""
    e = np.asarray(strain, np.float64).reshape(-1, 6)
    E = np.stack([e[:, 0], e[:, 5], e[:, 4], e[:, 5], e[:, 1], e[:, 3], e[:, 4], e[:, 3], e[:, 2]],
                 1).reshape(-1, 3, 3)
    E2 = np.einsum("nji,njk,nkl->nil", S, E, S)
    return np.stack([E2[:, 0, 0], E2[:, 1, 1], E2[:, 2, 2], E2[:, 1, 2], E2[:, 0, 2], E2[:, 0, 1]], 1)


def mean_orientation(R: np.ndarray, sym: np.ndarray, Rref: Optional[np.ndarray] = None) -> np.ndarray:
    """Mean orientation of ``R`` (n, 3, 3): symmetry-aligned to ``Rref`` (default the first), then
    the projection of the mean matrix onto SO(3)."""
    Rref = R[0] if Rref is None else Rref
    A = align_to(R, Rref, sym)
    U, _, Vt = np.linalg.svd(A.mean(0))
    D = np.diag([1.0, 1.0, np.linalg.det(U @ Vt)])
    return U @ D @ Vt


# --------------------------------------------------------------------------------------------
# grouping
# --------------------------------------------------------------------------------------------
def group_grains(
    R: np.ndarray,
    neighbours: Sequence[Sequence[int]],
    sym: np.ndarray,
    thresh_deg: float = 2.0,
    mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Grain label (n,) of every voxel: flood fill over neighbour pairs whose misorientation is
    below ``thresh_deg``. Labels run 0, 1, ... by decreasing grain size (ties: lowest first voxel).
    Voxels outside ``mask`` get -1."""
    n = len(R)
    R = np.asarray(R, dtype=np.float64)
    ok = np.ones(n, bool) if mask is None else np.asarray(mask, bool)
    lab = -np.ones(n, dtype=np.int64)
    groups: List[List[int]] = []
    for s in range(n):
        if lab[s] >= 0 or not ok[s]:
            continue
        g = len(groups)
        lab[s] = g
        members = [s]
        q = deque([s])
        while q:
            v = q.popleft()
            nb = [u for u in neighbours[v] if lab[u] < 0 and ok[u]]
            if not nb:
                continue
            mis = misorientation_deg(R[nb], np.broadcast_to(R[v], (len(nb), 3, 3)), sym)
            for u, m in zip(nb, np.atleast_1d(mis)):
                if m < thresh_deg and lab[u] < 0:
                    lab[u] = g
                    members.append(u)
                    q.append(u)
        groups.append(members)
    order = sorted(range(len(groups)), key=lambda g: (-len(groups[g]), min(groups[g])))
    remap = np.empty(len(groups), dtype=np.int64)
    remap[order] = np.arange(len(groups))
    out = lab.copy()
    out[lab >= 0] = remap[lab[lab >= 0]]
    return out


def shell_of(members: Sequence[int], neighbours: Sequence[Sequence[int]]) -> List[int]:
    """Voxels next to ``members`` that are not members (the one-voxel shell), sorted."""
    ms = set(int(v) for v in members)
    sh = {int(u) for v in ms for u in neighbours[v] if int(u) not in ms}
    return sorted(sh)


def boundary_members(
    members: Sequence[int], neighbours: Sequence[Sequence[int]]
) -> List[int]:
    """Members with at least one neighbour outside the grain."""
    ms = set(int(v) for v in members)
    return [int(v) for v in members if any(int(u) not in ms for u in neighbours[v])]


# --------------------------------------------------------------------------------------------
# spot model
# --------------------------------------------------------------------------------------------
@dataclass
class SpotRows:
    """Peak geometry of a set of voxels: one row per (voxel, reflection, omega branch) passing
    Bragg and eta. ``col``/``row`` (n, D, 3) are the continuous pixel coordinates of the projected
    triangle vertices on each detector; ``fpos`` the continuous frame position; ``frame`` the frame
    of the truncated bin (-1: outside the scan)."""

    vox: Tensor
    refl: Tensor
    branch: Tensor
    fpos: Tensor
    frame: Tensor
    col: Tensor
    row: Tensor
    hit: Tensor
    sin_eta: Optional[Tensor] = None  # |sin eta| of the diffracted ray (rocking window width)
    spread: Optional["SpotModel"] = None  # the model, only when its PSF or rocking is on


def strain_to_m(strain: Tensor) -> Tensor:
    """M = (I + eps)^-1 (n, 3, 3) from crystal-frame symmetric strain (n, 6)
    [e11 e22 e33 e23 e13 e12]; differentiable."""
    e = strain
    z = torch.stack(
        [e[:, 0], e[:, 5], e[:, 4], e[:, 5], e[:, 1], e[:, 3], e[:, 4], e[:, 3], e[:, 2]], dim=1
    ).reshape(-1, 3, 3)
    return torch.linalg.inv(torch.eye(3, dtype=z.dtype) + z)


def m_to_strain(M: np.ndarray) -> np.ndarray:
    """Inverse of ``strain_to_m``: eps = M^-1 - I as (n, 6) [e11 e22 e33 e23 e13 e12]."""
    E = np.linalg.inv(np.asarray(M, np.float64).reshape(-1, 3, 3)) - np.eye(3)
    E = 0.5 * (E + np.swapaxes(E, 1, 2))
    return np.stack([E[:, 0, 0], E[:, 1, 1], E[:, 2, 2], E[:, 1, 2], E[:, 0, 2], E[:, 0, 1]], 1)


class SpotModel:
    """``SmoothCost.predict`` vectorised over voxels (own orientation and triangle per voxel).

    Built from a ``SmoothCost`` (which holds the reflection list and detector constants).

    Two opt-in blurs of the predicted coverage points (``spread_points``), both default off
    (0 = off, the model is then bit for bit what it was):

    * ``psf_sigma_px``: detector point-spread, a Gaussian of this sigma (pixels) in column and
      row. Each coverage point's mass goes to 3 x 3 fixed offsets (the 3-point Gauss-Hermite
      nodes 0, +-sqrt(3) sigma with weights 2/3, 1/6, 1/6 per axis), which reproduce the
      Gaussian's variance and fourth moment; the weights sum to 1, so mass is conserved.
    * ``rocking_deg``: rocking-curve width (degrees of omega at |sin eta| = 1). A row's mass is
      spread over the omega frames that its window omega* +- (rocking_deg / |sin eta|) / 2
      overlaps (``docs/omega_peak_width_derivation.md``), weighted by the overlap. The spread
      points sit at the frame centres; mass in frames outside the scan is dropped."""

    def __init__(
        self, sc: Any, phase: int = 0, psf_sigma_px: float = 0.0, rocking_deg: float = 0.0
    ) -> None:
        if psf_sigma_px < 0 or rocking_deg < 0:
            raise ValueError("psf_sigma_px and rocking_deg must be >= 0")
        self.psf_sigma_px = float(psf_sigma_px)
        self.rocking_deg = float(rocking_deg)
        self.sc = sc
        self.g_hkl, self.g_mag = sc._g[phase]
        self.K = int(self.g_hkl.shape[0])
        self._consts: Dict[int, Any] = {}
        self.n_det = sc.n_det

    def _bragg(self, n: int) -> Any:
        c = self._consts.get(n)
        if c is None:
            c = self._consts[n] = BraggConstants(
                self.g_mag.repeat(n), self.sc.energy, self.sc.chi
            )
        return c

    def predict(self, R: Tensor, V: Tensor, strain: Optional[Tensor] = None) -> SpotRows:
        """Rows of orientations ``R`` (N, 3, 3) on triangles ``V`` (N, 3, 3) (sample frame).

        ``strain`` (optional, (N, 6) [e11 e22 e33 e23 e13 e12], crystal-frame symmetric strain):
        g = R M g0 with M = (I + eps)^-1, and the Bragg angle and omega come from the strained
        |g| (a pure dilation eps = d I gives |g| = |g0| / (1 + d)). Differentiable in ``strain``.
        ``None`` (default) is the unstrained model, bit for bit."""
        sc = self.sc
        R = R.to(DT)
        V = V.to(DT)
        N, K = int(R.shape[0]), self.K
        if strain is None:
            g_lab = torch.einsum("bij,kj->bki", R, self.g_hkl).reshape(N * K, 3)
            consts = self._bragg(N)
        else:
            A = R @ strain_to_m(strain.to(DT))
            g_lab = torch.einsum("bij,kj->bki", A, self.g_hkl).reshape(N * K, 3)
            consts = BraggConstants(torch.linalg.norm(g_lab, dim=1), sc.energy, sc.chi)
        br = scattering_omegas_fixed_beam(g_lab, consts)
        om = torch.stack([br.omega1, br.omega2], dim=1)
        obs = br.observable[:, None].expand(N * K, 2)
        nrm = g_lab / torch.linalg.norm(g_lab, dim=1, keepdim=True)
        nrm = (sc.base @ nrm.T).T
        c, s = torch.cos(om), torch.sin(om)
        nx, ny, nz = nrm[:, 0:1], nrm[:, 1:2], nrm[:, 2:3]
        ln = torch.stack([c * nx - s * ny, s * nx + c * ny, nz.expand_as(c)], dim=-1)
        dot = (ln * sc.beam).sum(-1, keepdim=True)
        ray = sc.beam - 2.0 * dot * ln
        rn = torch.linalg.norm(ray, dim=-1)
        eta = torch.atan2(ray[..., 1].abs() / rn, ray[..., 2].abs() / rn)
        valid = obs & (eta < sc.eta_limit) & (rn > 0)
        if sc.min_sin_eta > 0:
            valid = valid & (torch.sin(eta) >= sc.min_sin_eta)
        sel = torch.nonzero(valid.reshape(-1))[:, 0]
        bk = sel // 2
        branch = sel % 2
        b = bk // K
        refl = bk % K
        omega = om.reshape(-1)[sel]
        ray_s = ray.reshape(-1, 3)[sel] / rn.reshape(-1)[sel, None]
        fpos = (omega - sc.low) / sc.width
        fbin = torch.trunc(fpos).to(torch.long).clamp(-1, sc.n_bins)
        frame = sc.index_t[fbin]
        P = torch.einsum("ij,nsj->nsi", sc.base, V)[b]  # (n, 3, 3) rotated vertices
        cw, sw = torch.cos(omega)[:, None], torch.sin(omega)[:, None]
        px = cw * P[..., 0] - sw * P[..., 1] + sc.trans[0]
        py = sw * P[..., 0] + cw * P[..., 1] + sc.trans[1]
        pz = P[..., 2] + sc.trans[2]
        pts = torch.stack([px, py, pz], dim=-1)
        cols, rows, hits = [], [], []
        for g in sc._det:
            den = ray_s @ g["n"]
            num = -((pts @ g["n"]) + g["d"])
            par = den.abs() < 1e-12
            t = num / torch.where(par, torch.ones_like(den), den)[:, None]
            hit = (~par)[:, None] & (t > 0)
            q = pts + t[..., None] * ray_s[:, None, :] - g["off"]
            cols.append(((q @ g["bj"]) + g["hw"]) / g["pw"])
            rows.append(((q @ g["bk"]) + g["hh"]) / g["ph"])
            hits.append(hit)
        return SpotRows(
            vox=b,
            refl=refl,
            branch=branch,
            fpos=fpos,
            frame=frame,
            col=torch.stack(cols, 1),
            row=torch.stack(rows, 1),
            hit=torch.stack(hits, 1),
            sin_eta=torch.sin(eta).reshape(-1)[sel],
            spread=self if self.spreading else None,
        )

    @property
    def spreading(self) -> bool:
        """True when a PSF or a rocking curve is on."""
        return self.psf_sigma_px > 0 or self.rocking_deg > 0

    def rocking_weights(self, fpos: Tensor, sin_eta: Tensor) -> Tuple[Tensor, Tensor]:
        """(bins (n, K) long, weights (n, K)): the frame bins that the window fpos +- h overlaps,
        h = rocking_deg / |sin eta| / 2 in frames, and the overlap / (2 h) of each (rows of the
        weights sum to 1 before bins outside the scan are dropped). |sin eta| is floored at 0.01."""
        sc = self.sc
        h = (math.radians(self.rocking_deg) / sin_eta.abs().clamp_min(0.01) / 2.0) / sc.width
        h = h.clamp_min(1e-9)
        lo, hi = fpos - h, fpos + h
        b0 = torch.floor(lo).detach().to(torch.long)
        k = int((torch.floor(hi).detach().to(torch.long) - b0).max()) + 1 if len(b0) else 1
        bins = b0[:, None] + torch.arange(k)
        bf = bins.to(fpos.dtype)
        ov = (torch.minimum(hi[:, None], bf + 1.0) - torch.maximum(lo[:, None], bf)).clamp_min(0.0)
        w = ov / (2.0 * h)[:, None]
        live = sc.index_t[bins.clamp(-1, sc.n_bins)] >= 0
        return bins, w * live.to(w.dtype)

    def spread_points(
        self, rows: SpotRows, pc: Tensor, pr: Tensor, m: Tensor
    ) -> Tuple[Tensor, Tensor, Tensor, Tensor]:
        """(cols, rows, fpos, mass), each (n, P'), of the coverage points ``pc``, ``pr``, ``m``
        (n, P) of ``rows`` after the PSF and rocking spread (see the class doc). With both off the
        points are returned as they are, ``fpos`` the row's continuous frame position."""
        n, P = pc.shape
        fp = rows.fpos[:, None].expand_as(pc)
        if self.psf_sigma_px > 0:
            r3 = math.sqrt(3.0) * self.psf_sigma_px
            nodes = torch.tensor([-r3, 0.0, r3], dtype=pc.dtype)
            wts = torch.tensor([1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0], dtype=pc.dtype)
            dx = nodes.repeat_interleave(3)
            dy = nodes.repeat(3)
            ww = (wts[:, None] * wts[None, :]).reshape(-1)
            pc = (pc[:, :, None] + dx).reshape(n, -1)
            pr = (pr[:, :, None] + dy).reshape(n, -1)
            m = (m[:, :, None] * ww).reshape(n, -1)
            fp = fp[:, :, None].expand(n, P, 9).reshape(n, -1)
        if self.rocking_deg > 0:
            assert rows.sin_eta is not None
            bins, w = self.rocking_weights(rows.fpos, rows.sin_eta)
            k = bins.shape[1]
            Q = pc.shape[1]
            pc = pc[:, :, None].expand(n, Q, k).reshape(n, -1)
            pr = pr[:, :, None].expand(n, Q, k).reshape(n, -1)
            m = (m[:, :, None] * w[:, None, :]).reshape(n, -1)
            fp = (bins.to(pc.dtype) + 0.5)[:, None, :].expand(n, Q, k).reshape(n, -1)
        return pc, pr, fp, m


def subtriangle_bary(n_sub: int) -> np.ndarray:
    """Barycentric (u, v) (n_sub**2, 2) of the centroids of the ``n_sub**2`` equal sub-triangles of
    a triangle (point = V0 + u (V1 - V0) + v (V2 - V0))."""
    pts = []
    for i in range(n_sub):
        for j in range(n_sub - i):
            pts.append(((i + 1 / 3) / n_sub, (j + 1 / 3) / n_sub))
            if i + j <= n_sub - 2:
                pts.append(((i + 2 / 3) / n_sub, (j + 2 / 3) / n_sub))
    out = np.asarray(pts, dtype=np.float64)
    assert len(out) == n_sub * n_sub
    return out


def coverage_points(col: Tensor, row: Tensor, bary: Tensor) -> Tuple[Tensor, Tensor, Tensor]:
    """(cols (n, P), rows (n, P), mass (n, P)) of sub-triangle centroids of projected triangles
    with vertices ``col``/``row`` (n, 3); the mass is the projected area (pixels^2) / P."""
    u, v = bary[:, 0], bary[:, 1]
    c0, c1, c2 = col[:, 0:1], col[:, 1:2], col[:, 2:3]
    r0, r1, r2 = row[:, 0:1], row[:, 1:2], row[:, 2:3]
    pc = c0 + u * (c1 - c0) + v * (c2 - c0)
    pr = r0 + u * (r1 - r0) + v * (r2 - r0)
    area = 0.5 * ((c1 - c0) * (r2 - r0) - (c2 - c0) * (r1 - r0)).abs()  # (n, 1)
    mass = (area / bary.shape[0]).expand_as(pc)
    return pc, pr, mass


# --------------------------------------------------------------------------------------------
# region
# --------------------------------------------------------------------------------------------
@dataclass
class Spot:
    """One reflection of a grain on one detector: key, window and observed pixels."""

    key: Tuple[int, int, int]  # (reflection, branch, detector)
    c0: int
    c1: int  # window columns [c0, c1)
    r0: int
    r1: int
    f0: int
    f1: int  # window frames [f0, f1)
    obs: np.ndarray  # (m, 3) pixel centres (col + 0.5, row + 0.5, frame + 0.5)
    obs_val: np.ndarray  # (m,) clean intensity (1 when no store)
    ratio: float = float("nan")  # observed lit pixels / predicted area at the start


@dataclass
class GrainRegion:
    """The spots of one grain and the fixed shell prediction in their windows."""

    members: List[int]
    shell: List[int]
    spots: List[Spot]
    key_index: Dict[Tuple[int, int, int], int]
    shell_pts: List[np.ndarray] = field(default_factory=list)  # per spot (k, 3)
    shell_mass: List[np.ndarray] = field(default_factory=list)  # per spot (k,)


class ImageSource:
    """Binary stack (memmap, (n_det * n_omega, rows, cols)) and optional sparse intensity store."""

    def __init__(self, stack: np.ndarray, n_omega: int, store: Optional[Dict[str, np.ndarray]] = None):
        self.stack = stack
        self.n_omega = int(n_omega)
        self.store = store

    def lit(self, det: int, f0: int, f1: int, r0: int, r1: int, c0: int, c1: int) -> Tuple[np.ndarray, np.ndarray]:
        """(points (m, 3) pixel centres with frame + 0.5, values (m,)) of lit pixels in a box."""
        pts, vals = [], []
        for f in range(max(f0, 0), min(f1, self.n_omega)):
            idx = det * self.n_omega + f
            blk = np.asarray(self.stack[idx, r0:r1, c0:c1])
            rr, cc = np.nonzero(blk)
            if not len(rr):
                continue
            p = np.stack([cc + c0 + 0.5, rr + r0 + 0.5, np.full(len(rr), f + 0.5)], 1)
            pts.append(p)
            if self.store is not None:
                vals.append(self._values(idx, rr + r0, cc + c0))
            else:
                vals.append(np.ones(len(rr)))
        if not pts:
            return np.zeros((0, 3)), np.zeros(0)
        return np.concatenate(pts), np.concatenate(vals)

    def _values(self, idx: int, rr: np.ndarray, cc: np.ndarray) -> np.ndarray:
        st = self.store
        assert st is not None
        a, b = int(st["frame"][idx]), int(st["frame"][idx + 1])
        R = st["row"][a:b].astype(np.int64)
        C = st["col"][a:b].astype(np.int64)
        key = R * 4096 + C  # sorted by (row, col)
        q = rr.astype(np.int64) * 4096 + cc.astype(np.int64)
        pos = np.searchsorted(key, q)
        pos = np.clip(pos, 0, max(len(key) - 1, 0))
        ok = (len(key) > 0) & (key[pos] == q) if len(key) else np.zeros(len(q), bool)
        out = np.where(ok, st["val"][a:b][pos] if len(key) else 0.0, 0.0)
        return out.astype(np.float64)


def _rows_points(rows: SpotRows, bary: Tensor, det: int) -> Tuple[Tensor, Tensor, Tensor, Tensor]:
    """(cols, rows, mass, ok) per row on detector ``det``; ok: all vertices hit the plane and the
    row has a frame."""
    pc, pr, m = coverage_points(rows.col[:, det], rows.row[:, det], bary)
    ok = rows.hit[:, det].all(1) & (rows.frame >= 0)
    return pc, pr, m, ok


def _rows_cloud(
    rows: SpotRows, bary: Tensor, det: int
) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor]:
    """(cols, rows, fpos, mass, ok) of the coverage points per row on detector ``det``: as
    ``_rows_points`` plus the continuous frame position of each point; when the model that made
    ``rows`` has its PSF or rocking on, the points are the spread ones
    (``SpotModel.spread_points``)."""
    pc, pr, m, ok = _rows_points(rows, bary, det)
    if rows.spread is not None:
        pc, pr, fp, m = rows.spread.spread_points(rows, pc, pr, m)
    else:
        fp = rows.fpos[:, None].expand_as(pc)
    return pc, pr, fp, m, ok


def _touched_pixels(rows: SpotRows, ii: Tensor, det: int, bary: Tensor) -> int:
    """Distinct (frame, row, col) pixels touched by the coverage points of rows ``ii``."""
    pc, pr, _m = coverage_points(rows.col[ii, det], rows.row[ii, det], bary)
    f = rows.frame[ii][:, None].expand_as(pc)
    key = torch.stack([f.reshape(-1), torch.floor(pr).long().reshape(-1), torch.floor(pc).long().reshape(-1)], 1)
    return int(torch.unique(key, dim=0).shape[0])


def _touching_components(
    obs: np.ndarray,
    box: Tuple[int, int, int, int, int, int],
    pc: Tensor,
    pr: Tensor,
    frames: Tensor,
    keep_px: int,
    frame_dilation: int = 1,
) -> np.ndarray:
    """Mask (m,) of observed pixels whose component touches the dilated predicted footprint
    (``keep_px`` pixels, ``frame_dilation`` frames)."""
    from scipy import ndimage

    f0, f1, r0, r1, c0, c1 = box
    shape = (f1 - f0, r1 - r0, c1 - c0)
    grid = np.zeros(shape, dtype=bool)
    fi = (obs[:, 2] - 0.5).astype(np.int64) - f0
    ri = (obs[:, 1] - 0.5).astype(np.int64) - r0
    ci = (obs[:, 0] - 0.5).astype(np.int64) - c0
    grid[fi, ri, ci] = True
    st = np.zeros((3, 3, 3), bool)
    st[1] = True  # 8-connected in a frame; frames are not joined (another grain's spot in the
    # next frame at the same pixels would merge)
    lab, _n = ndimage.label(grid, structure=st)
    foot = np.zeros(shape, dtype=bool)
    pf = (frames[:, None].expand_as(pc).reshape(-1).numpy() - f0).astype(np.int64)
    pi = np.floor(pr.reshape(-1).numpy()).astype(np.int64) - r0
    pj = np.floor(pc.reshape(-1).numpy()).astype(np.int64) - c0
    ok = (pf >= 0) & (pf < shape[0]) & (pi >= 0) & (pi < shape[1]) & (pj >= 0) & (pj < shape[2])
    foot[pf[ok], pi[ok], pj[ok]] = True
    k = 2 * keep_px + 1
    foot = ndimage.maximum_filter(foot, size=(2 * frame_dilation + 1, k, k))
    hit = np.unique(lab[foot & (lab > 0)])
    return np.isin(lab[fi, ri, ci], hit)


def build_region(
    model: SpotModel,
    images: ImageSource,
    R_start: np.ndarray,
    V: np.ndarray,
    members: Sequence[int],
    shell: Sequence[int],
    bary: Tensor,
    pad_px: int = 6,
    pad_frames: int = 1,
    min_obs: int = 1,
    keep_px: Optional[int] = 3,
    ratio: Optional[Tuple[float, float]] = (0.8, 1.6),
    strain: Optional[np.ndarray] = None,
) -> GrainRegion:
    """Spots of the grain ``members`` at their start orientations ``R_start`` (n_all, 3, 3) over
    the layer's triangles ``V`` (n_all, 3, 3), and the shell prediction in their windows.

    ``keep_px``: keep only the observed connected components (8-connected in a frame) that come
    within ``keep_px`` pixels and one frame (1 + the rocking half-window in frames when the model
    spreads) of a start-predicted coverage point; other grains'
    spots in the window are dropped. None keeps every lit pixel.

    ``ratio``: (lo, hi) bounds on observed lit pixels / pixels touched by the start prediction
    (counted with 8 x 8 sub-triangles per voxel, as the stock render lights every pixel a spot
    touches); a
    spot outside them is dropped (overlap with another grain's spot, or a spot cut by a detector
    edge or a gap). The ratio of every spot is kept in ``Spot.ratio``. None keeps all. The
    denominator is the footprint of the unspread start prediction even when the model spreads
    (PSF or rocking), while the observed count covers the widened box; for data with PSF or
    rocking the ratio can exceed ``ratio[1]`` and drop spots, so widen ``ratio`` (or pass None)
    for such data.

    ``strain`` (optional, (n_all, 6), see ``SpotModel.predict``): the start prediction of members
    and shell is made with it (default None: unstrained, unchanged).

    With the model's PSF or rocking on, the pixel pad grows by ceil(sqrt(3) psf_sigma_px) and the
    frame pad by ceil(rocking half-window in frames) of the spot's narrowest |sin eta|, so the
    spread mass falls inside the observed box. Scan ends and gaps: the model keeps rows whose
    omega* is outside the scan or in a gap and spreads them into live frames, whereas the stock
    renderer and the strict cost drop such peaks first (``omega_valid``, ``valid_indices``)."""
    st = None if strain is None else torch.as_tensor(np.asarray(strain), dtype=DT)
    members = [int(v) for v in members]
    shell = [int(v) for v in shell]
    fine = torch.as_tensor(subtriangle_bary(8))
    with torch.no_grad():
        rows = model.predict(
            torch.as_tensor(R_start[members]), torch.as_tensor(V[members]),
            None if st is None else st[members],
        )
    spots: List[Spot] = []
    H, W = images.stack.shape[1], images.stack.shape[2]
    for det in range(model.n_det):
        pc, pr, _m, ok = _rows_points(rows, bary, det)
        cols_v = rows.col[:, det]
        rows_v = rows.row[:, det]
        inb = ok & (cols_v.min(1).values >= 0) & (cols_v.max(1).values < W)
        inb = inb & (rows_v.min(1).values >= 0) & (rows_v.max(1).values < H)
        keys: Dict[Tuple[int, int], List[int]] = {}
        for i in torch.nonzero(inb)[:, 0].tolist():
            keys.setdefault((int(rows.refl[i]), int(rows.branch[i])), []).append(i)
        for (refl, br), idx in sorted(keys.items()):
            ii = torch.as_tensor(idx)
            c = cols_v[ii]
            r = rows_v[ii]
            fr = rows.frame[ii]
            spread_f = 0
            pad_c = pad_px
            if model.spreading:
                # the spread mass reaches the rocking half-window and the PSF's outermost node
                # beyond the unspread rows: widen the box (and the frames) by that much
                pad_c += int(math.ceil(math.sqrt(3.0) * model.psf_sigma_px))
                if model.rocking_deg > 0 and rows.sin_eta is not None:
                    hmax = (
                        math.radians(model.rocking_deg)
                        / float(rows.sin_eta[ii].abs().clamp_min(0.01).min())
                        / 2.0
                        / model.sc.width
                    )
                    spread_f = int(math.ceil(hmax))
            c0 = max(int(math.floor(float(c.min()))) - pad_c, 0)
            c1 = min(int(math.floor(float(c.max()))) + pad_c + 1, W)
            r0 = max(int(math.floor(float(r.min()))) - pad_c, 0)
            r1 = min(int(math.floor(float(r.max()))) + pad_c + 1, H)
            f0 = int(fr.min()) - pad_frames - spread_f
            f1 = int(fr.max()) + pad_frames + spread_f + 1
            obs, val = images.lit(det, f0, f1, r0, r1, c0, c1)
            if keep_px is not None and len(obs):
                keep = _touching_components(
                    obs, (f0, f1, r0, r1, c0, c1), pc[ii], pr[ii], fr, keep_px,
                    1 + spread_f,
                )
                obs, val = obs[keep], val[keep]
            if len(obs) < min_obs:
                continue
            q = len(obs) / max(_touched_pixels(rows, ii, det, fine), 1)
            if ratio is not None and not (ratio[0] <= q <= ratio[1]):
                continue
            spots.append(Spot((refl, br, det), c0, c1, r0, r1, f0, f1, obs, val, q))
    key_index = {s.key: i for i, s in enumerate(spots)}
    region = GrainRegion(members, shell, spots, key_index)
    region.shell_pts = [np.zeros((0, 3)) for _ in spots]
    region.shell_mass = [np.zeros(0) for _ in spots]
    if shell and spots:
        with torch.no_grad():
            srows = model.predict(
                torch.as_tensor(R_start[shell]), torch.as_tensor(V[shell]),
                None if st is None else st[shell],
            )
            for det in range(model.n_det):
                pc, pr, fpall, m, ok = _rows_cloud(srows, bary, det)
                for j, s in enumerate(spots):
                    if s.key[2] != det:
                        continue
                    inside = (pc >= s.c0) & (pc < s.c1) & (pr >= s.r0) & (pr < s.r1)
                    rsel = ok & (srows.frame >= s.f0) & (srows.frame < s.f1) & inside.any(1)
                    if not bool(rsel.any()):
                        continue
                    ins = inside[rsel]
                    p = torch.stack([pc[rsel], pr[rsel], fpall[rsel]], -1)[ins].numpy()
                    region.shell_pts[j] = np.concatenate([region.shell_pts[j], p])
                    region.shell_mass[j] = np.concatenate(
                        [region.shell_mass[j], m[rsel][ins].numpy()]
                    )
    return region


# --------------------------------------------------------------------------------------------
# packed prediction of a grain over its region
# --------------------------------------------------------------------------------------------
@dataclass
class Packed:
    """Padded per-spot tensors: predicted points/mass and observed points/mass."""

    x: Tensor  # (S, P, 3) predicted (col, row, fpos)
    a: Tensor  # (S, P) predicted mass (unnormalised)
    y: Tensor  # (S, O, 3) observed pixel centres
    b: Tensor  # (S, O) observed mass (normalised per spot)
    yv: Tensor  # (S, O) observed intensity (raw)


def bin_observed(obs: np.ndarray, w: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """Observed pixels merged into k x k pixel bins per frame: mass-weighted centre, summed mass."""
    if k <= 1 or not len(obs):
        return obs, w
    key = np.stack([np.floor(obs[:, 0] / k), np.floor(obs[:, 1] / k), obs[:, 2]], 1)
    uk, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    ws = np.bincount(inv, weights=w, minlength=len(uk))
    cen = np.stack([np.bincount(inv, weights=w * obs[:, j], minlength=len(uk)) for j in range(3)], 1)
    return cen / ws[:, None], ws


def pack_observed(
    region: GrainRegion, intensity: bool = False, bin_px: int = 1
) -> Tuple[Tensor, Tensor, Tensor]:
    """(y (S, O, 3), b (S, O) normalised mass, raw mass) of the observed pixels, merged into
    ``bin_px`` square bins when ``bin_px`` > 1 (fewer points for large spots)."""
    S = len(region.spots)
    packs = []
    for s in region.spots:
        w = np.asarray(s.obs_val, float) if intensity else np.ones(len(s.obs))
        packs.append(bin_observed(np.asarray(s.obs, float), w, bin_px))
    O = max(len(p[0]) for p in packs)
    y = torch.zeros((S, O, 3), dtype=DT)
    b = torch.zeros((S, O), dtype=DT)
    yv = torch.zeros((S, O), dtype=DT)
    for i, (o, w) in enumerate(packs):
        m = len(o)
        y[i, :m] = torch.as_tensor(o)
        b[i, :m] = torch.as_tensor(w / w.sum())
        yv[i, :m] = torch.as_tensor(w)
    return y, b, yv


def pack_predicted(
    region: GrainRegion,
    rows: SpotRows,
    bary: Tensor,
    n_det: int,
    with_shell: bool = True,
) -> Tuple[Tensor, Tensor]:
    """(x (S, P, 3), a (S, P)) predicted coverage points of the grain's rows per spot, plus the
    fixed shell points. Rows whose key is not a spot of the region are dropped (a reflection that
    appears only away from the start is not scored; it has no window). With the PSF or rocking of
    the model on (``rows.spread``) the points are the spread ones."""
    S = len(region.spots)
    groups: List[List[Tuple[Tensor, Tensor]]] = [[] for _ in range(S)]
    for det in range(n_det):
        pc, pr, fpall, m, ok = _rows_cloud(rows, bary, det)
        keyv = (rows.refl * 2 + rows.branch).tolist()
        okl = ok.tolist()
        for i, kv in enumerate(keyv):
            if not okl[i]:
                continue
            j = region.key_index.get((kv // 2, kv % 2, det))
            if j is None:
                continue
            p = torch.stack([pc[i], pr[i], fpall[i]], -1)
            mi = m[i]
            if rows.spread is not None:
                # spread points past the scan or in a gap carry no mass: drop them so that one
                # wide-window row does not pad every spot (defaults keep every point)
                live = mi > 0
                p, mi = p[live], mi[live]
            groups[j].append((p, mi))
    sizes = []
    for j in range(S):
        n = sum(int(g[1].shape[0]) for g in groups[j])
        if with_shell:
            n += len(region.shell_mass[j])
        sizes.append(n)
    P = max(max(sizes), 1)
    x = torch.zeros((S, P, 3), dtype=DT)
    a = torch.zeros((S, P), dtype=DT)
    xs, as_ = [], []
    for j in range(S):
        parts = [g[0] for g in groups[j]]
        ms = [g[1] for g in groups[j]]
        if with_shell and len(region.shell_mass[j]):
            parts.append(torch.as_tensor(region.shell_pts[j]))
            ms.append(torch.as_tensor(region.shell_mass[j]))
        if parts:
            pj = torch.cat(parts, 0)
            mj = torch.cat(ms, 0)
            pj = torch.cat([pj, torch.zeros((P - len(mj), 3), dtype=DT)], 0)
            mj = torch.cat([mj, torch.zeros(P - len(mj), dtype=DT)], 0)
        else:
            pj = torch.zeros((P, 3), dtype=DT)
            mj = torch.zeros(P, dtype=DT)
        xs.append(pj)
        as_.append(mj)
    x = torch.stack(xs)
    a = torch.stack(as_)
    return x, a


# --------------------------------------------------------------------------------------------
# Sinkhorn
# --------------------------------------------------------------------------------------------
def frame_cost(x: Tensor, y: Tensor, lam: float) -> Tensor:
    """Squared ground cost (S, P, O) between predicted (col, row, fpos) and observed pixel centres
    (col, row, frame + 0.5): squared pixel distance plus ``(lam * d_f)**2``, with ``d_f`` the
    distance of ``fpos`` outside the observed frame [frame, frame + 1) (0 inside it)."""
    d2 = ((x[:, :, None, :2] - y[:, None, :, :2]) ** 2).sum(-1)
    df = torch.relu((x[:, :, None, 2] - y[:, None, :, 2]).abs() - 0.5)
    return d2 + (lam * df) ** 2


def _softmin(eps: float, C: Tensor, h: Tensor) -> Tensor:
    """-eps log sum_j exp(h_j - C_ij / eps), over the last axis; h (S, 1, O) log-weights + pot."""
    return -eps * torch.logsumexp(h - C / eps, dim=-1)


def sinkhorn_ot(
    a: Tensor,
    b: Tensor,
    Cab: Tensor,
    eps_list: Sequence[float],
    n_last: int = 1,
    reach2: Optional[float] = None,
) -> Tensor:
    """Entropic OT value per spot (S,) between weights ``a`` (S, P) and ``b`` (S, O) (each summing
    to 1 over positive entries) with cost ``Cab`` (S, P, O). Symmetric Sinkhorn with eps scaling
    (geomloss style): iterations run without gradient, then ``n_last`` updates with gradient so the
    value is differentiable through the cost (envelope theorem). ``reach2`` (cost units) makes it
    unbalanced (KL marginal penalty); None is balanced."""
    tiny = 1e-300
    la = torch.where(a > 0, torch.log(a.clamp_min(tiny)), torch.full_like(a, -1e30))
    lb = torch.where(b > 0, torch.log(b.clamp_min(tiny)), torch.full_like(b, -1e30))

    def damp(eps: float) -> float:
        return 1.0 if reach2 is None else 1.0 / (1.0 + eps / reach2)

    Cba = Cab.transpose(1, 2)
    with torch.no_grad():
        Cd = Cab.detach()
        Cdt = Cd.transpose(1, 2)
        eps0 = eps_list[0]
        f = damp(eps0) * _softmin(eps0, Cd, lb[:, None, :])
        g = damp(eps0) * _softmin(eps0, Cdt, la[:, None, :])
        for eps in eps_list:
            ft = damp(eps) * _softmin(eps, Cd, (lb + g / eps)[:, None, :])
            gt = damp(eps) * _softmin(eps, Cdt, (la + f / eps)[:, None, :])
            f, g = 0.5 * (f + ft), 0.5 * (g + gt)
    eps = eps_list[-1]
    for _ in range(n_last):
        f, g = (
            damp(eps) * _softmin(eps, Cab, (lb + g.detach() / eps)[:, None, :]),
            damp(eps) * _softmin(eps, Cba, (la + f.detach() / eps)[:, None, :]),
        )
    if reach2 is None:
        return (a * f).sum(-1) + (b * g).sum(-1)
    # unbalanced dual: reach2 * <a, 1 - exp(-f / reach2)> + same for b
    return reach2 * ((a * (1 - torch.exp(-f / reach2))).sum(-1) + (b * (1 - torch.exp(-g / reach2))).sum(-1))


def eps_schedule(blur_px: float, diameter_px: float, scaling: float = 0.7) -> List[float]:
    """Geometric eps list (cost units, px^2) from diameter^2 down to blur^2."""
    out = []
    e = diameter_px ** 2
    while e > blur_px ** 2:
        out.append(e)
        e *= scaling ** 2
    out.append(blur_px ** 2)
    return out


def sinkhorn_divergence(
    x: Tensor,
    a: Tensor,
    y: Tensor,
    b: Tensor,
    blur_px: float,
    lam: float = 10.0,
    diameter_px: float = 40.0,
    scaling: float = 0.7,
    reach_px: Optional[float] = None,
    debias: bool = True,
) -> Tensor:
    """Sinkhorn divergence per spot (S,) (px^2) between predicted points ``x`` (S, P, 3) with mass
    ``a`` (normalised here) and observed ``y`` (S, O, 3) with mass ``b`` (normalised). Debiased:
    OT(a, b) - OT(a, a) / 2 - OT(b, b) / 2."""
    empty = a.sum(-1) <= 0
    if bool(empty.any()):  # a reflection with no predicted mass: one dummy point, fixed penalty
        a = a.clone()
        a[empty, 0] = 1.0
    an = a / a.sum(-1, keepdim=True).clamp_min(1e-300)
    eps_list = eps_schedule(blur_px, diameter_px, scaling)
    r2 = None if reach_px is None else reach_px ** 2
    ot_ab = sinkhorn_ot(an, b, frame_cost(x, y, lam), eps_list, reach2=r2)
    penalty = (diameter_px ** 2) if reach_px is None else reach_px ** 2
    if not debias:
        return torch.where(empty, torch.full_like(ot_ab, penalty), ot_ab)
    ot_aa = sinkhorn_ot(an, an, frame_cost(x, x, lam), eps_list, reach2=r2)
    with torch.no_grad():
        ot_bb = sinkhorn_ot(b, b, frame_cost(y, y, lam), eps_list, reach2=r2)
    out = ot_ab - 0.5 * ot_aa - 0.5 * ot_bb
    return torch.where(empty, torch.full_like(out, penalty), out)


# --------------------------------------------------------------------------------------------
# least squares baseline (soft raster)
# --------------------------------------------------------------------------------------------
def soft_raster(
    region: GrainRegion, x: Tensor, a: Tensor, sigma_px: float, tau_frames: float
) -> List[Tensor]:
    """Predicted coverage on each spot's window grid (F, H, W): each point's mass spread by a
    Gaussian (sigma, pixels) in col and row and a smoothed box in frames (edge width tau)."""
    out = []
    for j, s in enumerate(region.spots):
        m = a[j]
        keep = m > 0
        if not bool(keep.any()):
            out.append(torch.zeros((s.f1 - s.f0, s.r1 - s.r0, s.c1 - s.c0), dtype=DT))
            continue
        xj, mj = x[j][keep], m[keep]
        cc = torch.arange(s.c0, s.c1, dtype=DT) + 0.5
        rr = torch.arange(s.r0, s.r1, dtype=DT) + 0.5
        ff = torch.arange(s.f0, s.f1, dtype=DT)
        wx = _pixel_weight(xj[:, 0:1], cc[None], sigma_px)  # (k, W)
        wy = _pixel_weight(xj[:, 1:2], rr[None], sigma_px)  # (k, H)
        z = xj[:, 2:3]
        wf = torch.sigmoid((z - ff[None]) / tau_frames) - torch.sigmoid((z - ff[None] - 1) / tau_frames)
        A = (mj[:, None, None] * wf[:, :, None] * wy[:, None, :]).reshape(len(mj), -1)
        img = (A.T @ wx).reshape(len(ff), len(rr), len(cc))
        out.append(img)
    return out


def _pixel_weight(p: Tensor, centres: Tensor, sigma: float) -> Tensor:
    """Integral over a unit pixel around ``centres`` of a Gaussian at ``p`` (sigma pixels)."""
    s = sigma * math.sqrt(2.0)
    return 0.5 * (torch.erf((centres + 0.5 - p) / s) - torch.erf((centres - 0.5 - p) / s))


def observed_grids(region: GrainRegion, intensity: bool = True) -> List[Tensor]:
    out = []
    for s in region.spots:
        g = torch.zeros((s.f1 - s.f0, s.r1 - s.r0, s.c1 - s.c0), dtype=DT)
        if len(s.obs):
            f = (s.obs[:, 2] - 0.5).astype(np.int64) - s.f0
            r = (s.obs[:, 1] - 0.5).astype(np.int64) - s.r0
            c = (s.obs[:, 0] - 0.5).astype(np.int64) - s.c0
            v = s.obs_val if intensity else np.ones(len(s.obs))
            g[f, r, c] = torch.as_tensor(v, dtype=DT)
        out.append(g)
    return out


def ls_objective(pred: List[Tensor], obs: List[Tensor]) -> Tensor:
    """Per spot (S,): least squares with the best scale per spot, relative to |obs|^2:
    1 - (P.O)^2 / (|P|^2 |O|^2) (0 for a perfect match up to scale, 1 for no overlap)."""
    vals = []
    for P, O in zip(pred, obs):
        po = (P * O).sum()
        pp = (P * P).sum().clamp_min(1e-300)
        oo = (O * O).sum().clamp_min(1e-300)
        vals.append(1.0 - po * po / (pp * oo))
    return torch.stack(vals)


# --------------------------------------------------------------------------------------------
# fields
# --------------------------------------------------------------------------------------------
def poly_features(d: Tensor, order: int) -> Tensor:
    """Features (n, k) of normalised offsets ``d`` (n, p): linear (p) [+ quadratic p(p+1)/2]."""
    feats = [d]
    if order >= 2:
        p = d.shape[1]
        quad = [d[:, i] * d[:, j] for i in range(p) for j in range(i, p)]
        feats.append(torch.stack(quad, 1))
    return torch.cat(feats, 1)


@dataclass
class FieldA:
    """Stage A field: ``R_v = exp([w0 + C f(d_v)]) R_ref`` with ``theta = [w0 (3), C (3 x k)]``."""

    R_ref: np.ndarray
    centre: np.ndarray
    radius: float
    feats: Tensor  # (n, k)

    @classmethod
    def build(cls, R_ref: np.ndarray, xy: np.ndarray, order: int = 1) -> "FieldA":
        c = xy.mean(0)
        rad = float(np.linalg.norm(xy - c, axis=1).max()) if len(xy) > 1 else 1.0
        rad = rad if rad > 0 else 1.0
        d = torch.as_tensor((xy - c) / rad, dtype=DT)
        return cls(R_ref, c, rad, poly_features(d, order))

    @property
    def n_params(self) -> int:
        return 3 + 3 * int(self.feats.shape[1])

    def rotvecs(self, theta: Tensor) -> Tensor:
        w0 = theta[:3]
        C = theta[3:].reshape(3, -1)
        return w0[None, :] + self.feats @ C.T

    def orientations(self, theta: Tensor) -> Tensor:
        return rotvec_to_matrix(self.rotvecs(theta)) @ torch.as_tensor(self.R_ref, dtype=DT)

    def fit_theta(self, R: np.ndarray, sym: np.ndarray) -> np.ndarray:
        """Least-squares theta of given voxel orientations (symmetry-aligned to R_ref)."""
        A = align_to(R, self.R_ref, sym)
        w = matrix_to_rotvec(np.einsum("nij,kj->nik", A, self.R_ref))  # A R_ref^T
        X = np.concatenate([np.ones((len(w), 1)), self.feats.numpy()], 1)
        coef, *_ = np.linalg.lstsq(X, w, rcond=None)  # (1 + k, 3)
        return np.concatenate([coef[0], coef[1:].T.reshape(-1)])


# --------------------------------------------------------------------------------------------
# one grain's fit
# --------------------------------------------------------------------------------------------
@dataclass
class ObjectiveSpec:
    """Objective settings. ``kind``: "sinkhorn" (main), "ls" (baseline) or "strict" (mean strict
    binary cost over the grain's voxels, through ``GrainProblem.strict_fn``; not differentiable)."""

    kind: str = "sinkhorn"
    blur_px: float = 1.0
    diameter_px: float = 12.0
    reach_px: Optional[float] = 3.0
    lam: float = 10.0
    sigma_px: float = 0.5  # ls: Gaussian splat
    tau_frames: float = 0.1  # ls: frame edge width
    intensity: bool = False  # sinkhorn: weight observed pixels by clean intensity
    dilate: float = 0.0  # sinkhorn: split each coverage point over 4 offsets (+-dilate px in col
    # and row): a pixel-box blur of the prediction, as lit pixels are those a spot touches


class GrainProblem:
    """The fit of one grain: region at the map's orientations, objective, stage A and stage B."""

    def __init__(
        self,
        model: SpotModel,
        images: ImageSource,
        R_map: np.ndarray,
        V: np.ndarray,
        xy: np.ndarray,
        members: Sequence[int],
        shell: Sequence[int],
        sym: np.ndarray,
        n_sub: int = 2,
        order: int = 1,
        spec: Optional[ObjectiveSpec] = None,
        region_kw: Optional[Dict[str, Any]] = None,
        max_points: int = 240,
        chunk_elems: int = 3_000_000,
        strict_fn: Optional[Any] = None,
        strain: Optional[np.ndarray] = None,
    ) -> None:
        self.model = model
        # optional crystal-frame strain (n_voxels_total, 6), used by the region and ``objective``
        self.strain = (
            None if strain is None else torch.as_tensor(np.asarray(strain)[list(members)], dtype=DT)
        )
        # strict_fn(R (n, 3, 3) numpy, members) -> per-voxel strict costs; used by kind "strict"
        self.strict_fn = strict_fn
        self.n_strict = 0
        self.members = [int(v) for v in members]
        self.shell = [int(v) for v in shell]
        self.sym = sym
        self.spec = spec or ObjectiveSpec()
        # at most about ``max_points`` predicted points per spot: fewer sub-triangles for big grains
        while n_sub > 1 and n_sub * n_sub * len(self.members) > max_points:
            n_sub -= 1
        self.n_sub = n_sub
        self.chunk_elems = int(chunk_elems)
        self.bary = torch.as_tensor(subtriangle_bary(n_sub))
        self.V = torch.as_tensor(V[self.members], dtype=DT)
        # one symmetry representative for the whole grain: reflection indices (the spot keys)
        # depend on it
        R_m = np.asarray(R_map[self.members], dtype=np.float64)
        self.R_start = align_to(R_m, R_m[0], sym)
        R_al = np.array(R_map, dtype=np.float64, copy=True)
        R_al[self.members] = self.R_start
        if strain is not None:  # the strain follows the members' new symmetry representative
            strain = np.array(strain, dtype=np.float64, copy=True)
            strain[self.members] = rotate_strain(strain[self.members], align_ops(R_m, R_m[0], sym))
            self.strain = torch.as_tensor(strain[self.members], dtype=DT)
        self.region = build_region(
            model, images, R_al, V, self.members, self.shell, self.bary,
            **dict(region_kw or {}, strain=strain),
        )
        o_max = max((len(sp.obs) for sp in self.region.spots), default=0)
        self.bin_px = 1
        while o_max / (self.bin_px ** 2) > max_points * 1.5 and self.bin_px < 4:
            self.bin_px += 1
        self.y, self.b, _ = pack_observed(self.region, self.spec.intensity, self.bin_px)
        self._obs_grids: Optional[List[Tensor]] = None
        R_ref = mean_orientation(R_m, sym, R_m[0])
        self.field = FieldA.build(R_ref, xy[self.members], order)
        self.theta0 = self.field.fit_theta(self.R_start, sym)
        self.n_evals = 0

    @property
    def n_spots(self) -> int:
        return len(self.region.spots)

    def objective(self, R: Tensor, spec: Optional[ObjectiveSpec] = None) -> Tensor:
        """Mean over spots of the objective at member orientations ``R`` (n, 3, 3)."""
        sp = spec or self.spec
        self.n_evals += 1
        if sp.kind == "strict":
            assert self.strict_fn is not None, "kind 'strict' needs strict_fn"
            self.n_strict += len(self.members)
            c = self.strict_fn(R.detach().numpy(), self.members)
            return torch.as_tensor(float(np.mean(c)), dtype=DT)
        rows = self.model.predict(R, self.V, self.strain)
        x, a = pack_predicted(self.region, rows, self.bary, self.model.n_det)
        if sp.kind == "sinkhorn" and sp.dilate > 0:
            off = torch.tensor([[1, 1], [1, -1], [-1, 1], [-1, -1]], dtype=DT) * sp.dilate
            xo = x[:, :, None, :].repeat(1, 1, 4, 1)
            xo[..., :2] = xo[..., :2] + off
            x = xo.reshape(x.shape[0], -1, 3)
            a = (a[:, :, None] / 4.0).expand(-1, -1, 4).reshape(a.shape[0], -1)
        if sp.kind == "sinkhorn":
            S, P, O = x.shape[0], x.shape[1], self.y.shape[1]
            step = max(1, self.chunk_elems // max(P * max(O, P), 1))
            parts = [
                sinkhorn_divergence(
                    x[i : i + step], a[i : i + step], self.y[i : i + step], self.b[i : i + step],
                    blur_px=sp.blur_px, lam=sp.lam, diameter_px=sp.diameter_px,
                    reach_px=sp.reach_px,
                )
                for i in range(0, S, step)
            ]
            v = torch.cat(parts)
        elif sp.kind == "ls":
            if self._obs_grids is None:
                self._obs_grids = observed_grids(self.region, intensity=True)
            v = ls_objective(soft_raster(self.region, x, a, sp.sigma_px, sp.tau_frames), self._obs_grids)
        else:
            raise ValueError(sp.kind)
        return v.mean()

    def value(self, R: np.ndarray, spec: Optional[ObjectiveSpec] = None) -> float:
        with torch.no_grad():
            return float(self.objective(torch.as_tensor(R, dtype=DT), spec))

    def stage_a(
        self,
        sigma0_deg: float = 0.05,
        max_evals: int = 200,
        seed: int = 0,
        theta0: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray, Dict[str, Any]]:
        """CMA over theta (grain rotation + field). Returns (theta, member orientations, info)."""
        import cma

        th0 = self.theta0 if theta0 is None else theta0
        e0 = self.n_evals

        def f(th: np.ndarray) -> float:
            with torch.no_grad():
                return float(self.objective(self.field.orientations(torch.as_tensor(th, dtype=DT))))

        f_start = f(th0)
        es = cma.CMAEvolutionStrategy(
            list(th0), math.radians(sigma0_deg),
            {"seed": seed + 1, "maxfevals": max_evals, "verbose": -9, "tolfun": 1e-6, "tolx": 1e-7},
        )
        best_th, best_f = np.asarray(th0, float), f_start
        while not es.stop():
            X = es.ask()
            F = [f(np.asarray(x)) for x in X]
            es.tell(X, F)
            i = int(np.argmin(F))
            if F[i] < best_f:
                best_f, best_th = F[i], np.asarray(X[i], float)
        with torch.no_grad():
            R = self.field.orientations(torch.as_tensor(best_th, dtype=DT)).numpy()
        info = dict(f_start=f_start, f_end=best_f, evals=self.n_evals - e0,
                    stop=str(dict(es.stop())))
        return best_th, R, info

    def stage_b(
        self,
        R_init: np.ndarray,
        steps: int = 60,
        lr_deg: float = 0.01,
        mu: float = 1.0,
        neighbours: Optional[Sequence[Sequence[int]]] = None,
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """Adam on a free rotation vector per voxel, ``R_v = exp([d_v]) R_init_v``, minimising
        objective + ``mu`` * mean over neighbour pairs of |d_u - d_v|^2 / (0.01 deg)^2 * 1e-2
        (the penalty is in units where a 0.01 deg step between neighbours costs ``mu`` * 1e-2).
        Returns (orientations, info); the best iterate by the full loss is kept."""
        n = len(self.members)
        pos = {v: i for i, v in enumerate(self.members)}
        pairs: List[Tuple[int, int]] = []
        if neighbours is not None:
            for v in self.members:
                for u in neighbours[v]:
                    if int(u) in pos and pos[int(u)] > pos[v]:
                        pairs.append((pos[v], pos[int(u)]))
        P = torch.as_tensor(pairs, dtype=torch.long).reshape(-1, 2)
        R0 = torch.as_tensor(R_init, dtype=DT)
        d = torch.zeros((n, 3), dtype=DT, requires_grad=True)
        opt = torch.optim.Adam([d], lr=math.radians(lr_deg))
        unit = math.radians(0.01) ** 2
        e0 = self.n_evals

        def loss_of(dd: Tensor) -> Tuple[Tensor, Tensor]:
            R = rotvec_to_matrix(dd) @ R0
            obj = self.objective(R)
            if len(P):
                smooth = ((dd[P[:, 0]] - dd[P[:, 1]]) ** 2).sum(-1).mean() / unit * 1e-2
            else:
                smooth = torch.zeros((), dtype=DT)
            return obj, smooth

        hist = []
        best = (float("inf"), d.detach().clone(), 0.0, 0.0)
        for it in range(steps):
            opt.zero_grad()
            obj, smooth = loss_of(d)
            loss = obj + mu * smooth
            lv = float(loss)
            if lv < best[0]:
                best = (lv, d.detach().clone(), float(obj), float(smooth))
            loss.backward()
            opt.step()
            hist.append((float(obj), float(smooth)))
        with torch.no_grad():
            obj, smooth = loss_of(d)
            lv = float(obj + mu * smooth)
            if lv < best[0]:
                best = (lv, d.detach().clone(), float(obj), float(smooth))
            R = (rotvec_to_matrix(best[1]) @ R0).numpy()
        info = dict(
            f_start=hist[0][0] if hist else float("nan"), f_end=best[2], smooth_end=best[3],
            evals=self.n_evals - e0, n_pairs=len(pairs),
            max_offset_deg=float(np.degrees(best[1].norm(dim=1).max())) if n else 0.0,
        )
        return R, info
