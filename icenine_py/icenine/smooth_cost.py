"""Smooth surrogate costs for ranking candidates at the coarse search levels (Task E, opt-in).

The strict cost (``cost_functions.VoxelCostFunction``) is a ratio of pixel counts, so it changes
in steps: a predicted spot moves to the neighbouring frame (image), a vertex of the predicted
triangle moves to another pixel, a peak enters or leaves the evaluated set. This module builds
costs that change gradually with the orientation instead. They are used to RANK candidates at the
coarse levels (``AdaptiveVoxelReconstructor`` with config key ``SmoothCostRank``); the strict cost
stays the final arbiter everywhere else (the two-costs-by-design rule of the ray rule applies).

Arms (combine with ``+``; ``S`` alone is the control: binary credit, hard peak set):

* **E1, distance-transform overlap.** Each predicted spot is sampled at its three projected
  vertices and its centroid. A sample point earns credit ``exp(-d^2 / 2 sigma^2)`` where ``d`` is
  its distance to the nearest lit pixel of the frame (Chebyshev distance to the lit pixel's square,
  0 inside a lit pixel). Without E1 the credit is binary (1 inside a lit pixel, else 0), which is a
  sampled version of the strict pixel overlap.
* **E2, omega blending.** A spot also takes credit from the neighbouring frames, weighted by the
  distance of its predicted omega to each frame's omega interval: weight
  ``1 - smoothstep(d / omega_half)`` (``d`` in frames; 1 for the spot's own frame, 1 at the
  boundary of the neighbour, 0 at ``omega_half`` frames past it). The credit is the largest
  weighted frame credit, so crossing a frame boundary does not switch images abruptly, and the true
  orientation keeps full credit on images rendered with one frame per spot (an overlap-weighted
  mean over frames would not: it was tried first and raised the cost at the truth from 0.053 to
  0.296 on voxel 15901).
* **E3, denominator stabilisation.** Each peak carries a weight that ramps smoothly to zero as its
  spot approaches a detector edge, as its eta approaches the eta limit and as its omega approaches
  either end of the scanned range. The quality is the weighted mean of the per-peak credits, so the
  ratio does not jump when the number of peaks changes. The per-peak credit is the mean over the
  detectors the spot reaches (weighted the same way), not over all detectors, which removes the
  fixed 1/2 of a peak whose spot cannot reach detector 1.
* **E4, observed-ray residual.** Spots are extracted from the images, paired across the two
  detectors (``spot_pairing.pair_frame``, the ``ring_tight`` preset of Task D) and attached to the
  voxel with Task D's footprint rule. Per predicted peak the residual is the angle between the
  predicted scattered ray and the nearest attached observed ray, plus the distance of the predicted
  omega to that ray's frame interval; the cost is the mean of ``r^2 / (r^2 + tau^2)``. Continuous
  in the orientation, with an autograd gradient (``SmoothCost.e4_cost_torch``).

Without E4 the cost is ``1 - quality`` with quality = mean per-peak credit (E1/E2/E3 as set);
with E4 alone it is the residual cost; with E4 and others the two are averaged.

Spec strings (config key ``SmoothCostRank``, no spaces): ``off`` or the arms followed by optional
``key=value`` items, separated by commas, e.g. ``E1+E2,sigma=8`` or ``E4,tau=0.25,levels=2``.
Keys: ``sigma`` (E1 kernel width in pixels; ``auto`` = from the level, see ``level_widths``),
``omega`` (E2 kernel reach in frames past the boundary; ``auto``), ``tau`` (E4 scale in degrees;
``auto``), ``levels`` (the ranking uses the smooth cost at levels ``0 .. levels-1``, default 3; with
fewer search levels than that the last level, and so the hand-off order, is ranked too), ``ramp``
(E3 detector ramp in pixels), ``ramp_eta`` (degrees), ``ramp_omega`` (frames).

Frame and beam assumption of E4 (as ``spot_pairing``): lab = Rz(omega) sample, beam along +x, no
deflection. E1 to E3 use the general sample-to-lab transform of the strict cost.
"""

import math
from dataclasses import dataclass, field
from typing import Any, Dict, FrozenSet, List, Optional, Sequence, Tuple

import numpy as np
import torch
from scipy import ndimage
from scipy.spatial import cKDTree

from .cost_functions import FrameTable, frame_index_array
from .diffraction_core import BraggConstants, scattering_omegas_fixed_beam

ARMS = ("S", "E1", "E2", "E3", "E4")  # S: sampled binary overlap (control, no smoothing)
FRAME_SPACING = 1.0e5  # frame coordinate scale of the lit-pixel tree (far above any pixel distance)
N_SAMPLES = 4  # sample points per spot: three vertices and the centroid


# --------------------------------------------------------------------------------------------
# spec
# --------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class SmoothSpec:
    """Parsed ``SmoothCostRank`` spec. A width of 0 means ``auto`` (taken from the level)."""

    arms: FrozenSet[str] = frozenset()
    sigma_px: float = 0.0
    omega_half: float = 0.0
    tau_deg: float = 0.0
    levels: int = 3
    ramp_px: float = 10.0
    ramp_eta_deg: float = 1.0
    ramp_omega: float = 0.5

    @property
    def active(self) -> bool:
        return bool(self.arms)

    def label(self) -> str:
        return "+".join(a for a in ARMS if a in self.arms) or "off"

    @classmethod
    def parse(cls, spec: str) -> "SmoothSpec":
        s = (spec or "off").strip()
        if s.lower() in ("off", "0", "none", ""):
            return cls()
        if " " in s:
            raise ValueError(f"SmoothCostRank spec must not contain spaces: {spec!r}")
        head, *items = s.split(",")
        arms = frozenset(a.strip().upper() for a in head.split("+") if a.strip())
        bad = arms - set(ARMS)
        if bad or not arms:
            raise ValueError(
                f"SmoothCostRank: unknown arm(s) {sorted(bad)} in {spec!r}; use {ARMS}"
            )
        kw: Dict[str, Any] = {}
        names = dict(
            sigma="sigma_px",
            omega="omega_half",
            tau="tau_deg",
            levels="levels",
            ramp="ramp_px",
            ramp_eta="ramp_eta_deg",
            ramp_omega="ramp_omega",
        )
        for it in items:
            if "=" not in it:
                raise ValueError(f"SmoothCostRank: item {it!r} is not key=value in {spec!r}")
            k, v = it.split("=", 1)
            if k not in names:
                raise ValueError(f"SmoothCostRank: unknown key {k!r}; use {sorted(names)}")
            if k == "levels":
                kw[names[k]] = int(v)
                if kw[names[k]] < 1:
                    raise ValueError("SmoothCostRank: levels must be >= 1")
            elif v == "auto":
                if k not in ("sigma", "omega", "tau"):
                    raise ValueError(f"SmoothCostRank: {k} has no auto value")
                kw[names[k]] = 0.0
            else:
                x = float(v)
                if not x > 0:
                    raise ValueError(f"SmoothCostRank: {k} must be > 0, got {v}")
                kw[names[k]] = x
        return cls(arms=arms, **kw)


# --------------------------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------------------------
def smoothstep(x: torch.Tensor) -> torch.Tensor:
    """0 below 0, 1 above 1, 3x^2 - 2x^3 in between (continuous with a continuous slope)."""
    t = torch.clamp(x, 0.0, 1.0)
    out: torch.Tensor = t * t * (3.0 - 2.0 * t)
    return out


def sample_hull_points(sample: Any) -> np.ndarray:
    """(n, 2) triangle vertices of every voxel of the sample's mic (sample frame, mm)."""
    vox = sample.get_mic().voxels
    s3 = math.sqrt(3.0)
    pts = []
    for v in vox:
        x, y = float(v.position[0]), float(v.position[1])
        s = float(v.side_length)
        if v.points_up:
            pts += [(x, y), (x + s, y), (x + s / 2, y + s / 2 * s3)]
        else:
            pts += [(x, y), (x + s / 2, y - s / 2 * s3), (x + s, y)]
    return np.asarray(pts, dtype=np.float64)


@dataclass
class Prediction:
    """Per-peak geometry of a batch of orientations (only the rows that pass Bragg + eta).

    ``b`` batch index, ``refl``/``branch`` the reflection and omega branch, ``omega`` (rad),
    ``fpos`` continuous frame position (omega - low) / width, ``bin`` the omega bin (truncation as
    the strict cost), ``frame`` the frame of that bin (-1: none), ``eta`` (rad), ``ray`` (n, 3)
    unit scattered direction, ``col``/``row`` (n, D, S) continuous pixel coordinates of the sample
    points on each detector, ``hit`` (n, D, S) the ray reaches the detector plane."""

    b: torch.Tensor
    refl: torch.Tensor
    branch: torch.Tensor
    omega: torch.Tensor
    fpos: torch.Tensor
    bin: torch.Tensor
    frame: torch.Tensor
    eta: torch.Tensor
    ray: torch.Tensor
    col: torch.Tensor
    row: torch.Tensor
    hit: torch.Tensor
    n_batch: int


@dataclass
class ObservedRays:
    """Paired observed spots (one per ray): frame, unit direction (lab), and both spots' boxes."""

    frame: np.ndarray
    dirn: np.ndarray
    box0: np.ndarray  # (n, 4) cmin, cmax, rmin, rmax on detector 0 (pixel indices)
    box1: np.ndarray
    n_spots: Tuple[int, int] = (0, 0)
    n_candidates: int = 0


# --------------------------------------------------------------------------------------------
# the evaluator
# --------------------------------------------------------------------------------------------
class SmoothCost:
    """Batched smooth costs of one image set.

    Built from the same objects as ``VoxelCostFunction``. The lit-pixel trees are built on first
    use from ``exp_data`` (through the S1 ``FrameTable``); the observed rays (E4) on the first E4
    call. The images must not change afterwards (``FrameTable`` contract); build a new
    ``SmoothCost`` for new images.
    """

    def __init__(
        self,
        simulator: Any,
        detector_list: Sequence[Any],
        range_map: Any,
        exp_data: Any,
        sample: Any,
        structure_list: Sequence[Any],
        eta_limit: float = math.pi / 2.0,
        max_q: float = 0.0,
        min_sin_eta: float = 0.0,
        hull_points: Optional[np.ndarray] = None,
    ) -> None:
        from .spot_pairing import geometry_from_detectors

        self.simulator = simulator
        self.detector_list = list(detector_list)
        self.range_map = range_map
        self.exp_data = exp_data
        self.sample = sample
        self.eta_limit = float(eta_limit)
        self.min_sin_eta = float(min_sin_eta)
        self.n_det = len(self.detector_list)
        self.dtype = torch.float64
        d = self.dtype
        self.geoms = geometry_from_detectors(self.detector_list)
        self._frame_index = frame_index_array(range_map)
        self.n_bins = len(self._frame_index) - 1
        self.n_frames = int(self._frame_index.max()) + 1
        self.low = float(range_map.low)
        self.width = float(range_map.width)
        assert self.width > 0, "SmoothCost assumes increasing omega bins"
        self.index_t = torch.as_tensor(self._frame_index, dtype=torch.long)
        m = sample.sample_to_lab_matrix.to(d)
        self.base = m[:3, :3].clone()
        self.trans = m[:3, 3].clone()
        self.beam = torch.as_tensor(np.asarray(simulator.beam_direction, dtype=np.float64))
        self.beam = self.beam / torch.linalg.norm(self.beam)
        self.energy = float(simulator.beam_energy)
        self.chi = float(simulator.beam_deflection_chi)
        self._g: Dict[int, Tuple[torch.Tensor, torch.Tensor]] = {}
        for p, st in enumerate(structure_list):
            rv = st.get_reflection_vectors()
            if max_q > 0:
                rv = [r for r in rv if r.q_mag <= max_q]
            if rv:
                self._g[p] = (
                    torch.tensor(np.stack([r.q_vec for r in rv]), dtype=d),
                    torch.tensor([r.q_mag for r in rv], dtype=d),
                )
        self._consts: Dict[Tuple[int, int], BraggConstants] = {}
        self._trees: Optional[List[Optional[cKDTree]]] = None
        self.rows = 0
        self.cols = 0
        self.n_lit = [0] * self.n_det
        self._rays: Optional[ObservedRays] = None
        self._hull_points = hull_points
        self._attached: Dict[Tuple[float, float], np.ndarray] = {}
        # detector constants (float64 torch)
        self._det: List[Dict[str, Any]] = [
            dict(
                n=torch.as_tensor(g.n, dtype=d),
                d=float(g.d),
                off=torch.as_tensor(g.pos + g.org, dtype=d),
                bj=torch.as_tensor(g.bj, dtype=d),
                bk=torch.as_tensor(g.bk, dtype=d),
                hw=float(g.hw),
                hh=float(g.hh),
                pw=float(g.pw),
                ph=float(g.ph),
            )
            for g in self.geoms
        ]

    # ---------------------------------------------------------------- geometry
    def px_per_rad(self) -> float:
        """Pixels per radian of ray turn on detector 0 (distance / pixel width)."""
        g = self.geoms[0]
        return float(g.distance_x() / g.pw)

    def predict(self, R: Any, voxel_vertices: Any, phase: int = 0) -> Prediction:
        """Peak geometry of orientations ``R`` (B, 3, 3) (numpy or torch) for one voxel
        (differentiable in a torch ``R``)."""
        d = self.dtype
        Rt = R if isinstance(R, torch.Tensor) else torch.as_tensor(np.asarray(R), dtype=d)
        Rt = Rt.to(d)
        if Rt.dim() == 2:
            Rt = Rt[None]
        B = int(Rt.shape[0])
        g_hkl, g_mag = self._g[phase]
        K = int(g_hkl.shape[0])
        key = (phase, B)
        consts = self._consts.get(key)
        if consts is None:
            consts = self._consts[key] = BraggConstants(g_mag.repeat(B), self.energy, self.chi)
        g_lab = torch.einsum("bij,kj->bki", Rt, g_hkl).reshape(B * K, 3)
        br = scattering_omegas_fixed_beam(g_lab, consts)
        om = torch.stack([br.omega1, br.omega2], dim=1)  # (BK, 2)
        obs = br.observable[:, None].expand(B * K, 2)
        nrm = g_lab / torch.linalg.norm(g_lab, dim=1, keepdim=True)
        nrm = (self.base @ nrm.T).T  # base rotation
        c, s = torch.cos(om), torch.sin(om)
        nx, ny, nz = nrm[:, 0:1], nrm[:, 1:2], nrm[:, 2:3]
        ln = torch.stack([c * nx - s * ny, s * nx + c * ny, nz.expand_as(c)], dim=-1)  # (BK,2,3)
        dot = (ln * self.beam).sum(-1, keepdim=True)
        ray = self.beam - 2.0 * dot * ln  # (BK, 2, 3)
        rn = torch.linalg.norm(ray, dim=-1)
        eta = torch.atan2(ray[..., 1].abs() / rn, ray[..., 2].abs() / rn)
        valid = obs & (eta < self.eta_limit) & (rn > 0)
        if self.min_sin_eta > 0:
            valid = valid & (torch.sin(eta) >= self.min_sin_eta)
        sel = torch.nonzero(valid.reshape(-1))[:, 0]
        bk = sel // 2
        branch = sel % 2
        b = bk // K
        refl = bk % K
        omega = om.reshape(-1)[sel]
        ray_s = (ray.reshape(-1, 3)[sel]) / rn.reshape(-1)[sel, None]
        eta_s = eta.reshape(-1)[sel]
        fpos = (omega - self.low) / self.width
        # truncation toward zero, as the strict cost's lenient path (``frame_of_omegas``); with
        # CppEdgeRules on the strict cost floors instead, which differs only for an omega just
        # below the scanned range (fpos in (-1, 0): truncation gives bin 0, floor gives -1)
        fbin = torch.trunc(fpos).to(torch.long)
        fbin_c = fbin.clamp(-1, self.n_bins)
        frame = self.index_t[fbin_c]
        # sample points: three vertices + centroid, sample frame -> lab (Rz(omega) @ base @ v + t)
        V = voxel_vertices
        V = V if isinstance(V, torch.Tensor) else torch.as_tensor(np.asarray(V))
        V = V.to(d)
        P = torch.cat([V, V.mean(0, keepdim=True)], dim=0)  # (S, 3)
        P = (self.base @ P.T).T  # (S, 3)
        cw, sw = torch.cos(omega)[:, None], torch.sin(omega)[:, None]
        px = cw * P[None, :, 0] - sw * P[None, :, 1] + self.trans[0]
        py = sw * P[None, :, 0] + cw * P[None, :, 1] + self.trans[1]
        pz = P[None, :, 2].expand_as(px) + self.trans[2]
        pts = torch.stack([px, py, pz], dim=-1)  # (n, S, 3)
        cols, rows, hits = [], [], []
        for g in self._det:
            den = ray_s @ g["n"]  # (n,)
            num = -((pts @ g["n"]) + g["d"])  # (n, S)
            par = den.abs() < 1e-12
            t = num / torch.where(par, torch.ones_like(den), den)[:, None]
            hit = (~par)[:, None] & (t > 0)
            q = pts + t[..., None] * ray_s[:, None, :] - g["off"]
            cols.append(((q @ g["bj"]) + g["hw"]) / g["pw"])
            rows.append(((q @ g["bk"]) + g["hh"]) / g["ph"])
            hits.append(hit)
        return Prediction(
            b=b,
            refl=refl,
            branch=branch,
            omega=omega,
            fpos=fpos,
            bin=fbin,
            frame=frame,
            eta=eta_s,
            ray=ray_s,
            col=torch.stack(cols, dim=1),
            row=torch.stack(rows, dim=1),
            hit=torch.stack(hits, dim=1),
            n_batch=B,
        )

    # ---------------------------------------------------------------- lit-pixel trees
    def _ensure_trees(self) -> None:
        if self._trees is not None:
            return
        table = FrameTable(self.exp_data, self.n_frames, self.n_det)
        table.ensure(np.arange(self.n_frames))
        self.rows, self.cols = table.rows, table.cols
        trees: List[Optional[cKDTree]] = []
        for det in range(self.n_det):
            parts = []
            for w in range(self.n_frames):
                img = table.images[w * self.n_det + det]
                if not img.any():
                    continue
                r, c = np.nonzero(img)
                if len(r):
                    parts.append(
                        np.stack(
                            [c + 0.5, r + 0.5, np.full(len(r), w * FRAME_SPACING)], axis=1
                        ).astype(np.float64)
                    )
            if parts:
                pts = np.concatenate(parts)
                self.n_lit[det] = len(pts)
                trees.append(cKDTree(pts))
            else:
                trees.append(None)
        self._table = table
        self._trees = trees

    def _distance(
        self, det: int, col: np.ndarray, row: np.ndarray, frame: np.ndarray, ub: float
    ) -> np.ndarray:
        """Chebyshev distance (pixels) of each point to the centre of the nearest lit pixel of its
        frame (inf beyond ``ub`` or with no frame)."""
        self._ensure_trees()
        out = np.full(col.shape, np.inf)
        tree = self._trees[det]  # type: ignore[index]
        ok = frame >= 0
        if tree is None or not ok.any():
            return out
        q = np.stack([col[ok], row[ok], frame[ok] * FRAME_SPACING], axis=-1)
        dist, _ = tree.query(q, k=1, p=np.inf, distance_upper_bound=ub)
        out[ok] = dist
        return out

    # ---------------------------------------------------------------- E1 / E2 / E3 quality
    def quality_terms(
        self,
        pr: Prediction,
        spec: SmoothSpec,
        sigma_px: float,
        omega_half: float,
        dcache: Optional[Dict[Any, Any]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """(per-peak weight, per-peak credit) of the image arms (numpy, length n).

        ``dcache`` (optional, shared by several specs on the same ``pr``): the nearest-lit
        distances per (detector, frame offset), queried once with the bound ``dcache["ub"]``,
        which must be at least this spec's bound; the result equals a call without it."""
        self._ensure_trees()
        n = int(pr.b.shape[0])
        D = self.n_det
        col = pr.col.detach().numpy()
        row = pr.row.detach().numpy()
        hit = pr.hit.detach().numpy()
        fpos = pr.fpos.detach().numpy()
        fbin = pr.bin.detach().numpy()
        e1 = "E1" in spec.arms
        ub = 0.5 + 1e-9 if not e1 else 0.5 + 4.0 * sigma_px
        inb = hit & (col >= 0) & (col < self.cols) & (row >= 0) & (row < self.rows)  # (n, D, S)
        # frame offsets and their kernel (E2: credit of frame k' times K(distance of the predicted
        # omega to that frame's interval), K = 1 - smoothstep(d / h); the max over frames is kept)
        if "E2" in spec.arms:
            h = float(omega_half)
            kmax = int(math.ceil(h))
            offs = list(range(-kmax, kmax + 1))
            wts = []
            for o in offs:
                lo = fbin + o
                dist = np.maximum(np.maximum(lo - fpos, fpos - (lo + 1.0)), 0.0)
                wts.append(1.0 - smoothstep(torch.as_tensor(dist / h)).numpy())
        else:
            offs, wts = [0], [np.ones(n)]
        credit_d = np.zeros((n, D))
        for o, w in zip(offs, wts):
            use = w > 0
            if not use.any():
                continue
            bo = fbin + o
            fr = np.where(
                (bo >= 0) & (bo < self.n_bins), self._frame_index[bo.clip(0, self.n_bins)], -1
            )
            for det in range(D):
                c = col[:, det, :]
                r = row[:, det, :]
                if dcache is not None and dcache["ub"] >= ub:
                    ck = (det, o)
                    if ck not in dcache:
                        fr_all = np.broadcast_to(fr[:, None], c.shape)
                        dcache[ck] = self._distance(
                            det, c.ravel(), r.ravel(), fr_all.ravel(), float(dcache["ub"])
                        ).reshape(c.shape)
                    dist = np.where(use[:, None], dcache[ck], np.inf)
                else:
                    fr_s = np.broadcast_to(np.where(use, fr, -1)[:, None], c.shape)
                    dist = self._distance(det, c.ravel(), r.ravel(), fr_s.ravel(), ub).reshape(
                        c.shape
                    )
                if e1:
                    ex = np.clip(dist - 0.5, 0.0, None)
                    F = np.where(dist <= ub, np.exp(-(ex * ex) / (2.0 * sigma_px * sigma_px)), 0.0)
                else:
                    F = (dist <= 0.5).astype(np.float64)
                m = inb[:, det, :]
                cnt = m.sum(1)
                cr = np.where(cnt > 0, (F * m).sum(1) / np.maximum(cnt, 1), 0.0)
                credit_d[:, det] = np.maximum(credit_d[:, det], w * cr)
        if "E3" in spec.arms:
            # presence weights: detector margin of the centroid, eta and omega margins
            cc = col[:, :, -1]
            rr = row[:, :, -1]
            margin = np.minimum.reduce([cc, self.cols - cc, rr, self.rows - rr])
            wd = smoothstep(torch.as_tensor(margin / spec.ramp_px)).numpy() * hit[:, :, -1]
            eta = pr.eta.detach().numpy()
            we = smoothstep(
                torch.as_tensor((self.eta_limit - eta) / math.radians(spec.ramp_eta_deg))
            ).numpy()
            wo = smoothstep(
                torch.as_tensor(np.minimum(fpos, self.n_bins - fpos) / spec.ramp_omega)
            ).numpy()
            sw = wd.sum(1)
            q = np.where(sw > 0, (wd * credit_d).sum(1) / np.maximum(sw, 1e-300), 0.0)
            W = we * wo * wd.max(1)
            return W, q
        present = (inb.any(2).any(1)) & (pr.frame.numpy() >= 0)
        q = credit_d.sum(1) / D
        return present.astype(np.float64), q

    # ---------------------------------------------------------------- E4
    def observed_rays(self) -> ObservedRays:
        """Spots of every frame (8-connected components of the binary images, binary centroids),
        paired across the two detectors with Task D's ``ring_tight`` preset."""
        if self._rays is not None:
            return self._rays
        from .constants import KEV_OVER_HBAR_C_IN_ANG
        from . import spot_pairing as sp

        self._check_e4_frame()
        self._ensure_trees()
        tab = self._table
        hull_pts = (
            self._hull_points if self._hull_points is not None else sample_hull_points(self.sample)
        )
        hull = sp.SampleHull(hull_pts)
        side = self._voxel_side()
        k = KEV_OVER_HBAR_C_IN_ANG * self.energy
        qs = np.concatenate([g[1].numpy() for g in self._g.values()])
        rings = np.unique(np.round(2.0 * np.arcsin(qs / (2.0 * k)), 9))
        side_px = side / self.geoms[0].pw
        params = sp.PairingParams(
            centroid_tol_px=1.0,
            registration_mm=0.0,
            voxel_size_mm=side,
            omega_step=self.width,
            voxel_term=False,
            arc_term=True,
            ring_two_theta=rings,
            multi_npix=int(2.5 * 0.433 * side_px**2),
            max_partners=1,
            multi_partners=3,
        )
        struct = np.ones((3, 3), dtype=bool)
        fr_l, dir_l, b0_l, b1_l = [], [], [], []
        n_sp = [0, 0]
        n_cand = 0
        for w in range(self.n_frames):
            bins = np.nonzero(self._frame_index[:-1] == w)[0]
            if len(bins) == 0:
                continue
            om_c = self.low + self.width * (0.5 * (bins.min() + bins.max() + 1))
            comps: List[Optional[Dict[str, Any]]] = []
            for det in (0, 1):
                img = tab.images[w * self.n_det + det]
                if not img.any():
                    comps.append(None)
                    continue
                lab, nlab = ndimage.label(img, structure=struct)
                n_sp[det] += nlab
                if nlab == 0:
                    comps.append(None)
                    continue
                rr, cc = np.nonzero(lab)
                ll = lab[rr, cc]
                npix = np.bincount(ll, minlength=nlab + 1)[1:].astype(np.float64)
                cm = np.stack(
                    [
                        np.bincount(ll, weights=rr, minlength=nlab + 1)[1:] / npix,
                        np.bincount(ll, weights=cc, minlength=nlab + 1)[1:] / npix,
                    ],
                    axis=1,
                )
                sl = ndimage.find_objects(lab)
                box = np.array(
                    [[s[1].start, s[1].stop - 1, s[0].start, s[0].stop - 1] for s in sl], float
                )
                edge = (
                    (box[:, 0] <= 0)
                    | (box[:, 2] <= 0)
                    | (box[:, 1] >= img.shape[1] - 1)
                    | (box[:, 3] >= img.shape[0] - 1)
                )
                comps.append(
                    dict(
                        c_col=cm[:, 1],
                        c_row=cm[:, 0],
                        npix=npix,
                        edge=edge,
                        w=box[:, 1] - box[:, 0],
                        h=box[:, 3] - box[:, 2],
                        box=box,
                    )
                )
            if comps[0] is None or comps[1] is None:
                continue
            r = sp.pair_frame(self.geoms, hull, float(om_c), comps[0], comps[1], params)
            n_cand += len(r["i0"])
            ch = r["chosen"].astype(bool)
            i0 = r["i0"][ch].astype(int)
            i1 = r["i1"][ch].astype(int)
            fr_l.append(np.full(len(i0), w))
            dir_l.append(r["dirn"][ch])
            b0_l.append(comps[0]["box"][i0])
            b1_l.append(comps[1]["box"][i1])
        cat = lambda L, s: np.concatenate(L) if L else np.zeros(s)  # noqa: E731
        self._rays = ObservedRays(
            frame=cat(fr_l, (0,)).astype(np.int64),
            dirn=cat(dir_l, (0, 3)),
            box0=cat(b0_l, (0, 4)),
            box1=cat(b1_l, (0, 4)),
            n_spots=(n_sp[0], n_sp[1]),
            n_candidates=n_cand,
        )
        return self._rays

    def _check_e4_frame(self) -> None:
        m = np.asarray(self.sample.sample_to_lab_matrix.double().numpy())
        if not np.allclose(m, np.eye(4), atol=1e-9):
            raise ValueError("E4 needs lab = Rz(omega) sample (identity sample-to-lab transform)")
        if not np.allclose(self.beam.numpy(), [1.0, 0.0, 0.0], atol=1e-9) or self.chi != 0.0:
            raise ValueError("E4 needs the beam along +x without deflection")
        if self.n_det != 2:
            raise ValueError("E4 needs two detectors")

    def _voxel_side(self) -> float:
        vox = self.sample.get_mic().voxels
        return float(vox[0].side_length)

    def attached_rays(self, voxel_vertices: Any, margin_px: float = 1.0) -> np.ndarray:
        """Indices of the observed rays attached to the voxel by Task D's footprint rule: the line
        from the voxel centroid (at the frame-centre omega) along the ray's direction hits the
        bounding boxes of both its spots (plus ``margin_px`` and the arc motion r x d_omega / 2)."""
        from . import spot_pairing as sp

        V = np.asarray(
            (
                voxel_vertices.detach().numpy()
                if isinstance(voxel_vertices, torch.Tensor)
                else voxel_vertices
            ),
            dtype=np.float64,
        )
        cen = V.mean(0)[:2]
        key = (round(float(cen[0]), 9), round(float(cen[1]), 9))
        if key in self._attached:
            return self._attached[key]
        rays = self.observed_rays()
        if len(rays.frame) == 0:
            self._attached[key] = np.zeros(0, dtype=np.int64)
            return self._attached[key]
        bins_c = np.array(
            [np.nonzero(self._frame_index[:-1] == w)[0].mean() + 0.5 for w in range(self.n_frames)]
        )
        om = self.low + self.width * bins_c[rays.frame]
        o = np.zeros((len(om), 3))
        o[:, :2] = sp.sample_to_lab_xy(np.broadcast_to(cen, (len(om), 2)), om)
        rad = float(np.linalg.norm(cen))
        margin = margin_px + rad * self.width / 2.0 / self.geoms[0].pw
        ok = np.ones(len(om), dtype=bool)
        for det, box in ((0, rays.box0), (1, rays.box1)):
            g = self.geoms[det]
            den = rays.dirn @ g.n
            t = -(o @ g.n + g.d) / np.where(np.abs(den) > 1e-12, den, 1e-12)
            pt = o + t[:, None] * rays.dirn
            c, r = sp.lab_to_pixel(g, pt)
            ok &= (
                (c >= box[:, 0] - margin)
                & (c <= box[:, 1] + 1 + margin)
                & (r >= box[:, 2] - margin)
                & (r <= box[:, 3] + 1 + margin)
            )
        idx = np.nonzero(ok)[0]
        self._attached[key] = idx
        return idx

    def _e4_terms(
        self, pr: Prediction, voxel_vertices: Any, tau_rad: float, spec: SmoothSpec
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """(per-peak weight, per-peak robust residual) of E4 as torch tensors (differentiable in
        the orientation through ``pr``)."""
        rays = self.observed_rays()
        att = self.attached_rays(voxel_vertices)
        n = int(pr.b.shape[0])
        # predicted peaks seen on both detectors (centroid sample point in bounds)
        col = pr.col[:, :, -1]
        row = pr.row[:, :, -1]
        hit = pr.hit[:, :, -1]
        if "E3" in spec.arms:
            margin = torch.minimum(
                torch.minimum(col, self.cols - col), torch.minimum(row, self.rows - row)
            )
            wd = smoothstep(margin / spec.ramp_px) * hit.to(self.dtype)
            we = smoothstep((self.eta_limit - pr.eta) / math.radians(spec.ramp_eta_deg))
            wo = smoothstep(torch.minimum(pr.fpos, self.n_bins - pr.fpos) / spec.ramp_omega)
            W = (wd.min(1).values * we * wo).detach()
        else:
            inb = hit & (col >= 0) & (col < self.cols) & (row >= 0) & (row < self.rows)
            W = (inb.all(1) & (pr.frame >= 0)).to(self.dtype)
        if len(att) == 0 or n == 0:
            return W, torch.ones(n, dtype=self.dtype)
        fa = torch.as_tensor(rays.frame[att])
        da = torch.as_tensor(rays.dirn[att], dtype=self.dtype)  # (A, 3)
        # omega interval of each attached ray's frame
        lo_bin = np.array(
            [np.nonzero(self._frame_index[:-1] == w)[0].min() for w in range(self.n_frames)]
        )
        hi_bin = np.array(
            [np.nonzero(self._frame_index[:-1] == w)[0].max() + 1 for w in range(self.n_frames)]
        )
        lo = torch.as_tensor(self.low + self.width * lo_bin[rays.frame[att]], dtype=self.dtype)
        hi = torch.as_tensor(self.low + self.width * hi_bin[rays.frame[att]], dtype=self.dtype)
        chunk = max(1, 2_000_000 // max(len(att), 1))
        out = []
        for s0 in range(0, n, chunk):
            rs = pr.ray[s0 : s0 + chunk]
            om = pr.omega[s0 : s0 + chunk]
            cosang = torch.clamp(rs @ da.T, -1.0, 1.0)
            # angle via atan2 of |cross| and dot: accurate and differentiable near 0
            cross = torch.linalg.norm(
                torch.cross(
                    rs[:, None, :].expand(-1, len(att), -1),
                    da[None].expand(len(rs), -1, -1),
                    dim=-1,
                ),
                dim=-1,
            )
            ang = torch.atan2(cross, cosang)
            dom = torch.clamp(
                torch.maximum(lo[None] - om[:, None], om[:, None] - hi[None]), min=0.0
            )
            r2 = ang * ang + dom * dom
            near = (fa[None] - pr.frame[s0 : s0 + chunk, None]).abs() <= 1
            r2 = torch.where(near, r2, torch.full_like(r2, float("inf")))
            r2min = r2.min(dim=1).values
            out.append(
                torch.where(
                    torch.isfinite(r2min), r2min / (r2min + tau_rad**2), torch.ones_like(r2min)
                )
            )
        rho = torch.cat(out)
        return W, rho

    # ---------------------------------------------------------------- public costs
    def level_widths(
        self, spec: SmoothSpec, diameter_rad: Optional[float]
    ) -> Tuple[float, float, float]:
        """(sigma_px, omega_half_frames, tau_rad) of a spec at a level of the given diameter.

        ``auto`` widths follow the quick MC's trial radius r = max(diameter / 3, 0.2 deg) of the
        level: sigma = r x (pixels per radian on detector 0) / 4, omega = max(1, r / frame width)
        frames, tau = r / 4. Without a diameter (no level): sigma 3 px, omega 1 frame, tau
        0.05 deg."""
        if diameter_rad is None:
            r = None
        else:
            r = max(diameter_rad / 3.0, math.radians(0.2))
        sig = spec.sigma_px or (3.0 if r is None else r * self.px_per_rad() / 4.0)
        omh = spec.omega_half or (1.0 if r is None else max(1.0, r / self.width))
        tau = (
            math.radians(spec.tau_deg)
            if spec.tau_deg
            else (math.radians(0.05) if r is None else r / 4.0)
        )
        return float(sig), float(omh), float(tau)

    def evaluate_many(
        self,
        R: Any,
        voxel_vertices: Any,
        phase: int = 0,
        spec: Optional[SmoothSpec] = None,
        diameter_rad: Optional[float] = None,
        chunk: int = 256,
    ) -> np.ndarray:
        """Smooth cost (lower = better, in [0, 1]) of each orientation of ``R`` (B, 3, 3)."""
        assert spec is not None and spec.active, "evaluate_many needs an active SmoothSpec"
        return self.evaluate_specs(R, voxel_vertices, phase, [spec], diameter_rad, chunk)[0]

    def evaluate_specs(
        self,
        R: Any,
        voxel_vertices: Any,
        phase: int,
        specs: Sequence[SmoothSpec],
        diameter_rad: Optional[float] = None,
        chunk: int = 256,
    ) -> List[np.ndarray]:
        """``evaluate_many`` for several specs at once: the geometry and the lit-pixel queries are
        shared, the results equal separate calls."""
        R = np.asarray(R, dtype=np.float64)
        if R.ndim == 2:
            R = R[None]
        outs = [np.ones(len(R)) for _ in specs]
        if phase not in self._g:
            return outs
        widths = [self.level_widths(sp_, diameter_rad) for sp_ in specs]
        ub_all = max(
            [0.5 + 1e-9] + [0.5 + 4.0 * w[0] for sp_, w in zip(specs, widths) if "E1" in sp_.arms]
        )
        with torch.no_grad():
            for s0 in range(0, len(R), chunk):
                pr = self.predict(R[s0 : s0 + chunk], voxel_vertices, phase)
                B = pr.n_batch
                b = pr.b.numpy()
                dcache: Dict[Any, Any] = {"ub": ub_all}
                for out, sp_, (sig, omh, tau) in zip(outs, specs, widths):
                    img_arms = bool(sp_.arms & {"S", "E1", "E2", "E3"})
                    e4 = "E4" in sp_.arms
                    parts = []
                    if img_arms or not e4:
                        W, q = self.quality_terms(pr, sp_, sig, omh, dcache)
                        sw = np.bincount(b, weights=W, minlength=B)
                        sq = np.bincount(b, weights=W * q, minlength=B)
                        parts.append(np.where(sw > 0, 1.0 - sq / np.maximum(sw, 1e-300), 1.0))
                    if e4:
                        W4, rho = self._e4_terms(pr, voxel_vertices, tau, sp_)
                        W4n, rn = W4.numpy(), rho.numpy()
                        sw = np.bincount(b, weights=W4n, minlength=B)
                        sr = np.bincount(b, weights=W4n * rn, minlength=B)
                        parts.append(np.where(sw > 0, sr / np.maximum(sw, 1e-300), 1.0))
                    out[s0 : s0 + B] = np.mean(parts, axis=0)
        return outs

    def e4_cost_torch(
        self,
        R: torch.Tensor,
        voxel_vertices: Any,
        phase: int = 0,
        tau_deg: float = 0.05,
        spec: Optional[SmoothSpec] = None,
    ) -> torch.Tensor:
        """E4 cost of ONE orientation ``R`` (3, 3) float64 torch tensor, differentiable (autograd)
        through the predicted rays and omegas; the assignment to observed rays is held fixed at
        the evaluated point (a piecewise-smooth function with a gradient almost everywhere)."""
        sp_ = spec or SmoothSpec(arms=frozenset({"E4"}))
        pr = self.predict(R[None], voxel_vertices, phase)
        W, rho = self._e4_terms(pr, voxel_vertices, math.radians(tau_deg), sp_)
        sw = W.sum()
        return torch.where(sw > 0, (W * rho).sum() / torch.clamp(sw, min=1e-300), torch.ones(()))


# --------------------------------------------------------------------------------------------
# ranking hook
# --------------------------------------------------------------------------------------------
@dataclass
class SmoothRanker:
    """Ranking key of the coarse levels: smooth cost of every candidate (lower = better)."""

    cost: SmoothCost
    spec: SmoothSpec
    voxel_vertices: Any
    phase: int = 0
    n_evals: int = 0
    history: List[Dict[str, Any]] = field(default_factory=list)

    def applies(self, level: int) -> bool:
        return level < self.spec.levels

    def __call__(self, level: int, candidates: Sequence[Any], diameter_rad: float) -> np.ndarray:
        R = np.stack([np.asarray(c.orientation, dtype=np.float64) for c in candidates])
        key = self.cost.evaluate_many(
            R, self.voxel_vertices, self.phase, self.spec, diameter_rad=diameter_rad
        )
        self.n_evals += len(R)
        self.history.append(dict(level=level, n=len(R)))
        return key


def make_ranker(
    setup: Any,
    spec: SmoothSpec,
    voxel_vertices: Any,
    phase: int,
    eta_limit: float,
    min_sin_eta: float,
) -> SmoothRanker:
    """``SmoothRanker`` of a reconstruction setup (all reflections, the local cost's peak set)."""
    sc = SmoothCost(
        simulator=setup.simulator,
        detector_list=setup.detector_list,
        range_map=setup.range_map,
        exp_data=setup.exp_data,
        sample=setup.sample,
        structure_list=setup.structure_list,
        eta_limit=eta_limit,
        min_sin_eta=min_sin_eta,
    )
    return SmoothRanker(cost=sc, spec=spec, voxel_vertices=voxel_vertices, phase=phase)
