"""
Stage 0 evaluation tools for the toy orientation NN.

Provides
- a batched, float64 re-implementation of the simulator's per-peak ray tracing
  (`BatchedObserver`) that returns, for many orientation offsets at once, which
  ROI peaks are observed, in which frame, and where their spot vertices land;
- an exact Bayes baseline (`ExactBayes`): the posterior over the orientation
  offset delta given the thresholded data, computed by importance sampling
  restricted to the set of offsets that reproduce the observed data exactly;
- per-axis error metrics (rotation about the stage axis vs. perpendicular).

Conventions follow docs/nn_inverse_problem_formulation.md: the orientation is
R(delta) = exp([delta]x) R_nom with delta a rotation vector in the sample frame;
angles in this module's public API are degrees, internals use radians.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import torch
from scipy.spatial.transform import Rotation

from .constants import KEV_OVER_HBAR_C_IN_ANG
from .diffraction_core import get_scattering_omegas_torch
from .image_data import ImageData

DEG = np.pi / 180.0
ABSENT = -2  # signature value for a peak that is not observed


# ---------------------------------------------------------------------------
# Rotation helpers
# ---------------------------------------------------------------------------


def _norm(x: torch.Tensor, keepdim: bool = False) -> torch.Tensor:
    """Euclidean norm over the last dimension. Much faster than torch.linalg.norm
    for length-3 vectors on large batches."""
    return torch.sqrt((x * x).sum(dim=-1, keepdim=keepdim))


def rotvec_to_matrix(delta_rad: torch.Tensor) -> torch.Tensor:
    """Batched Rodrigues formula. (B, 3) rotation vectors -> (B, 3, 3) matrices."""
    theta = _norm(delta_rad, keepdim=True)  # (B, 1)
    safe = torch.where(theta > 1e-12, theta, torch.ones_like(theta))
    k = delta_rad / safe
    zero = torch.zeros_like(k[..., 0])
    K = torch.stack(
        [zero, -k[..., 2], k[..., 1], k[..., 2], zero, -k[..., 0], -k[..., 1], k[..., 0], zero],
        dim=-1,
    ).reshape(*k.shape[:-1], 3, 3)
    s = torch.sin(theta)[..., None]
    c = torch.cos(theta)[..., None]
    eye = torch.eye(3, dtype=delta_rad.dtype, device=delta_rad.device)
    return eye + s * K + (1.0 - c) * (K @ K)


def quaternions_to_offsets_deg(q_wxyz: np.ndarray, R_nom: np.ndarray) -> np.ndarray:
    """Rotation-vector offsets (degrees) of orientations given as [w,x,y,z]
    quaternions relative to a nominal matrix, delta = rotvec(R R_nom^T)."""
    q = np.asarray(q_wxyz, dtype=np.float64)
    R = Rotation.from_quat(np.concatenate([q[:, 1:], q[:, :1]], axis=1)).as_matrix()
    return Rotation.from_matrix(R @ R_nom.T).as_rotvec() / DEG


def offsets_to_matrices(delta_deg: np.ndarray, R_nom: np.ndarray) -> np.ndarray:
    """(N, 3) offsets in degrees -> (N, 3, 3) orientation matrices exp([delta]x) R_nom."""
    return Rotation.from_rotvec(np.asarray(delta_deg, dtype=np.float64) * DEG).as_matrix() @ R_nom


def offsets_to_quaternions(delta_deg: np.ndarray, R_nom: np.ndarray) -> np.ndarray:
    """(N, 3) offsets in degrees -> (N, 4) [w,x,y,z] quaternions with w >= 0."""
    xyzw = Rotation.from_matrix(offsets_to_matrices(delta_deg, R_nom)).as_quat()
    q = np.concatenate([xyzw[:, 3:], xyzw[:, :3]], axis=1)
    return np.where(q[:, :1] < 0, -q, q)


def sample_prior_offsets(n: int, radius_deg: float, rng: np.random.Generator) -> np.ndarray:
    """Training prior: uniform in the ball |delta| <= radius_deg. (n, 3) degrees."""
    directions = rng.normal(size=(n, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    radii = radius_deg * rng.uniform(size=(n, 1)) ** (1.0 / 3.0)
    return directions * radii


def sample_fixed_magnitude_offsets(
    n: int, magnitude_deg: float, rng: np.random.Generator
) -> np.ndarray:
    """Offsets of a fixed magnitude in uniformly random directions. (n, 3) degrees."""
    directions = rng.normal(size=(n, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    return directions * magnitude_deg


# ---------------------------------------------------------------------------
# Batched observer
# ---------------------------------------------------------------------------


@dataclass
class Observation:
    """Result of observing B candidate offsets on the M ROI peaks."""

    present: torch.Tensor  # (B, M) bool: peak observed on its home detector
    frame: torch.Tensor  # (B, M) long: frame index, -1 where the peak is not in range
    omega: torch.Tensor  # (B, M) float: Bragg crossing (rad), valid where Bragg-observable
    verts: torch.Tensor  # (B, M, 3, 2) float: spot vertices (col, row) on the home detector


class BatchedObserver:
    """Vectorised float64 version of the simulator's per-peak ray tracing.

    For a fixed voxel, nominal orientation and fixed ROI peak set, `observe`
    returns everything the thresholded data depend on: for every ROI peak,
    whether it is observed, its frame index and its three spot vertices in pixel
    coordinates. The maths mirrors ForwardSimulation._simulate_peaks (serial
    path): Bragg solve, Rz(omega) @ base rotation, beam reflection, eta filter,
    ray-plane intersection with the peak's home detector.
    """

    def __init__(
        self, R_nom, voxel_vertices, sample, detector_list, range_map, exp_setup, roi_list
    ):
        self.dtype = torch.float64
        d = self.dtype
        self.roi_list = roi_list
        self.M = len(roi_list)
        self.R_nom = torch.as_tensor(np.asarray(R_nom, dtype=np.float64), dtype=d)
        self.g_hkl = torch.stack([p.g_hkl.to(d) for p in roi_list])  # (M, 3)
        self.branch = torch.tensor([p.omega_branch for p in roi_list])  # (M,)
        self.det_idx = torch.tensor([p.detector_index for p in roi_list])  # (M,)
        self.vertices = torch.as_tensor(np.asarray(voxel_vertices), dtype=d)  # (3, 3) sample frame

        bm = sample.sample_to_lab_matrix.to(d)
        self.base = bm[:3, :3].clone()
        self.translation = bm[:3, 3].clone()
        self.beam_dir = torch.as_tensor(np.asarray(exp_setup.get_beam_direction()), dtype=d)
        self.energy = float(exp_setup.beam_energy)
        self.chi_laue = float(exp_setup.get_beam_deflection_chi_laue())
        self.eta_limit = float(exp_setup.get_eta_limit())

        self.range_low = float(range_map.low)
        self.range_width = float(range_map.width)
        self.range_n = int(range_map.num_intervals)
        self.range_index = torch.tensor(
            [-1 if i is None else int(i) for i in range_map.index_list], dtype=torch.long
        )
        self.frame_width_rad = abs(self.range_width)
        # The Jacobian uses |range_width| while nominal_offsets / measurement_features use the
        # signed width and assume the frame index of a bin equals its position; both hold only if:
        assert self.range_width > 0, f"range_width must be positive, got {self.range_width}"
        _valid = self.range_index[self.range_index >= 0]
        assert len(_valid) == 0 or bool(
            (_valid[1:] - _valid[:-1] == 1).all()
        ), "range_index must be contiguous and increasing over the valid omega bins"

        # Per-peak home-detector geometry, gathered once.
        geo = []
        for det in detector_list:
            plane = det.detector_plane
            geo.append(
                dict(
                    normal=plane.normal.to(d),
                    plane_d=float(plane.d),
                    offset=(det._position + det._lab_frame_coord_origin).to(d),
                    basis_j=det._lab_frame_basis_j.to(d),
                    basis_k=det._lab_frame_basis_k.to(d),
                    half_w=float(det.pixel_half_width),
                    half_h=float(det.pixel_half_height),
                    pix_w=float(det.pixel_width),
                    pix_h=float(det.pixel_height),
                )
            )
        # All detectors' planes, for the simulator's rule that a vertex missing any
        # detector plane drops the peak on every detector.
        self.all_normals = torch.stack([g["normal"] for g in geo])  # (D, 3)
        self.all_plane_d = torch.tensor([g["plane_d"] for g in geo], dtype=d)  # (D,)
        di = self.det_idx.tolist()
        self.d_normal = torch.stack([geo[i]["normal"] for i in di])  # (M, 3)
        self.d_plane = torch.tensor([geo[i]["plane_d"] for i in di], dtype=d)  # (M,)
        self.d_offset = torch.stack([geo[i]["offset"] for i in di])  # (M, 3)
        self.d_bj = torch.stack([geo[i]["basis_j"] for i in di])  # (M, 3)
        self.d_bk = torch.stack([geo[i]["basis_k"] for i in di])  # (M, 3)
        self.d_hw = torch.tensor([geo[i]["half_w"] for i in di], dtype=d)  # (M,)
        self.d_hh = torch.tensor([geo[i]["half_h"] for i in di], dtype=d)
        self.d_pw = torch.tensor([geo[i]["pix_w"] for i in di], dtype=d)
        self.d_ph = torch.tensor([geo[i]["pix_h"] for i in di], dtype=d)
        self.d_ncols = torch.tensor([detector_list[i].num_cols for i in di])  # (M,)
        self.d_nrows = torch.tensor([detector_list[i].num_rows for i in di])

    # -- helpers ---------------------------------------------------------

    def _bragg(self, delta_rad: torch.Tensor):
        """Bragg solve for all peaks. Returns g_s (B,M,3), omega (B,M), observable (B,M)."""
        B = delta_rad.shape[0]
        Rm = rotvec_to_matrix(delta_rad) @ self.R_nom  # (B, 3, 3)
        g_s = torch.einsum("bij,mj->bmi", Rm, self.g_hkl)  # (B, M, 3)
        flat = g_s.reshape(-1, 3)
        res = get_scattering_omegas_torch(flat, _norm(flat), self.energy, self.chi_laue)
        omega1 = res.omega1.reshape(B, self.M)
        omega2 = res.omega2.reshape(B, self.M)
        omega = torch.where(self.branch[None, :] == 1, omega1, omega2)
        return g_s, omega, res.observable.reshape(B, self.M)

    def _frame_index(self, omega: torch.Tensor) -> torch.Tensor:
        f = (omega - self.range_low) / self.range_width
        n = torch.floor(f).long()
        valid = (f >= 0) & (n < self.range_n)
        idx = self.range_index[torch.clamp(n, 0, self.range_n - 1)]
        return torch.where(valid, idx, torch.full_like(idx, -1))

    @staticmethod
    def _rotz(v: torch.Tensor, c: torch.Tensor, s: torch.Tensor) -> torch.Tensor:
        """Rotate vectors (..., 3) about z by an angle with cos c, sin s (broadcastable)."""
        x = c * v[..., 0] - s * v[..., 1]
        y = s * v[..., 0] + c * v[..., 1]
        return torch.stack([x, y, v[..., 2].expand_as(x)], dim=-1)

    # -- main entry points ------------------------------------------------

    def observe_frames(self, delta_deg: torch.Tensor):
        """Cheap first stage: frame index and eta acceptance only (no ray tracing).

        Returns (frame (B,M) long, ok (B,M) bool) where ok = Bragg-observable, in
        a measured frame, and passing the eta filter. A peak that is `ok` may
        still be lost later if a spot vertex misses its detector.
        """
        delta_rad = delta_deg.to(self.dtype) * DEG
        g_s, omega, observable = self._bragg(delta_rad)
        frame = self._frame_index(omega)
        c, s = torch.cos(omega), torch.sin(omega)
        sd = g_s / _norm(g_s, keepdim=True)
        ln = self._rotz(sd @ self.base.T, c, s)  # (B, M, 3) lab-frame scattering direction
        dot = (ln * self.beam_dir).sum(-1, keepdim=True)
        rd = self.beam_dir - 2.0 * dot * ln
        rd_n = rd / _norm(rd, keepdim=True)
        eta = torch.atan2(rd_n[..., 1].abs(), rd_n[..., 2].abs())
        ok = observable & (frame >= 0) & (eta < self.eta_limit)
        return frame, ok, (omega, c, s, rd)

    def observe(self, delta_deg: torch.Tensor) -> Observation:
        """Full observation of the ROI peaks at B candidate offsets (degrees)."""
        frame, ok, (omega, c, s, rd) = self.observe_frames(delta_deg)
        B = frame.shape[0]
        # Spot vertices in the lab frame: (B, M, 3 verts, 3)
        base_v = self.vertices @ self.base.T  # (3, 3)
        lab_v = self._rotz(base_v[None, None, :, :], c[..., None], s[..., None]) + self.translation
        denom = (self.d_normal[None] * rd).sum(-1)  # (B, M)
        denom_ok = denom.abs() > 1e-8
        safe = torch.where(denom_ok, denom, torch.ones_like(denom))
        numer = -(
            (self.d_normal[None, :, None, :] * lab_v).sum(-1) + self.d_plane[None, :, None]
        )  # (B,M,3)
        t = numer / safe[..., None]
        hit = denom_ok[..., None] & (t > 0)
        inter = lab_v + t[..., None] * rd[:, :, None, :]
        loc = inter - self.d_offset[None, :, None, :]
        j = (loc * self.d_bj[None, :, None, :]).sum(-1)
        k = (loc * self.d_bk[None, :, None, :]).sum(-1)
        col = (j + self.d_hw[None, :, None]) / self.d_pw[None, :, None]
        row = (k + self.d_hh[None, :, None]) / self.d_ph[None, :, None]
        verts = torch.stack([col, row], dim=-1)
        # A spot is only recorded if it overlaps the pixel grid (same rule as
        # orientation_nn.spot_overlaps_grid: bounding box of truncated vertices).
        trunc = torch.where(verts < 0, torch.full_like(verts, -1.0), torch.floor(verts))
        cols, rows = trunc[..., 0], trunc[..., 1]
        ncols = self.d_ncols[None, :].to(cols.dtype)
        nrows = self.d_nrows[None, :].to(cols.dtype)
        on_grid = (
            (cols.amax(-1) >= 0)
            & (cols.amin(-1) <= ncols - 1)
            & (rows.amax(-1) >= 0)
            & (rows.amin(-1) <= nrows - 1)
        )
        # Every vertex must hit every detector plane (ForwardSimulation._simulate_peaks).
        all_planes = torch.ones_like(ok)
        for dn, dd in zip(self.all_normals, self.all_plane_d):
            den = (rd * dn).sum(-1)  # (B, M)
            den_ok = den.abs() > 1e-8
            num = -((lab_v * dn).sum(-1) + dd)  # (B, M, 3)
            t_d = num / torch.where(den_ok, den, torch.ones_like(den))[..., None]
            all_planes = all_planes & den_ok & (t_d > 0).all(dim=-1)
        present = ok & hit.all(dim=-1) & all_planes & on_grid
        return Observation(present=present, frame=frame, omega=omega, verts=verts)

    def sin_eta(self, delta_deg: torch.Tensor) -> torch.Tensor:
        """|sin eta| of every ROI peak at the given offsets, (B, M). eta is the azimuth of
        the diffracted beam about the incident beam, measured from the rotation axis
        (small |sin eta| = near-axis peaks, which drift many frames per degree)."""
        _frame, _ok, (_omega, _c, _s, rd) = self.observe_frames(delta_deg)
        rd_n = rd / _norm(rd, keepdim=True)
        return torch.sin(torch.atan2(rd_n[..., 1].abs(), rd_n[..., 2].abs())).abs()

    def peak_context(self, h_deg: float = 0.01) -> torch.Tensor:
        """Per-peak context features at the nominal orientation, (M, 14 + ndet) float32.

        For each ROI spot: its spot-motion Jacobian Gamma (col, row vs. rotation about sample
        x, y, z; px per degree, /20), the frame gradient d omega*/d delta (dimensionless),
        |sin eta|, sin theta, a one-hot of its detector, its nominal centroid (col, row) as a
        fraction of the detector size, and its nominal frame scaled to [-1, 1]. This is
        everything a shared per-peak encoder needs to know about how the peak responds to an
        orientation offset, so it can be applied to peaks of any orientation.
        Central differences of the observer at +-h_deg; spots missing at a probe get zeros.
        """
        pts = [np.zeros(3)]
        for i in range(3):
            e = np.zeros(3)
            e[i] = h_deg
            pts += [e, -e]
        obs = self.observe(torch.as_tensor(np.array(pts), dtype=self.dtype))
        cent = obs.verts.mean(dim=2)  # (7, M, 2)
        om = obs.omega  # (7, M)
        ok = obs.present.all(dim=0)  # (M,)
        gamma = torch.stack(
            [(cent[1 + 2 * i] - cent[2 + 2 * i]) / (2 * h_deg) for i in range(3)], dim=-1
        )
        gomega = torch.stack(
            [(om[1 + 2 * i] - om[2 + 2 * i]) / (2 * h_deg * DEG) for i in range(3)], -1
        )
        gamma = torch.where(ok[:, None, None], gamma, torch.zeros_like(gamma))
        gomega = torch.where(ok[:, None], gomega, torch.zeros_like(gomega))
        sin_eta = self.sin_eta(torch.zeros(1, 3, dtype=self.dtype))[0]
        wavenumber = KEV_OVER_HBAR_C_IN_ANG * self.energy
        sin_theta = _norm(self.g_hkl) / (2.0 * wavenumber)
        ndet = len(self.all_normals)
        onehot = torch.nn.functional.one_hot(self.det_idx, ndet).to(self.dtype)
        centroid = cent[0] / torch.stack([self.d_ncols, self.d_nrows], dim=-1).to(self.dtype)
        frame = 2.0 * obs.frame[0].to(self.dtype) / self.range_n - 1.0
        ctx = torch.cat(
            [
                (gamma.reshape(-1, 6) / 20.0),
                gomega,
                sin_eta[:, None],
                sin_theta[:, None],
                onehot,
                centroid,
                frame[:, None],
            ],
            dim=-1,
        )
        return ctx.float()

    @staticmethod
    def vertex_keys(obs: Observation) -> torch.Tensor:
        """(B, M, 6) long: spot vertices (col0,row0,col1,row1,col2,row2) truncated the
        way the C++ rasteriser does (negative -> -1, else int())."""
        v = obs.verts.reshape(*obs.verts.shape[:2], 6)
        return torch.where(v < 0, torch.full_like(v, -1.0), torch.floor(v)).long()


def nominal_offsets(observer: "BatchedObserver") -> np.ndarray:
    """Exact nominal sub-window offsets per ROI peak, (M, 3) float64.

    Column 0/1: the nominal spot centroid's fractional position inside its pixel
    (col - floor(col), row - floor(row)), i.e. the offset between the window centre
    (WindowSpec.from_nominal cuts at floor(centroid) - W/2) and the exact nominal centroid.
    Column 2: the nominal crossing omega relative to the centre of its nominal frame, in
    frames: (omega_nom - omega_centre(frame0)) / range_width, in [-0.5, 0.5).
    measurement_features gives (lit centroid - window centre, frame - frame0); subtracting
    these columns turns them into (measurement - exact nominal prediction).
    """
    out = observer.observe(torch.zeros(1, 3, dtype=observer.dtype))
    cen = out.verts[0].mean(dim=1).numpy()
    frac = cen - np.floor(cen)
    bin_of = {int(w): b for b, w in enumerate(observer.range_index.tolist()) if w >= 0}
    om = out.omega[0].numpy()
    fr = out.frame[0].numpy()
    ff = np.zeros(observer.M)
    for m in range(observer.M):
        b = bin_of.get(int(fr[m]))
        if b is not None:
            centre = observer.range_low + (b + 0.5) * observer.range_width
            ff[m] = (om[m] - centre) / observer.range_width
    return np.concatenate([frac, ff[:, None]], axis=1)


def pair_index(roi_list: Sequence[Any]) -> np.ndarray:
    """(M,) int: index of the other-detector entry of the same diffracted ray, or -1.

    Two ROI entries are the same ray when they share (reflection_index, omega_branch) and
    differ in detector_index (define_roi_set with detectors="all" makes one entry per
    detector the spot overlaps). With more than two detectors the partner is the entry on
    the next detector index (cyclic) that is present.
    """
    groups = {}
    for i, p in enumerate(roi_list):
        groups.setdefault((p.reflection_index, p.omega_branch), []).append(i)
    out = np.full(len(roi_list), -1, dtype=np.int64)
    for members in groups.values():
        if len(members) < 2:
            continue
        members = sorted(members, key=lambda i: roi_list[i].detector_index)
        for a, i in enumerate(members):
            out[i] = members[(a + 1) % len(members)]
    return out


def lit_pixel_set(key: Tuple[int, ...], num_cols: int, num_rows: int) -> frozenset:
    """Pixels the simulator's rasteriser lights for a spot with truncated vertex
    key (col0,row0,col1,row1,col2,row2): the same truncate -> Sutherland-Hodgman
    clip -> round -> scanline fill pipeline as ImageData.add_triangle_scanline."""
    polygon = [
        (float(key[0]), float(key[1])),
        (float(key[2]), float(key[3])),
        (float(key[4]), float(key[5])),
    ]
    clipped = ImageData._sutherland_hodgman_clip(
        polygon, 0.0, float(num_cols - 1), 0.0, float(num_rows - 1)
    )
    if len(clipped) < 3:
        return frozenset()
    pixels = ImageData._scanline_fill([(round(x), round(y)) for x, y in clipped])
    return frozenset((c, r) for c, r in pixels if 0 <= c < num_cols and 0 <= r < num_rows)


@dataclass
class WindowSpec:
    """Fixed per-spot window placement for a dataset.

    Each spot's window is window_size x window_size pixels with integer top-left
    corner (col0, row0) on its detector, centred on the spot's nominal centroid, and
    covers frames frame0 - K .. frame0 + K around its nominal frame frame0.
    """

    window_size: int
    frame_half_width: int  # K
    col0: np.ndarray  # (M,) int
    row0: np.ndarray  # (M,) int
    frame0: np.ndarray  # (M,) int

    @classmethod
    def from_nominal(
        cls, observer: "BatchedObserver", window_size: int, frame_half_width: int
    ) -> "WindowSpec":
        out = observer.observe(torch.zeros(1, 3, dtype=observer.dtype))
        if not bool(out.present.all()):
            raise ValueError("every ROI spot must be present at the nominal orientation")
        centroid = out.verts[0].mean(dim=1).numpy()  # (M, 2) col, row
        half = window_size // 2
        return cls(
            window_size=window_size,
            frame_half_width=frame_half_width,
            col0=np.floor(centroid[:, 0]).astype(np.int64) - half,
            row0=np.floor(centroid[:, 1]).astype(np.int64) - half,
            frame0=out.frame[0].numpy().astype(np.int64),
        )


def render_windows(
    observer: "BatchedObserver", spec: WindowSpec, deltas_deg: np.ndarray, chunk: int = 256
) -> Tuple[torch.Tensor, torch.Tensor]:
    """Frame-coded thresholded windows, rendered from the observer's exact spot data.

    For each offset and ROI spot, the window holds the pixels the simulator's
    rasteriser lights for that spot (lit_pixel_set, i.e. exactly the thresholded
    detector data restricted to the window), each set to 1 + (frame - frame0 + K);
    unlit pixels are 0. With alpha = 0 a spot lights a single frame, so this is a
    lossless encoding of the (frame, pixels) data inside the window. A spot that is
    absent, or whose frame lies outside frame0 +- K, gives an all-zero window.

    Returns (windows uint8 (N, M, W, W), status uint8 (N, M)) where status is
    0 = inside the window, 1 = absent, 2 = frame outside +-K, 3 = lit pixels
    partly or entirely outside the window.
    """
    N, M, W, K = len(deltas_deg), observer.M, spec.window_size, spec.frame_half_width
    windows = torch.zeros(N, M, W, W, dtype=torch.uint8)
    status = torch.zeros(N, M, dtype=torch.uint8)
    ncols = observer.d_ncols.tolist()
    nrows = observer.d_nrows.tolist()
    cache: Dict[Tuple[int, ...], frozenset] = {}
    for a in range(0, N, chunk):
        d = torch.as_tensor(np.asarray(deltas_deg[a : a + chunk], dtype=np.float64))
        obs = observer.observe(d)
        keys = BatchedObserver.vertex_keys(obs).tolist()
        present = obs.present.tolist()
        frames = obs.frame.tolist()
        for i in range(len(d)):
            n = a + i
            for m in range(M):
                if not present[i][m]:
                    status[n, m] = 1
                    continue
                offset = frames[i][m] - int(spec.frame0[m])
                if abs(offset) > K:
                    status[n, m] = 2
                    continue
                key = tuple(keys[i][m])
                ck = (m,) + key
                pix = cache.get(ck)
                if pix is None:
                    pix = lit_pixel_set(key, ncols[m], nrows[m])
                    cache[ck] = pix
                code = 1 + offset + K
                c0, r0 = int(spec.col0[m]), int(spec.row0[m])
                outside = False
                for c, r in pix:
                    x, y = c - c0, r - r0
                    if 0 <= x < W and 0 <= y < W:
                        windows[n, m, y, x] = code
                    else:
                        outside = True
                if outside:
                    status[n, m] = 3
        if len(cache) > 500_000:
            cache.clear()
    return windows, status


def render_distractor_windows(
    observer: "BatchedObserver",
    spec: WindowSpec,
    sources: List["BatchedObserver"],
    source_deltas_deg: List[np.ndarray],
    chunk: int = 128,
    source_active: Optional[List[np.ndarray]] = None,
) -> torch.Tensor:
    """Spots of *other* scatterers (neighbour voxels, twins) that land in the target's windows.

    sources: one BatchedObserver per distractor source (its own orientation, voxel vertices and
    ROI peak set); source_deltas_deg[s] (N, 3) is that source's orientation offset for each of
    the N samples. Every present source spot whose lit pixels overlap a target window of the
    same detector and whose frame lies within the window's +-K frames is drawn into that window
    with the same frame coding as render_windows. Returns uint8 (N, M, W, W), zero where no
    distractor pixel falls. source_active[s] (N,) bool switches source s on/off per sample
    (default: always on). It is a separate layer: combine_windows overlays it on the target's
    own windows, where the target's pixels always win, so the target's spots are never altered.
    """
    N, M, W, K = len(source_deltas_deg[0]), observer.M, spec.window_size, spec.frame_half_width
    out = torch.zeros(N, M, W, W, dtype=torch.uint8)
    tdet = observer.det_idx.numpy()
    col0, row0, frame0 = spec.col0, spec.row0, spec.frame0
    cache: Dict[Tuple[int, ...], frozenset] = {}
    for si, (src, dsrc) in enumerate(zip(sources, source_deltas_deg)):
        active = None if source_active is None else np.asarray(source_active[si], dtype=bool)
        sdet = src.det_idx.numpy()
        ncols, nrows = src.d_ncols.tolist(), src.d_nrows.tolist()
        same_det = sdet[:, None] == tdet[None, :]  # (Q, M)
        for a in range(0, N, chunk):
            o = src.observe(torch.as_tensor(np.asarray(dsrc[a : a + chunk], dtype=np.float64)))
            keys = BatchedObserver.vertex_keys(o).numpy()  # (B, Q, 6)
            present = o.present.numpy()
            if active is not None:
                present = present & active[a : a + chunk, None]
            frames = o.frame.numpy()  # (B, Q)
            cs, rs = keys[..., 0::2], keys[..., 1::2]
            # bounding box of the lit pixels (rounding can add one pixel: margin 1)
            c_lo, c_hi = cs.min(-1) - 1, cs.max(-1) + 1
            r_lo, r_hi = rs.min(-1) - 1, rs.max(-1) + 1
            hit = (
                present[:, :, None]
                & same_det[None]
                & (c_hi[:, :, None] >= col0[None, None])
                & (c_lo[:, :, None] < (col0 + W)[None, None])
                & (r_hi[:, :, None] >= row0[None, None])
                & (r_lo[:, :, None] < (row0 + W)[None, None])
                & (np.abs(frames[:, :, None] - frame0[None, None]) <= K)
            )
            for i, q, m in zip(*np.nonzero(hit)):
                key = tuple(int(k) for k in keys[i, q])
                ck = (int(q),) + key
                pix = cache.get(ck)
                if pix is None:
                    pix = lit_pixel_set(key, ncols[q], nrows[q])
                    cache[ck] = pix
                code = 1 + int(frames[i, q]) - int(frame0[m]) + K
                n = a + i
                for c, r in pix:
                    x, y = c - int(col0[m]), r - int(row0[m])
                    if 0 <= x < W and 0 <= y < W and out[n, m, y, x] == 0:
                        out[n, m, y, x] = code
            if len(cache) > 500_000:
                cache.clear()
    return out


def combine_windows(windows: torch.Tensor, distractors: Optional[torch.Tensor]) -> torch.Tensor:
    """Overlay a distractor layer under the target's windows: target pixels always win."""
    if distractors is None:
        return windows
    return torch.where(windows > 0, windows, distractors)


@dataclass
class CorruptionConfig:
    """Random corruptions applied to frame-coded windows (per entry, independently).

    p_miss: the whole spot is not recorded (window zeroed); p_flip: each lit pixel is dropped
    and each 4-neighbour of a lit pixel is lit with this probability (threshold jitter at the
    spot edge); p_hot: a window gets one isolated hot pixel with this probability; p_blob: a
    window gets a spurious 2-4 px blob (one random frame) with this probability. neighbours:
    overlay the dataset's distractor layer (neighbour voxels / twin spots) before the rest.
    """

    neighbours: bool = True
    p_miss: float = 0.1
    p_flip: float = 0.05
    p_hot: float = 0.05
    p_blob: float = 0.1

    @classmethod
    def named(cls, name: str) -> Optional["CorruptionConfig"]:
        if name in ("none", "clean"):
            return None
        if name == "neighbours":
            return cls(True, 0.0, 0.0, 0.0, 0.0)
        if name == "noise":
            return cls(False)
        if name == "all":
            return cls(True)
        raise ValueError(f"unknown corruption {name!r}")


def corrupt_windows(
    windows: torch.Tensor,
    distractors: Optional[torch.Tensor],
    cfg: Optional[CorruptionConfig],
    frame_half_width: int,
    gen: Optional[torch.Generator] = None,
    valid: Optional[torch.Tensor] = None,
) -> torch.Tensor:
    """Apply cfg to uint8 windows (..., W, W) (leading dims free). Random draws come from gen.

    valid: optional bool mask over the leading dims (e.g. entry index < n_peaks). Entries marked
    invalid (zero padding) are returned all-zero, so corruption cannot create "present" spots in
    them. The random draws are unchanged. Default None: every entry, padding included, is
    corrupted (the behaviour of all reported runs).
    """
    if cfg is None:
        return windows
    x = combine_windows(windows, distractors if cfg.neighbours else None).clone()
    shape = x.shape
    lead = shape[:-2]
    H, W = shape[-2:]

    def rnd(*size: int) -> torch.Tensor:
        return torch.rand(*size, generator=gen, device="cpu").to(x.device)

    if cfg.p_flip > 0:
        lit = x > 0
        padded = torch.nn.functional.pad(x.float(), (1, 1, 1, 1))
        # brightest (max) neighbour code = a neighbouring lit pixel's frame code
        nb = torch.stack(
            [
                padded[..., 0:-2, 1:-1],
                padded[..., 2:, 1:-1],
                padded[..., 1:-1, 0:-2],
                padded[..., 1:-1, 2:],
            ]
        ).amax(0)
        grow = (~lit) & (nb > 0) & (rnd(*shape) < cfg.p_flip)
        drop = lit & (rnd(*shape) < cfg.p_flip)
        x = torch.where(drop, torch.zeros_like(x), x)
        x = torch.where(grow, nb.to(x.dtype), x)
    n_codes = 2 * frame_half_width + 1
    if cfg.p_hot > 0:
        has = rnd(*lead) < cfg.p_hot
        pos = (rnd(*lead, 2) * torch.tensor([H, W], device=x.device)).long()
        code = (rnd(*lead) * n_codes).long().clamp(max=n_codes - 1) + 1
        yy = torch.arange(H, device=x.device).view(*([1] * len(lead)), H, 1)
        xx = torch.arange(W, device=x.device).view(*([1] * len(lead)), 1, W)
        hot = (
            (yy == pos[..., 0, None, None]) & (xx == pos[..., 1, None, None]) & has[..., None, None]
        )
        x = torch.where(hot & (x == 0), code[..., None, None].to(x.dtype), x)
    if cfg.p_blob > 0:
        has = rnd(*lead) < cfg.p_blob
        size = (rnd(*lead, 2) * 3).long().clamp(max=2) + 2  # 2..4
        pos = (rnd(*lead, 2) * torch.tensor([H - 4, W - 4], device=x.device)).long()
        code = (rnd(*lead) * n_codes).long().clamp(max=n_codes - 1) + 1
        yy = torch.arange(H, device=x.device).view(*([1] * len(lead)), H, 1)
        xx = torch.arange(W, device=x.device).view(*([1] * len(lead)), 1, W)
        blob = (
            (yy >= pos[..., 0, None, None])
            & (yy < (pos[..., 0] + size[..., 0])[..., None, None])
            & (xx >= pos[..., 1, None, None])
            & (xx < (pos[..., 1] + size[..., 1])[..., None, None])
            & has[..., None, None]
        )
        x = torch.where(blob & (x == 0), code[..., None, None].to(x.dtype), x)
    if cfg.p_miss > 0:
        miss = rnd(*lead) < cfg.p_miss
        x = torch.where(miss[..., None, None], torch.zeros_like(x), x)
    if valid is not None:
        x = torch.where(valid.to(x.device)[..., None, None], x, torch.zeros_like(x))
    return x


def corrupt_dataset(
    windows: torch.Tensor,
    distractors: Optional[torch.Tensor],
    name: str,
    frame_half_width: int,
    seed: int = 12345,
    chunk: int = 100,
    valid: Optional[torch.Tensor] = None,
) -> torch.Tensor:
    """Deterministically corrupted copy of a whole window array (the fixed test sets).

    valid: optional bool mask over the leading dims of windows (see corrupt_windows); default
    None corrupts padded entries too (as in the reported runs).

    name: CorruptionConfig.named ("none", "neighbours", "noise", "all"); the same seed and
    chunking always give the same corrupted windows, so the Gauss-Newton baseline and the
    networks are evaluated on identical inputs.
    """
    cfg = CorruptionConfig.named(name)
    if cfg is None:
        return windows
    gen = torch.Generator().manual_seed(seed)
    out = torch.empty_like(windows)
    for a in range(0, len(windows), chunk):
        d = None if distractors is None else distractors[a : a + chunk]
        v = None if valid is None else valid[a : a + chunk]
        out[a : a + chunk] = corrupt_windows(
            windows[a : a + chunk], d, cfg, frame_half_width, gen, valid=v
        )
    return out


def decode_windows(windows: torch.Tensor, frame_half_width: int) -> torch.Tensor:
    """Frame-coded uint8 windows (..., W, W) -> float channels (..., 2, W, W):
    channel 0 = lit (0/1), channel 1 = frame offset / K in [-1, 1] at lit pixels, 0 elsewhere."""
    w = windows.float()
    lit = (w > 0).float()
    frame = torch.where(w > 0, (w - 1.0 - frame_half_width) / frame_half_width, torch.zeros_like(w))
    return torch.stack([lit, frame], dim=-3)


# ---------------------------------------------------------------------------
# Exact Bayes baseline
# ---------------------------------------------------------------------------


class ExactBayes:
    """Posterior over the offset given thresholded data, by importance sampling.

    With no noise, the posterior is the prior restricted to the consistent set
    C(D) = {delta : data(delta) = D}. "Data" is, for each ROI peak, whether it is
    observed, its frame, and (if use_pixels) the exact set of lit pixels, taken
    from the simulator's own rasteriser. Overlaps between different peaks' pixels
    are ignored (they are rare: about one spot per few million pixels).

    The consistent set is tiny compared with the prior, so it is found
    adaptively: start from a Gaussian around delta_true (always a member) scaled
    by the linearised information matrix, shrink until members are found, then
    re-centre and re-scale on the members and estimate the mean and covariance
    with self-normalised importance weights 1/q(delta).
    """

    def __init__(
        self,
        observer: BatchedObserver,
        prior_radius_deg: float,
        use_pixels: bool = True,
        chunk: int = 2048,
    ):
        self.obs = observer
        self.prior_radius = float(prior_radius_deg)
        self.use_pixels = use_pixels
        self.chunk = chunk
        self._equal_cache: Dict[Tuple[int, ...], bool] = {}
        self._fill_cache: Dict[Tuple[int, ...], frozenset] = {}

    # -- pixel sets --------------------------------------------------------

    def _lit_pixels(self, m: int, key: Tuple[int, ...]) -> frozenset:
        ck = (int(self.obs.d_ncols[m]), int(self.obs.d_nrows[m])) + key
        out = self._fill_cache.get(ck)
        if out is None:
            out = lit_pixel_set(key, ck[0], ck[1])
            self._fill_cache[ck] = out
        return out

    def _same_pixels(self, m: int, key: Tuple[int, ...], truth_key: Tuple[int, ...]) -> bool:
        if key == truth_key:
            return True
        ck = (m,) + key + truth_key
        hit = self._equal_cache.get(ck)
        if hit is None:
            hit = self._lit_pixels(m, key) == self._lit_pixels(m, truth_key)
            self._equal_cache[ck] = hit
        return hit

    # -- linearised information (used only to scale the first proposal) ----

    def information_matrix(self, delta_true_deg: np.ndarray, h_deg: float = 0.005) -> np.ndarray:
        """Independent-quantisation information matrix J (per rad^2) from central
        finite differences at delta_true (docs 3.6.4): frames at variance
        dw^2/12, pixel coordinates at 1/12."""
        o = self.obs
        pts = [delta_true_deg]
        for i in range(3):
            e = np.zeros(3)
            e[i] = h_deg
            pts += [delta_true_deg + e, delta_true_deg - e]
        obs = o.observe(torch.as_tensor(np.array(pts), dtype=o.dtype))
        h = h_deg * DEG
        present = obs.present.all(dim=0)  # peaks present at every probe point
        om = obs.omega[:, present].numpy()  # (7, Mp)
        cent = obs.verts[:, present].mean(dim=2).numpy()  # (7, Mp, 2)
        grad_w = np.stack(
            [(om[1 + 2 * i] - om[2 + 2 * i]) / (2 * h) for i in range(3)], axis=-1
        )  # (Mp, 3)
        grad_u = np.stack(
            [(cent[1 + 2 * i] - cent[2 + 2 * i]) / (2 * h) for i in range(3)], axis=-1
        )  # (Mp,2,3)
        J = (12.0 / o.frame_width_rad**2) * grad_w.T @ grad_w
        if self.use_pixels:
            J = J + 12.0 * np.einsum("pai,paj->ij", grad_u, grad_u)
        return J

    # -- membership --------------------------------------------------------

    def _members(self, deltas_deg: torch.Tensor, truth) -> torch.Tensor:
        """Boolean (B,) mask of candidates that reproduce the observed data."""
        o = self.obs
        frame_true, present_true, keys_true = truth
        inside = _norm(deltas_deg) <= self.prior_radius
        member = torch.zeros(len(deltas_deg), dtype=torch.bool)
        need = present_true[None, :]
        for a in range(0, len(deltas_deg), self.chunk):
            sl = slice(a, a + self.chunk)
            frame, ok, _ = o.observe_frames(deltas_deg[sl])
            # Stage 1 (cheap): every truly present peak must be ok and in the same frame.
            stage1 = (((frame == frame_true[None, :]) & ok) | ~need).all(dim=1) & inside[sl]
            idx = torch.where(stage1)[0]
            if len(idx) == 0:
                continue
            obs = o.observe(deltas_deg[sl][idx])
            same = (obs.present == present_true[None, :]).all(dim=1)
            if self.use_pixels:
                keys = o.vertex_keys(obs)  # (S, M, 6)
                differs = (
                    (keys != keys_true[None]).any(dim=-1) & present_true[None, :] & obs.present
                )
                for s_i in torch.where(same & differs.any(dim=1))[0].tolist():
                    for m in torch.where(differs[s_i])[0].tolist():
                        if not self._same_pixels(
                            m, tuple(keys[s_i, m].tolist()), tuple(keys_true[m].tolist())
                        ):
                            same[s_i] = False
                            break
            sub = member[sl]
            sub[idx] = same
            member[sl] = sub
        return member

    def posterior(
        self,
        delta_true_deg: np.ndarray,
        rng: np.random.Generator,
        n_per_round: int = 20000,
        min_ess: float = 300.0,
        max_rounds: int = 12,
        init_scale: float = 1.0 / 256.0,
    ) -> Dict[str, object]:
        """Posterior mean/covariance (degrees, degrees^2) given data at delta_true.

        init_scale multiplies the linearised independent-quantisation covariance to
        give the first proposal; the noise-free cell is far smaller than that
        covariance (about 1/P in variance for P peaks), so starting small saves rounds.
        """
        o = self.obs
        dt = torch.as_tensor(np.asarray(delta_true_deg, dtype=np.float64)[None], dtype=o.dtype)
        obs_true = o.observe(dt)
        truth = (obs_true.frame[0], obs_true.present[0], o.vertex_keys(obs_true)[0])
        n_present = int(truth[1].sum())

        self._equal_cache.clear()  # keyed on the truth of the current case
        self._fill_cache.clear()
        try:
            cov_lin = np.linalg.inv(self.information_matrix(np.asarray(delta_true_deg))) / DEG**2
            np.linalg.cholesky(cov_lin)  # must be positive definite to draw from
        except np.linalg.LinAlgError:
            cov_lin = np.eye(3) * 0.01**2
        mean = np.asarray(delta_true_deg, dtype=np.float64)
        cov = cov_lin.copy()
        scale = init_scale
        log = []

        def draw(mean, cov, n):
            L = np.linalg.cholesky(cov + 1e-30 * np.eye(3))
            z = rng.normal(size=(n, 3))
            x = mean + z @ L.T
            logdet = 2.0 * np.log(np.diag(L)).sum()
            logq = -0.5 * ((z**2).sum(axis=1) + logdet + 3 * np.log(2 * np.pi))
            return x, logq

        result = None
        for rnd in range(max_rounds):
            x, logq = draw(mean, cov * scale, n_per_round)
            m = self._members(torch.as_tensor(x, dtype=o.dtype), truth).numpy()
            n_acc = int(m.sum())
            log.append((rnd, scale, n_acc))
            if n_acc < 30:
                scale /= 4.0
                continue
            w = np.exp(-(logq[m] - logq[m].min()))  # shifted for numerical stability
            w /= w.sum()
            xm = x[m]
            new_mean = (w[:, None] * xm).sum(axis=0)
            diff = xm - new_mean
            new_cov = (w[:, None, None] * (diff[:, :, None] * diff[:, None, :])).sum(axis=0)
            ess = 1.0 / (w**2).sum()
            result = dict(
                mean=new_mean,
                cov=new_cov,
                ess=float(ess),
                n_accept=n_acc,
                rounds=log,
                n_present=n_present,
            )
            if ess >= min_ess and rnd >= 1:
                break
            mean, cov, scale = new_mean, new_cov * 4.0 + 1e-30 * np.eye(3), 1.0
        if result is None:
            # No members were found: report failure explicitly instead of leaking delta_true.
            nan3 = np.full(3, np.nan)
            result = dict(
                mean=nan3,
                cov=np.full((3, 3), np.nan),
                ess=0.0,
                n_accept=0,
                rounds=log,
                n_present=n_present,
            )
        return result


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------


def error_summary(delta_hat_deg: np.ndarray, delta_true_deg: np.ndarray) -> Dict[str, float]:
    """Per-axis error summary. z = rotation about the stage axis, perp = x and y.

    Assumes the base sample rotation maps sample z to lab z (true for Example2, where
    it is the identity), so z is the stage axis. Errors are differences of rotation
    vectors, which is accurate for small offsets. Reports RMS of each, the RMS total, the median misorientation angle between the
    estimated and true orientations, and the fraction of cases whose misorientation
    is below 0.5 deg (the bench_hp_sweep.py success criterion) and below 0.1 deg.
    """
    err = np.asarray(delta_hat_deg, float) - np.asarray(delta_true_deg, float)
    ang = (
        np.linalg.norm(
            Rotation.from_matrix(
                Rotation.from_rotvec(np.asarray(delta_hat_deg) * DEG).as_matrix()
                @ Rotation.from_rotvec(np.asarray(delta_true_deg) * DEG)
                .as_matrix()
                .transpose(0, 2, 1)
            ).as_rotvec(),
            axis=1,
        )
        / DEG
    )
    return dict(
        n=len(err),
        rms_z=float(np.sqrt(np.mean(err[:, 2] ** 2))),
        rms_perp=float(np.sqrt(np.mean((err[:, 0] ** 2 + err[:, 1] ** 2) / 2.0))),
        rms_total=float(np.sqrt(np.mean((err**2).sum(axis=1)))),
        median_angle=float(np.median(ang)),
        success_0p5=float(np.mean(ang < 0.5)),
        success_0p1=float(np.mean(ang < 0.1)),
    )
