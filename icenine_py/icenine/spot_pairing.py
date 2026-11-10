"""Geometric pairing of diffraction spots across two detectors (no search over orientations).

Two detector distances see one diffracted ray twice. A ray is a straight line, so a spot on
detector 0 and a spot on detector 1 that lie on one line through the illuminated sample plane
(lab z = ``z_plane``, the layer perpendicular to the rotation axis) are one peak, and the line's
intersection with the plane is where the peak came from. This module does that once per frame, in
the lab frame, with batched numpy over all (detector-0, detector-1) component pairs of a frame.

Pure geometry. It takes component centroids (continuous pixel coordinates, pixel centre = index
+ 0.5, the convention of ``ray_overlap``), detector geometry and tolerances, and returns candidate
pairs with their origin on the sample plane, the origin's linearised uncertainty, a cost and an
assignment. It does not read images and does not use the cost function.

Pieces
------
* ``DetectorGeometry`` / ``geometry_from_detectors``: detector plane, pixel basis and size
  (the same fields as ``ray_overlap._det_np``).
* ``pixel_to_lab``: lab point of a continuous pixel coordinate (inverse of the forward map in
  ``ForwardSimulation._simulate_peaks_batched`` stage 5).
* ``line_origin``: intersection of the line through two lab points with the plane z = z_plane,
  and its linearised uncertainty from the pixel tolerances of the four centroid coordinates.
* ``SampleHull``: convex hull of the sample in the sample plane; signed distance of a point.
* ``PairingParams``: tolerance terms. They follow the arm A window of the ray rule (voxel size
  x (L1/L0 - 1), pixel size x L1/L0, arc motion r x d_omega, detector registration) where they
  apply; see the class docstring for which term enters where.
* ``pair_frame``: candidate pairs of one frame, scored and assigned.
* ``attribute_peaks``: peaks (origins in the sample frame) within a radius of voxel centroids,
  with a ``scipy.spatial.cKDTree``.

Frame and beam assumptions
--------------------------
The lab frame is the sample frame rotated about the lab z axis by omega: lab = Rz(omega) sample.
This holds only when the sample-to-lab transform of the configuration is the identity (zero
SampleLocation and SampleOrientation) and the beam runs along +x with no deflection
(``BEAM_DIR``; BeamDeflectionChiLaue 0). Example2.ManyGrains satisfies this, and the Task D loader
(``scripts/eval_reduction/pairing_common.load_geometry``) asserts it. Any other set-up needs the
general sample-to-lab transform here first.

Optional score terms beyond the line/hull test: the scattering angle 2 theta of the line against
the ring angles of the crystal structure, and the pixel-count ratio of the two spots (a planar
region seen along one direction has the same area on both detectors up to the detector tilts).
They are used only when ``ring_two_theta`` / ``shape_weight`` are set.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

BEAM_DIR = np.array([1.0, 0.0, 0.0])
SHAPE_OFFSET_PX = 4.0


@dataclass(frozen=True)
class DetectorGeometry:
    """Detector plane and pixel grid in the lab frame (float64 numpy)."""

    pos: np.ndarray
    org: np.ndarray
    bj: np.ndarray
    bk: np.ndarray
    pw: float
    ph: float
    hw: float
    hh: float
    n: np.ndarray
    d: float

    def distance_x(self) -> float:
        """Detector distance along the beam: the lab x where the beam axis (lab x) meets the
        detector plane n . x + d = 0, i.e. -d / n_x. For Example2 (StandardGeometry.2Det) this is
        the LabFrameLocation x, 3.362 and 5.387 mm. It is not the lab x of the pixel-(0, 0)
        corner (3.263 and 5.290 mm there)."""
        return float(-self.d / self.n[0])


def geometry_from_detectors(detector_list: Sequence[Any]) -> List[DetectorGeometry]:
    """DetectorGeometry of each ``icenine.detector.Detector`` (via the ray rule's cached arrays)."""
    from .ray_overlap import _det_np

    out = []
    for det in detector_list:
        g = _det_np(det)
        out.append(
            DetectorGeometry(
                pos=g["pos"],
                org=g["org"],
                bj=g["bj"],
                bk=g["bk"],
                pw=g["pw"],
                ph=g["ph"],
                hw=g["hw"],
                hh=g["hh"],
                n=g["n"],
                d=g["d"],
            )
        )
    return out


def pixel_to_lab(g: DetectorGeometry, col: np.ndarray, row: np.ndarray) -> np.ndarray:
    """Lab points (..., 3) of continuous pixel coordinates (col, row); centre = index + 0.5."""
    j = np.asarray(col, dtype=np.float64) * g.pw - g.hw
    k = np.asarray(row, dtype=np.float64) * g.ph - g.hh
    return (g.pos + g.org) + j[..., None] * g.bj + k[..., None] * g.bk


def lab_to_pixel(g: DetectorGeometry, pts: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Continuous (col, row) of lab points, projected on the detector basis (stage 5 of the sim)."""
    rel = np.asarray(pts, dtype=np.float64) - (g.pos + g.org)
    return (rel @ g.bj + g.hw) / g.pw, (rel @ g.bk + g.hh) / g.ph


def rot_z(omega: np.ndarray) -> np.ndarray:
    """Rotation matrices (..., 3, 3) about the lab z axis: lab = Rz(omega) @ sample."""
    c, s = np.cos(omega), np.sin(omega)
    z, o = np.zeros_like(c), np.ones_like(c)
    return np.stack([c, -s, z, s, c, z, z, z, o], axis=-1).reshape(np.shape(c) + (3, 3))


def sample_to_lab_xy(q: np.ndarray, omega: Any) -> np.ndarray:
    """Rotate in-plane sample points (..., 2) to the lab frame at ``omega`` (broadcasts)."""
    c, s = np.cos(omega), np.sin(omega)
    return np.stack([c * q[..., 0] - s * q[..., 1], s * q[..., 0] + c * q[..., 1]], axis=-1)


def lab_to_sample_xy(p: np.ndarray, omega: Any) -> np.ndarray:
    """Inverse of :func:`sample_to_lab_xy`."""
    return sample_to_lab_xy(p, -np.asarray(omega))


def _origin_xyz(p0: np.ndarray, p1: np.ndarray, z_plane: float) -> Tuple[np.ndarray, np.ndarray]:
    """Intersection (..., 3) of the line p0 -> p1 with z = z_plane, and the parameter s
    (origin = p0 + s (p1 - p0); s < 0 for a physical ray that starts before detector 0)."""
    u = p1 - p0
    uz = u[..., 2]
    safe = np.where(np.abs(uz) > 1e-12, uz, 1e-12)
    s = (z_plane - p0[..., 2]) / safe
    return p0 + s[..., None] * u, s


def line_origin(
    g0: DetectorGeometry,
    g1: DetectorGeometry,
    c0: np.ndarray,
    r0: np.ndarray,
    c1: np.ndarray,
    r1: np.ndarray,
    tol0_px: float,
    tol1_px: float,
    z_plane: float = 0.0,
) -> Dict[str, np.ndarray]:
    """Origin on the sample plane of the line through two spot centroids, with its uncertainty.

    All inputs broadcast to one shape. Returns lab ``origin`` (..., 3), ray direction ``dirn``
    (unit, ..., 3), ``s`` (see :func:`_origin_xyz`), ``sigma_mm`` (the in-plane length of the
    origin shift caused by moving each of the four centroid coordinates by its tolerance, added in
    quadrature; large when the ray is nearly parallel to the plane) and ``cond`` (|u_z| / |u|,
    the z component of the unit direction: small means ill-conditioned).
    """
    c0, r0, c1, r1 = np.broadcast_arrays(
        *(np.asarray(a, dtype=np.float64) for a in (c0, r0, c1, r1))
    )
    p0 = pixel_to_lab(g0, c0, r0)
    p1 = pixel_to_lab(g1, c1, r1)
    org, s = _origin_xyz(p0, p1, z_plane)
    u = p1 - p0
    nu = np.linalg.norm(u, axis=-1)
    dirn = u / nu[..., None]
    var = np.zeros(c0.shape)
    for which, tol in ((0, tol0_px), (1, tol1_px)):
        for axis in (0, 1):
            dc = np.zeros(4)
            dc[2 * which + axis] = 1.0
            q0 = pixel_to_lab(g0, c0 + dc[0], r0 + dc[1])
            q1 = pixel_to_lab(g1, c1 + dc[2], r1 + dc[3])
            o2, _ = _origin_xyz(q0, q1, z_plane)
            shift = np.linalg.norm((o2 - org)[..., :2], axis=-1) * tol
            var += shift**2
    return {
        "origin": org,
        "dirn": dirn,
        "s": s,
        "sigma_mm": np.sqrt(var),
        "cond": np.abs(dirn[..., 2]),
    }


@dataclass
class SampleHull:
    """Convex hull of the sample in the sample plane (sample frame, mm)."""

    points: np.ndarray
    eq: np.ndarray = field(
        init=False
    )  # (F, 3): unit normal (x, y) and offset; inside: n.q + b <= 0

    def __post_init__(self) -> None:
        self.points = np.asarray(self.points, dtype=np.float64).reshape(-1, 2)
        self.eq = ConvexHull(self.points).equations

    def outside_distance(self, q: np.ndarray) -> np.ndarray:
        """Distance (mm) of points (..., 2) outside the hull; 0 inside (the largest facet excess,
        which is the exact distance for a point nearest a facet and a lower bound near a vertex)."""
        v = q @ self.eq[:, :2].T + self.eq[:, 2]
        return np.clip(v.max(axis=-1), 0.0, None)

    @property
    def radius(self) -> float:
        """Largest distance of a hull point from the rotation axis."""
        return float(np.linalg.norm(self.points, axis=1).max())


@dataclass
class PairingParams:
    """Tolerance and scoring parameters of :func:`pair_frame`.

    Tolerance terms, all in the units named, following the ray rule's arm A window:
      centroid_tol_px        half-width of a centroid coordinate on each detector, in that
                             detector's own pixels (pixel quantisation; the centroid of a planar
                             region seen along one direction is exact up to the pixel grid, so
                             this is one pixel by default). Both centroids are measured, so the
                             detector-0 error is not scaled by L1/L0 onto detector 1: the line fit
                             (``line_origin``) propagates each detector's error on its own.
      registration_mm        extra half-width on detector 1 only: registration (config
                             DetectorSpacingDeviation) plus voxel_size x (L1/L0 - 1) when
                             ``voxel_term`` is set (arm A isotropic part, read as half-widths).
      voxel_size_mm          voxel edge length.
      omega_step             omega step of the range map (radians); the hull test is made at the
                             frame-centre omega, so the origin may move by r x omega_step / 2 along
                             the arc: that enters the hull margin (``arc_term``).
    Hull margin (mm) = voxel_size + arc motion of the point (if ``arc_term``).
    """

    centroid_tol_px: float = 1.0
    registration_mm: float = 0.0
    voxel_size_mm: float = 0.0
    omega_step: float = 0.0
    voxel_term: bool = False
    arc_term: bool = True
    z_plane: float = 0.0
    ring_two_theta: Optional[np.ndarray] = None  # radians, sorted
    ring_tol_scale: float = 1.0  # multiplies the propagated angular tolerance
    shape_weight: float = 0.0  # 0 = off
    size_ratio_gate: float = 1.6  # reject pixel-count ratios beyond this (not on edge spots)
    sigma_gate_mm: float = 1.0e9  # reject pairs whose origin sigma is above this
    max_partners: int = 1  # a spot of unknown multiplicity pairs with at most this many
    multi_npix: int = 10**9  # spots with more pixels than this may have ``multi_partners``
    multi_partners: int = 3


def pair_frame(
    geoms: Sequence[DetectorGeometry],
    hull: SampleHull,
    omega_c: float,
    comp0: Dict[str, np.ndarray],
    comp1: Dict[str, np.ndarray],
    params: PairingParams,
    return_gates: bool = False,
) -> Dict[str, np.ndarray]:
    """Candidate pairs of one frame and their assignment.

    ``comp0`` / ``comp1`` hold per-component arrays ``c_col``, ``c_row`` (mean pixel index, the
    +0.5 pixel-centre shift is added here), ``npix`` and optionally ``edge`` (touches the detector
    border; the size gate is skipped for those) and ``w``, ``h`` (bounding-box width and height in
    pixels, used by the shape term). Returns arrays over the gated candidate pairs:
    ``i0``, ``i1``, ``cost``, ``origin_lab`` (n, 3), ``origin_sample`` (n, 2), ``dirn`` (n, 3,
    unit ray direction), ``sigma_mm``,
    ``cond``, ``dist_out`` (distance outside the hull), ``chosen`` (bool, the assignment) and the
    counts ``n_cand_i0`` (how many candidates each detector-0 component had). With
    ``return_gates`` the (n0, n1) bool matrices ``gate_ray`` (physical ray direction), ``gate_hull``
    (origin within the hull margin, uncertainty included), ``gate_sigma``, ``gate_ring`` and
    ``gate_size`` are added (True = passes; gates not in use are all True).
    """
    g0, g1 = geoms[0], geoms[1]
    n0, n1 = len(comp0["c_col"]), len(comp1["c_col"])
    empty = {k: np.zeros(0) for k in ("i0", "i1", "cost", "sigma_mm", "cond", "dist_out", "chosen")}
    empty["origin_lab"] = np.zeros((0, 3))
    empty["origin_sample"] = np.zeros((0, 2))
    empty["dirn"] = np.zeros((0, 3))
    empty["n_cand_i0"] = np.zeros(n0, dtype=np.int64)
    if n0 == 0 or n1 == 0:
        return empty
    c0 = comp0["c_col"][:, None] + 0.5
    r0 = comp0["c_row"][:, None] + 0.5
    c1 = comp1["c_col"][None, :] + 0.5
    r1 = comp1["c_row"][None, :] + 0.5
    ratio = g1.distance_x() / g0.distance_x()
    # detector 1 is measured, not predicted from detector 0: its own pixel quantisation plus the
    # detector-1 terms (the x L1/L0 factor of the ray rule's arm A would count detector 0 twice)
    tol1_px = (
        params.centroid_tol_px
        + (
            params.registration_mm
            + (params.voxel_size_mm * (ratio - 1.0) if params.voxel_term else 0.0)
        )
        / g1.pw
    )
    lo = line_origin(g0, g1, c0, r0, c1, r1, params.centroid_tol_px, tol1_px, params.z_plane)
    org = lo["origin"]
    q = lab_to_sample_xy(org[..., :2], omega_c)
    rad = np.linalg.norm(org[..., :2], axis=-1)
    margin = params.voxel_size_mm + (rad * params.omega_step / 2.0 if params.arc_term else 0.0)
    dist_out = hull.outside_distance(q)
    # a physical ray starts before detector 0 (s < 0) and moves away from the plane
    g_ray = (lo["s"] < 0) & (lo["dirn"][..., 0] > 0)
    g_hull = dist_out <= margin + lo["sigma_mm"]
    g_sigma = lo["sigma_mm"] <= params.sigma_gate_mm
    g_ring = np.ones((n0, n1), dtype=bool)
    g_size = np.ones((n0, n1), dtype=bool)
    cost = np.zeros((n0, n1))
    if params.ring_two_theta is not None:
        two_theta = np.arccos(np.clip(lo["dirn"] @ BEAM_DIR, -1.0, 1.0))
        rings = np.asarray(params.ring_two_theta)
        d_ring = np.abs(two_theta[..., None] - rings).min(axis=-1)
        # angular tolerance: pixel tolerance of both centroids over the separation
        sep = np.linalg.norm(pixel_to_lab(g1, c1, r1) - pixel_to_lab(g0, c0, r0), axis=-1)
        tol_ang = params.ring_tol_scale * (params.centroid_tol_px * g0.pw + tol1_px * g1.pw) / sep
        g_ring = d_ring <= tol_ang
        cost += d_ring / tol_ang
    if params.shape_weight > 0:
        # the same planar region seen along one direction: equal area and bounding box up to the
        # detector tilts; the worst log-ratio of pixel count and of width and height
        # (each plus 4 pixels)
        # (an offset of SHAPE_OFFSET_PX keeps one-pixel quantisation of thin spots from counting)
        off = SHAPE_OFFSET_PX
        lr = np.abs(np.log((comp1["npix"][None, :] + off) / (comp0["npix"][:, None] + off)))
        for key in ("w", "h"):
            if key in comp0 and key in comp1:
                lr = np.maximum(
                    lr,
                    np.abs(np.log((comp1[key][None, :] + off) / (comp0[key][:, None] + off))),
                )
        edge = np.zeros((n0, n1), dtype=bool)
        if "edge" in comp0:
            edge |= comp0["edge"][:, None]
        if "edge" in comp1:
            edge |= comp1["edge"][None, :]
        lr = np.where(edge, 0.0, lr)
        g_size = lr <= np.log(params.size_ratio_gate)
        cost += params.shape_weight * lr / np.log(params.size_ratio_gate)
    # position preference among equals: the origin's distance outside the hull, in margins
    cost += 0.05 * dist_out / np.maximum(margin, 1e-9)
    ok = g_ray & g_hull & g_sigma & g_ring & g_size
    i0, i1 = np.nonzero(ok)
    out = {
        "i0": i0,
        "i1": i1,
        "cost": cost[i0, i1],
        "origin_lab": org[i0, i1],
        "origin_sample": q[i0, i1],
        "dirn": lo["dirn"][i0, i1],
        "sigma_mm": lo["sigma_mm"][i0, i1],
        "cond": lo["cond"][i0, i1],
        "dist_out": dist_out[i0, i1],
        "n_cand_i0": np.bincount(i0, minlength=n0),
    }
    out["chosen"] = assign_pairs(i0, i1, out["cost"], comp0["npix"], comp1["npix"], params)
    if return_gates:
        out.update(
            gate_ray=g_ray, gate_hull=g_hull, gate_sigma=g_sigma, gate_ring=g_ring, gate_size=g_size
        )
    return out


def assign_pairs(
    i0: np.ndarray,
    i1: np.ndarray,
    cost: np.ndarray,
    npix0: np.ndarray,
    npix1: np.ndarray,
    params: PairingParams,
) -> np.ndarray:
    """Greedy assignment by ascending cost: each spot takes at most ``max_partners`` partners
    (``multi_partners`` when it has more than ``multi_npix`` pixels, i.e. may be several sources).
    Ties are broken by input order (stable sort). Returns a bool mask over the candidates."""
    chosen = np.zeros(len(i0), dtype=bool)
    cap0 = np.where(npix0 > params.multi_npix, params.multi_partners, params.max_partners)
    cap1 = np.where(npix1 > params.multi_npix, params.multi_partners, params.max_partners)
    used0 = np.zeros(len(npix0), dtype=np.int64)
    used1 = np.zeros(len(npix1), dtype=np.int64)
    for k in np.argsort(cost, kind="stable"):
        a, b = int(i0[k]), int(i1[k])
        if used0[a] < cap0[a] and used1[b] < cap1[b]:
            chosen[k] = True
            used0[a] += 1
            used1[b] += 1
    return chosen


def attribute_peaks(
    origins: np.ndarray, voxel_xy: np.ndarray, radius: np.ndarray
) -> List[np.ndarray]:
    """For each voxel, indices of the peaks whose sample-frame origin lies within ``radius`` mm of
    its centroid. ``radius`` is a scalar or one value per peak (spot footprint added by the caller).

    Uses a cKDTree of the origins; with per-peak radii it queries the maximum radius and filters.
    """
    origins = np.asarray(origins, dtype=np.float64).reshape(-1, 2)
    voxel_xy = np.asarray(voxel_xy, dtype=np.float64).reshape(-1, 2)
    rad = np.broadcast_to(np.asarray(radius, dtype=np.float64), (len(origins),))
    tree = cKDTree(origins)
    cand = tree.query_ball_point(voxel_xy, float(rad.max()) if len(rad) else 0.0)
    out = []
    for v, idx in enumerate(cand):
        idx = np.asarray(idx, dtype=np.int64)
        if len(idx):
            d = np.linalg.norm(origins[idx] - voxel_xy[v], axis=1)
            idx = idx[d <= rad[idx]]
        out.append(idx)
    return out
