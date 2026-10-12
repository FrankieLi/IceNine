"""
Ray-consistency overlap rule (opt-in; default off). Used only by the COARSE filtering.

TWO COSTS DIFFER BY DESIGN (owner decision 2026-10-11)
------------------------------------------------------
The rule is applied only in the coarse levels 0 .. ``levels-1`` (default 0-2) of the adaptive
search: the pixel-radius screen, the re-cost of the screened candidates, the quick MC of those
levels and the ranking cut that decides which candidates go to the next level. FindOptimal, the
VarianceMinimizing stage, the final cost and the confidence use the strict C++ predicted-spot rule
(``CppOverlapRule``), so a rule-passing orientation is still judged by the C++ cost at the end.
The reconstructor level loop (``AdaptiveVoxelReconstructor.reconstruct_voxel``) selects the cost
function per level; ``VoxelCostFunction(ray_rule=...)`` carries the rule.

The rule
--------
Per peak, after the PREDICTED spot is bright on detector 0, cast the OBSERVED ray from the voxel
centre (at the peak's predicted omega) through the bright detector-0 pixel(s) the predicted spot
overlaps, to every later detector d (n_det = 2 in the studied geometry; for n_det > 2 the same
straight line is intersected with each later detector, its tolerance scaled by that detector's
L_d/L_0, not re-aimed through bright pixels). Cases on detector d (C++ lit = predicted spot in
bounds, bright = overlap):

  a. predicted spot bright                                  -> valid (nothing changes)
  b. predicted in bounds and dark, observed ray bright      -> valid (counted as overlap)
  c. predicted and observed both in bounds and dark         -> invalid (the C++ rule)
  d. predicted and observed both out of bounds              -> valid (the C++ rule tolerates unlit)
  e. predicted in bounds and dark, observed out of bounds   -> valid (treated as unlit)
  f. predicted out of bounds, observed in bounds and dark   -> valid (the C++ rule forgives the
     unlit detector). Owner decision 2026-10-12: the variant f1 (invalid, a lit non-overlap) was
     dropped; at the truth it wrongly rejected about 3% of real peaks (0.966 kept vs 0.994,
     step 1) and steps 1-2 showed no benefit. May revisit.
  (predicted out of bounds, observed ray bright: left as the C++ rule, counted as ``f_bright``.)

"Observed ray bright" means a bright detector-d pixel inside the tolerance window around where
the ray meets the detector (in bounds: any pixel of the window is on the detector).

Tolerance window on detector d (half-widths in detector-d pixels)
  arm R: an isotropic square of ``radius_px`` pixels around the pixel containing the ray hit (the
         same convention as the pixel-radius screen).
  arm A: an axis-aligned anisotropic window per peak. Half-width (mm, isotropic part) =
         voxel_size*(L_d/L_0 - 1) + pixel_0*(L_d/L_0) + registration, plus the arc motion
         r*d_omega*(L_d/L_0 - 1) along the arc direction (z x P, projected on the detector axes
         j and k). r is the voxel distance from the rotation axis, d_omega the omega step of the
         range map, registration the DetectorSpacingDeviation. Each term is taken as a half-width
         (a conservative reading). Blurring from orientation gradients is ignored.

Approximations: the detector-0 "bright pixels the spot overlaps" are the bright pixels of the
+-pixel_radius window around the projected centre (screen, radius > 0) or of the bounding box of
the projected triangle (local cost, radius 0); at most ``MAX_RAYS`` of them are used.
"""

from dataclasses import dataclass
from typing import Any, List, Optional, Tuple

import numpy as np

MAX_RAYS = 64


@dataclass(frozen=True)
class RayRule:
    """Parameters of the ray-consistency rule (see module docstring).

    Attributes:
        arm: "R" isotropic radius window or "A" anisotropic window.
        radius_px: Window half-width in pixels for arm R.
        registration_mm: Registration term of arm A (config DetectorSpacingDeviation).
        levels: The rule is used for search levels 0 .. levels-1 (default 3 = levels 0-2).
    """

    arm: str = "R"
    radius_px: int = 2
    registration_mm: float = 0.0
    levels: int = 3

    def __post_init__(self) -> None:
        if self.arm not in ("R", "A"):
            raise ValueError(f"arm must be 'R' or 'A', got {self.arm!r}")
        if self.radius_px < 0:
            raise ValueError("radius_px must be >= 0")

    def label(self) -> str:
        t = f"R{self.radius_px}" if self.arm == "R" else "A"
        return t


CASE_NAMES = ("a", "b", "c", "d", "e", "f_dark", "f_bright", "det0_bright")


def _det_np(detector: Any) -> dict:
    """Detector geometry as float64 numpy (cached on the detector object)."""
    g = getattr(detector, "_ray_geom", None)
    if g is None:
        plane = detector._detector_plane
        g = dict(
            n=plane.normal.double().numpy(),
            d=float(plane.d),
            pos=detector._position.double().numpy(),
            org=detector._lab_frame_coord_origin.double().numpy(),
            bj=detector._lab_frame_basis_j.double().numpy(),
            bk=detector._lab_frame_basis_k.double().numpy(),
            pw=float(detector.pixel_width),
            ph=float(detector.pixel_height),
            hw=float(detector.pixel_half_width),
            hh=float(detector.pixel_half_height),
        )
        detector._ray_geom = g
    return g


def _box_flags(img: np.ndarray, c_lo: int, c_hi: int, r_lo: int, r_hi: int) -> Tuple[bool, bool]:
    """(any pixel of the box on the detector, any bright pixel) for the inclusive box."""
    rows, cols = img.shape
    if c_hi < 0 or r_hi < 0 or c_lo >= cols or r_lo >= rows:
        return False, False
    sub = img[max(r_lo, 0) : min(r_hi, rows - 1) + 1, max(c_lo, 0) : min(c_hi, cols - 1) + 1]
    return True, bool(sub.any())


def ray_flags(
    rule: RayRule,
    detector_list: List[Any],
    exp_data: Any,
    wedge: np.ndarray,
    lab_centres: np.ndarray,
    det_hit0: np.ndarray,
    det0_px: np.ndarray,
    pixel_radius: int,
    voxel_size: float,
    omega_step: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Observed-ray flags for the peaks, shape (M, n_det) uint8 each: obs_in, obs_bright.

    Args:
        wedge: (M,) image index of each peak.
        lab_centres: (M, 3) voxel centre in the lab frame at the peak's omega.
        det_hit0: (M,) bool, the predicted ray reaches detector 0's plane for all 3 vertices.
        det0_px: (M, 3, 2) projected triangle vertices on detector 0 as [col, row] (float).
        pixel_radius: radius of the detector-0 window (> 0), or 0 for the triangle bounding box.
        voxel_size: voxel edge length (mm), used by arm A.
        omega_step: omega step of the range map (radians), used by arm A.

    Only peaks whose detector-0 window holds a bright pixel get non-zero flags; the C stage
    applies them only where the predicted detector-0 spot overlaps.
    """
    m = len(wedge)
    n_det = len(detector_list)
    obs_in = np.zeros((m, n_det), dtype=np.uint8)
    obs_bright = np.zeros((m, n_det), dtype=np.uint8)
    if n_det < 2:
        return obs_in, obs_bright
    g0 = _det_np(detector_list[0])
    geoms = [_det_np(d) for d in detector_list]
    for p in np.nonzero(det_hit0)[0]:
        w = int(wedge[p])
        img0 = np.asarray(exp_data.get_image(w, 0).get_binary_numpy())
        rows, cols = img0.shape
        if pixel_radius > 0:
            cx, cy = int(det0_px[p, 0, 0]), int(det0_px[p, 0, 1])
            c_lo, c_hi = cx - pixel_radius, cx + pixel_radius
            r_lo, r_hi = cy - pixel_radius, cy + pixel_radius
        else:
            c_lo = int(np.floor(det0_px[p, :, 0].min()))
            c_hi = int(np.floor(det0_px[p, :, 0].max()))
            r_lo = int(np.floor(det0_px[p, :, 1].min()))
            r_hi = int(np.floor(det0_px[p, :, 1].max()))
        c_lo, c_hi = max(c_lo, 0), min(c_hi, cols - 1)
        r_lo, r_hi = max(r_lo, 0), min(r_hi, rows - 1)
        if c_lo > c_hi or r_lo > r_hi:
            continue
        rr, cc = np.nonzero(img0[r_lo : r_hi + 1, c_lo : c_hi + 1])
        if len(rr) == 0:
            continue
        if len(rr) > MAX_RAYS:
            rr, cc = rr[:MAX_RAYS], cc[:MAX_RAYS]
        cc = cc + c_lo
        rr = rr + r_lo
        # detector-0 pixel centres in the lab frame (inverse of the forward map in Stage C)
        j = (cc + 0.5) * g0["pw"] - g0["hw"]
        k = (rr + 0.5) * g0["ph"] - g0["hh"]
        q0 = g0["pos"] + g0["org"] + j[:, None] * g0["bj"] + k[:, None] * g0["bk"]
        origin = lab_centres[p].astype(np.float64)
        dirs = q0 - origin
        for d in range(1, n_det):
            g = geoms[d]
            denom = dirs @ g["n"]
            ok = np.abs(denom) > 1e-12
            if not ok.any():
                continue
            t = np.where(ok, -(origin @ g["n"] + g["d"]) / np.where(ok, denom, 1.0), 0.0)
            ok &= t > 0
            if not ok.any():
                continue
            hit = origin + t[:, None] * dirs
            rel = hit - g["pos"] - g["org"]
            col = ((rel @ g["bj"]) + g["hw"]) / g["pw"]
            row = ((rel @ g["bk"]) + g["hh"]) / g["ph"]
            if rule.arm == "R":
                lo_c = np.trunc(col).astype(np.int64) - rule.radius_px
                hi_c = lo_c + 2 * rule.radius_px
                lo_r = np.trunc(row).astype(np.int64) - rule.radius_px
                hi_r = lo_r + 2 * rule.radius_px
            else:
                ratio = t  # L_d / L_0 (the ray reaches detector 0 at t = 1)
                iso = (
                    voxel_size * (ratio - 1.0)
                    + 0.5 * (g0["pw"] + g0["ph"]) * ratio
                    + rule.registration_mm
                )
                rad = float(np.hypot(origin[0], origin[1]))
                arc = np.array([-origin[1], origin[0], 0.0])
                nrm = np.linalg.norm(arc)
                arc = arc / nrm if nrm > 0 else arc
                arc_mm = rad * omega_step * (ratio - 1.0)
                hc = (iso + arc_mm * abs(arc @ g["bj"])) / g["pw"]
                hr = (iso + arc_mm * abs(arc @ g["bk"])) / g["ph"]
                lo_c = np.floor(col - hc).astype(np.int64)
                hi_c = np.floor(col + hc).astype(np.int64)
                lo_r = np.floor(row - hr).astype(np.int64)
                hi_r = np.floor(row + hr).astype(np.int64)
            img = np.asarray(exp_data.get_image(w, d).get_binary_numpy())
            any_in = False
            any_bright = False
            for i in np.nonzero(ok)[0]:
                inb, br = _box_flags(img, int(lo_c[i]), int(hi_c[i]), int(lo_r[i]), int(hi_r[i]))
                any_in |= inb
                if br:
                    any_bright = True
                    break
            obs_in[p, d] = any_in
            obs_bright[p, d] = any_bright
    return obs_in, obs_bright


def _geom_array(detector_list: List[Any]) -> np.ndarray:
    """(n_det, 20) float64 geometry for the C ray check: n d pos org bj bk pw ph hw hh."""
    rows = []
    for det in detector_list:
        g = _det_np(det)
        rows.append(
            np.concatenate(
                [g["n"], [g["d"]], g["pos"], g["org"], g["bj"], g["bk"]]
                + [[g["pw"], g["ph"], g["hw"], g["hh"]]]
            )
        )
    return np.ascontiguousarray(np.array(rows, dtype=np.float64))


def ray_flags_c(
    rule: RayRule,
    detector_list: List[Any],
    image_list: List[np.ndarray],
    local_wedge: np.ndarray,
    lab_centres: np.ndarray,
    det_hit0: np.ndarray,
    det0_px: np.ndarray,
    pixel_radius: int,
    voxel_size: float,
    omega_step: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """C version of :func:`ray_flags` (same flags; ``ray_flags`` stays the reference).

    ``image_list`` is indexed ``[local_wedge * n_det + d]`` as for ``stage_d_overlap``;
    ``local_wedge`` is the (M,) int32 index into the unique wedges.
    """
    from ._rasterize import ray_flags as _c_ray_flags

    return _c_ray_flags(
        image_list,
        np.ascontiguousarray(local_wedge, dtype=np.int32),
        np.ascontiguousarray(det_hit0, dtype=np.uint8),
        np.ascontiguousarray(det0_px, dtype=np.float64),
        np.ascontiguousarray(lab_centres, dtype=np.float64),
        _geom_array(detector_list),
        len(detector_list),
        int(pixel_radius),
        0 if rule.arm == "R" else 1,
        int(rule.radius_px),
        float(voxel_size),
        float(omega_step),
        float(rule.registration_mm),
    )


def apply_cases_python(
    detector_lit: List[bool],
    spot_overlap: List[bool],
    obs_in_p: np.ndarray,
    obs_bright_p: np.ndarray,
    cases: Optional[np.ndarray],
) -> Tuple[int, int]:
    """Python twin of the C stage's ray block for one peak (modifies the lists in place).

    Returns (delta_pixel_overlap, delta_pixel_on_detector). ``cases`` (length 8) is incremented
    when given. Does nothing unless the predicted detector-0 spot overlaps.
    """
    d_ov = d_on = 0
    if not spot_overlap[0]:
        return 0, 0
    if cases is not None:
        cases[7] += 1
    for d in range(1, len(detector_lit)):
        oi, ob = bool(obs_in_p[d]), bool(obs_bright_p[d])
        if spot_overlap[d]:
            if cases is not None:
                cases[0] += 1
            continue
        if detector_lit[d]:
            if ob:
                spot_overlap[d] = True
                d_ov += 1
                if cases is not None:
                    cases[1] += 1
            elif not oi:
                detector_lit[d] = False
                d_on -= 1
                if cases is not None:
                    cases[4] += 1
            elif cases is not None:
                cases[2] += 1
        elif not oi:
            if cases is not None:
                cases[3] += 1
        elif ob:
            if cases is not None:
                cases[6] += 1
        else:
            if cases is not None:
                cases[5] += 1
    return d_ov, d_on


def ray_rule_from_config(config: Any) -> Optional[RayRule]:
    """RayRule from the config keys RayOverlapRule (off | R<n> | A), or None.

    The registration term of arm A is the config's DetectorSpacingDeviation (mm).
    """
    spec = str(getattr(config, "ray_overlap_rule", "off"))
    if spec == "off":
        return None
    reg = float(getattr(config, "detector_dist_deviation", 0.0))
    if spec == "A":
        return RayRule(arm="A", registration_mm=reg)
    return RayRule(arm="R", radius_px=int(spec[1:]), registration_mm=reg)
