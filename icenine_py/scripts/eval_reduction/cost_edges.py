"""Task B: the edges of the cost function, attributed to events (frame crossings, pixel crossings,
detector edge, eta cut, rule flips, denominator changes).

This file holds the measurement core (no library change): a per-peak trace of one cost evaluation,
the classification of the change of every peak between two orientations, and the attribution of the
change of the cost to those events. ``cost_edges_run.py`` drives it over cases and rays,
``cost_edges_summary.py`` turns the results into tables and figures.

Per-peak trace. ``trace_eval`` calls ``VoxelCostFunction.evaluate`` once with three module functions
of
``icenine.cost_functions`` wrapped for the duration of the call (the Bragg solution, the batched
overlap, the C overlap stage), and then calls the C overlap stage once per peak so that every peak
has
its own overlap counts and quality. The library is not changed. A peak is a reflection with one of
its
two omega branches: key = 2 * reflection + branch, stable across orientations.

Cost. cost = 1 - quality, quality = mean over the peaks with at least one pixel on a detector of
(pixel_overlap / pixel_on_detector) * (detectors_overlap / n_detectors). The per-peak values are
recombined here and checked against the library's own aggregate in every trace.

Events between two orientations (A, B), per peak, first match wins:
  presence:   the peak is evaluated in one and not in the other. Reason from the side where it is
  not
              evaluated: bragg_limit (the reflection is not observable: no Bragg solution), eta_cut
              (the eta window removes it), no_frame (its omega has no frame: outside the scanned
              range)
  frame:      the frame index (the image the spot is compared with) changed
  detector:   the pattern of vertices inside the detector area, or the plane-hit flag, changed
  pixel:      a vertex (or the screen centre) moved to another pixel
  subpixel:   same pixels, but the counts the rasteriser returns changed
  rule:       only the ray-rule case counts changed (ray rule only)
"""

import contextlib
import io
import math
import os
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(ICE / "scripts" / "common"))

ABSENT = {0: "bragg_limit", 1: "eta_cut", 2: "no_frame"}
PRESENT = 3
EVENT_TYPES = [
    "bragg_limit",
    "eta_cut",
    "no_frame",
    "frame",
    "detector",
    "pixel",
    "subpixel",
    "rule",
]
# how the plan's event table classifies each type
KIND = {
    "frame": "physics, hard cut",
    "pixel": "artefact",
    "subpixel": "artefact",
    "detector": "geometry, hard cut",
    "eta_cut": "geometry, hard cut",
    "bragg_limit": "geometry, hard cut",
    "no_frame": "geometry, hard cut",
    "rule": "rule design",
    "denominator": "artefact of the metric",
    "unattributed": "unknown",
}


# --------------------------------------------------------------------------------------------
# cost functions
# --------------------------------------------------------------------------------------------
def make_cf(base: SimpleNamespace, setting: str, rule: str) -> Any:
    """setting: 'full' (pixel radius 0, all reflections: the complete cost) or 'screen' (pixel
    radius 3, Q-max 5, centroid: the level-0 screen). rule: 'old' (lenient, library default),
    'cpp' (C++ overlap rule), 'cppE' (C++ rule + CppEdgeRules), 'r2' (ray rule, 2-pixel window)."""
    from icenine.cost_functions import VoxelCostFunction
    from icenine.ray_overlap import RayRule

    kw: Dict[str, Any] = dict(
        simulator=base.simulator,
        detector_list=base.detector_list,
        range_map=base.range_map,
        exp_data=base.exp_data,
        sample=base.sample,
        structure_list=base.structure_list,
        mode="hard",
        eta_limit=base.eta_limit,
        min_sin_eta=0.0,
        center_screen=True,
    )
    if setting == "screen":
        kw.update(pixel_radius=3, max_q=5.0)
    else:
        kw.update(pixel_radius=0)
    if rule == "cpp":
        kw["cpp_overlap_rule"] = True
    elif rule == "cppE":
        kw.update(cpp_overlap_rule=True, cpp_edge_rules=True)
    elif rule == "r2":
        kw["ray_rule"] = RayRule(arm="R", radius_px=2, registration_mm=base.reg_mm)
    elif rule != "old":
        raise ValueError(rule)
    return VoxelCostFunction(**kw)


# --------------------------------------------------------------------------------------------
# per-peak trace
# --------------------------------------------------------------------------------------------
def trace_eval(cf: Any, R: np.ndarray, verts: Any, phase: int) -> Dict[str, Any]:
    """Per-peak record of one evaluation (arrays over the 2 * n_reflections peak keys)."""
    from icenine import cost_functions as CF

    g_hkl = cf._phase_recip_vecs[phase][0]
    nrefl = int(g_hkl.shape[0])
    P = 2 * nrefl
    n_det = len(cf.detector_list)
    radius = int(cf.pixel_radius)
    edge = bool(getattr(cf, "cpp_edge_rules", False))
    cap: Dict[str, Any] = {}
    o_bragg, o_batched, o_stage = (
        CF.get_scattering_omegas_torch,
        CF.calculate_diffraction_overlap_batched,
        CF._c_stage_d_overlap,
    )

    def bragg_w(*a: Any, **k: Any) -> Any:
        r = o_bragg(*a, **k)
        cap["bragg"] = r
        return r

    def batched_w(**kw: Any) -> Any:
        cap["kw"] = kw
        return o_batched(**kw)

    def stage_w(*a: Any) -> Any:
        cap["stage"] = a
        return o_stage(*a)

    CF.get_scattering_omegas_torch, CF.calculate_diffraction_overlap_batched = bragg_w, batched_w
    CF._c_stage_d_overlap = stage_w
    try:
        info = cf.evaluate(np.asarray(R, dtype=np.float32), verts, phase)
    finally:
        CF.get_scattering_omegas_torch, CF.calculate_diffraction_overlap_batched = (
            o_bragg,
            o_batched,
        )
        CF._c_stage_d_overlap = o_stage

    K_sig = 1 if radius > 0 else 3
    out: Dict[str, Any] = dict(
        status=np.zeros(P, np.int8),
        omega=np.full(P, np.nan, np.float64),
        frame=np.full(P, -1, np.int32),
        hit=np.zeros((P, n_det), np.uint8),
        sig=np.zeros((P, n_det, K_sig, 2), np.int32),
        cnt=np.zeros((P, 5), np.int32),  # pixov, pixon, pkov, pkon, ndet
        q=np.zeros(P, np.float64),
        rays=np.zeros((P, 8), np.int32),
        agg=dict(
            pixov=int(info.pixel_overlap),
            pixon=int(info.pixel_on_detector),
            pkov=int(info.peak_overlap),
            pkon=int(info.peak_on_detector),
            quality=float(info.quality),
            nq=int(info.n_quality_points),
        ),
    )
    br = cap.get("bragg")
    if br is None:
        return out
    obs = br.observable.numpy().astype(bool)
    k_obs = np.nonzero(obs)[0]
    K = len(k_obs)
    if K == 0:
        return out
    rows_omega = np.concatenate([br.omega1.numpy()[obs], br.omega2.numpy()[obs]])
    row_key = np.concatenate([2 * k_obs, 2 * k_obs + 1])  # peak key of each of the 2K rows
    out["status"][row_key] = 1  # observable; stays 1 if the eta filter removes it
    out["omega"][row_key] = rows_omega
    kw = cap.get("kw")
    if kw is None:
        return out
    omegas = kw["peak_omegas"].numpy()
    M0 = len(omegas)
    keys = np.zeros(M0, np.int64)
    j = 0
    for p in range(M0):  # the eta filter keeps a subsequence of the rows: match in order
        while j < len(rows_omega) and rows_omega[j] != omegas[p]:
            j += 1
        keys[p] = row_key[j]
        j += 1
    rng = kw["range_map"]
    bins_f = (omegas - rng.low) / rng.width
    bins = (np.floor(bins_f) if edge else bins_f).astype(int)
    wedge = np.full(M0, -1, np.int64)
    for p, b in enumerate(bins):
        if 0 <= b < len(rng.index_list) and rng.index_list[b] is not None:
            wedge[p] = rng.index_list[b]
    out["status"][keys] = 2  # passed eta, no frame unless set below
    valid = np.nonzero(wedge >= 0)[0]
    out["frame"][keys[valid]] = wedge[valid]
    out["status"][keys[valid]] = PRESENT
    st = cap.get("stage")
    if st is None or len(valid) == 0:
        return out
    image_list, local_wedge, det_hit, centers, vverts = st[0], st[1], st[2], st[3], st[4]
    tail = st[5:]
    ray = len(tail) == 8
    q_sum, n_q = 0.0, 0
    for m, p in enumerate(valid):
        args = [
            image_list,
            local_wedge[m : m + 1],
            det_hit[m : m + 1],
            None if centers is None else np.ascontiguousarray(centers[m : m + 1]),
            None if vverts is None else np.ascontiguousarray(vverts[m : m + 1]),
            tail[0],
            tail[1],
            1,
            tail[3],
            tail[4],
        ]
        if len(tail) > 5:
            args.append(tail[5])
        if ray:
            args += [
                np.ascontiguousarray(tail[6][m : m + 1]),
                np.ascontiguousarray(tail[7][m : m + 1]),
            ]
        r = o_stage(*args)
        k = int(keys[p])
        out["cnt"][k] = (int(r[0]), int(r[1]), int(r[2]), int(r[3]), int(r[6]))
        out["q"][k] = float(r[4])
        if ray:
            out["rays"][k] = np.asarray(r[7], dtype=np.int32)
        out["hit"][k] = det_hit[m].astype(np.uint8)
        if centers is not None:
            out["sig"][k, :, 0, :] = centers[m]
        else:
            out["sig"][k] = np.floor(vverts[m]).astype(np.int32)
        if r[1] > 0:
            q_sum += float(r[4])
            n_q += 1
    out["q_sum"], out["n_q"] = q_sum, n_q
    return out


def check_trace(tr: Dict[str, Any], tol: float = 1e-6) -> float:
    """Largest difference between the quality recombined from the peaks and the library's own."""
    n_q = int(tr.get("n_q", 0))
    q = tr.get("q_sum", 0.0) / n_q if n_q else 0.0
    agg = tr["agg"]
    assert n_q == agg["nq"], (n_q, agg["nq"])
    assert int(tr["cnt"][:, 2].sum()) == agg["pkov"], (int(tr["cnt"][:, 2].sum()), agg["pkov"])
    assert int(tr["cnt"][:, 3].sum()) == agg["pkon"]
    assert int(tr["cnt"][:, 0].sum()) == agg["pixov"]
    d = abs(q - agg["quality"])
    assert d < max(tol, 1e-5), (q, agg["quality"])
    return d


# --------------------------------------------------------------------------------------------
# events and attribution
# --------------------------------------------------------------------------------------------
def inbounds_pattern(tr: Dict[str, Any], dims: Tuple[int, int], radius: int) -> np.ndarray:
    """(P, n_det, K) bool: the vertex (or screen centre) lies inside the detector area (for the
    screen, within the search radius of it)."""
    sig = tr["sig"]  # (P, n_det, K, 2) = col, row
    nr, nc = dims
    c, r = sig[..., 0], sig[..., 1]
    return (c >= -radius) & (c < nc + radius) & (r >= -radius) & (r < nr + radius)


def classify(
    a: Dict[str, Any], b: Dict[str, Any], dims: Tuple[int, int], radius: int
) -> np.ndarray:
    """Event type code per peak between traces a and b: index into EVENT_TYPES, or -1 for no change
    (and for peaks absent in both)."""
    P = len(a["status"])
    ev = np.full(P, -1, np.int8)
    sa, sb = a["status"], b["status"]
    pa, pb = sa == PRESENT, sb == PRESENT
    # presence changes: reason from the absent side
    for code, name in ABSENT.items():
        idx = (pa & (sb == code)) | (pb & (sa == code))
        ev[idx] = EVENT_TYPES.index(name)
    both = pa & pb
    fr = both & (a["frame"] != b["frame"])
    ev[fr] = EVENT_TYPES.index("frame")
    rest = both & ~fr
    ia, ib = inbounds_pattern(a, dims, radius), inbounds_pattern(b, dims, radius)
    det = rest & ((ia != ib).any(axis=(1, 2)) | (a["hit"] != b["hit"]).any(axis=1))
    ev[det] = EVENT_TYPES.index("detector")
    rest = rest & ~det
    px = rest & (a["sig"] != b["sig"]).any(axis=(1, 2, 3))
    ev[px] = EVENT_TYPES.index("pixel")
    rest = rest & ~px
    sp = rest & (a["cnt"] != b["cnt"]).any(axis=1)
    ev[sp] = EVENT_TYPES.index("subpixel")
    rest = rest & ~sp
    ru = rest & (a["rays"] != b["rays"]).any(axis=1)
    ev[ru] = EVENT_TYPES.index("rule")
    return ev


def attribute(
    a: Dict[str, Any], b: Dict[str, Any], ev: np.ndarray
) -> Tuple[float, Dict[str, float], Dict[str, int]]:
    """Change of the cost from a to b and its split by cause.

    quality = S / N (S = sum of the per-peak quality, N = number of peaks with a pixel on a
    detector).
    dQ = (S_b - S_a) / N_b + S_a (1 / N_b - 1 / N_a): the first term is a sum over peaks (assigned
    to
    the event type of each peak), the second is the change of the denominator at fixed numerator.
    Returns (dcost, signed contribution to the cost per cause, flags): cost = 1 - quality, so every
    contribution is minus the quality term. flags: n_den_only (1 if the denominator changed and the
    numerator sum did not), n_flip (peaks whose qualified flag flipped)."""
    Sa, Na = a.get("q_sum", 0.0), int(a.get("n_q", 0))
    Sb, Nb = b.get("q_sum", 0.0), int(b.get("n_q", 0))
    Qa = Sa / Na if Na else 0.0
    Qb = Sb / Nb if Nb else 0.0
    dcost = -(Qb - Qa)
    parts: Dict[str, float] = {}
    if Nb:
        dq = b["q"] - a["q"]
        for i in np.nonzero(dq)[0]:
            name = EVENT_TYPES[ev[i]] if ev[i] >= 0 else "unattributed"
            parts[name] = parts.get(name, 0.0) - float(dq[i]) / Nb
        if Na:
            parts["denominator"] = -(Sa * (1.0 / Nb - 1.0 / Na))
    elif Na:  # everything gone
        parts["denominator"] = -(Qb - Qa)
    parts = {k: v for k, v in parts.items() if v != 0.0}
    flags = dict(
        n_den_only=int(Na != Nb and abs(Sb - Sa) < 1e-12),
        n_flip=int(((a["cnt"][:, 2] > 0) != (b["cnt"][:, 2] > 0)).sum()),
        dN=Nb - Na,
    )
    return dcost, parts, flags


def event_counts(ev: np.ndarray) -> Dict[str, int]:
    return {
        EVENT_TYPES[t]: int((ev == t).sum()) for t in range(len(EVENT_TYPES)) if (ev == t).any()
    }


def flip_causes(a: Dict[str, Any], b: Dict[str, Any], ev: np.ndarray) -> Dict[str, int]:
    """Peaks whose qualified flag (peak_overlap) flipped, counted by event type."""
    fl = (a["cnt"][:, 2] > 0) != (b["cnt"][:, 2] > 0)
    out: Dict[str, int] = {}
    for i in np.nonzero(fl)[0]:
        name = EVENT_TYPES[ev[i]] if ev[i] >= 0 else "unattributed"
        out[name] = out.get(name, 0) + 1
    return out


# --------------------------------------------------------------------------------------------
# orientations
# --------------------------------------------------------------------------------------------
def unit_directions(n: int, seed: int) -> np.ndarray:
    ax = np.random.default_rng(seed).normal(size=(n, 3))
    return ax / np.linalg.norm(ax, axis=1, keepdims=True)


def rotate_about(R: np.ndarray, radii_deg: np.ndarray, d: np.ndarray) -> np.ndarray:
    """R rotated (left-multiplied) by radius x direction d (unit 3-vector): (len(radii), 3, 3)."""
    from scipy.spatial.transform import Rotation

    rv = np.radians(radii_deg)[:, None] * d[None]
    return Rotation.from_rotvec(rv).as_matrix() @ R


def quiet() -> Any:
    return contextlib.redirect_stdout(io.StringIO())


def ray_radii(window: float = 0.05, step: float = 0.0005) -> Dict[str, np.ndarray]:
    """Radii (deg) of the fine windows along a ray. Each window is ``window`` wide at ``step``.
    Regimes: finisher (<0.05 deg), quick MC (0.05-0.5 deg), coarse grid (0.5-5 deg)."""
    starts = {
        "fin": [0.0],
        "mc": [0.05, 0.2, 0.45 - window],
        "grid": [0.5, 1.5, 3.0, 5.0 - window],
    }
    n = int(round(window / step)) + 1
    return {k: np.stack([np.linspace(s0, s0 + window, n) for s0 in ss]) for k, ss in starts.items()}


def log_radii(n: int = 120, lo: float = 0.0005, hi: float = 5.0) -> np.ndarray:
    return np.concatenate([[0.0], np.geomspace(lo, hi, n)])
