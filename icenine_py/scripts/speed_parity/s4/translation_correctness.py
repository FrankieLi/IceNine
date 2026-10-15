"""S4 Step 1: is a translated voxel's spot set the original spot set moved by A_p dx?

Subcommands (from icenine_py/; each writes benchmarks/s4_translation/step1_<name>.json):

``law``    1a and 1b on Phase D geometry. Orientations: 100 drawn from the level-0 set plus 20
           within 1 degree of the truth of 3 reference voxels. Per orientation a base voxel and
           same-type grid voxels from one voxel step to about 1 mm away (plus two out-of-layer
           moves for 1a). 1a: continuous coordinates of the local cost function's peaks (vertices
           and centroid): a float64 re-implementation of stages B and C, and the float32 path
           (the stage-D input of ``evaluate``), against u(x) + A dx; peak sets compared exactly.
           1b: lit pixels from the simulator (``ForwardSimulation._simulate_peaks_batched``
           vertices, ``ImageData.add_triangle_scanline`` raster): (i) integer shift of one
           chosen peak, (ii) sub-pixel phase re-rasterised plus the integer part, (iii) pure
           pixel-set translation by round(A dx).
``bitmap`` 1e: the dilated bitmap and centre bounds test against the C ``stage_d_overlap``
           pixel-radius branch (``pixel_radius_overlap_uint8``) for every centre pixel of 4
           frames, both detectors, including a band outside the image.
``cost``   1c: re-cost (local cost function) of every screen passer of the reference voxels,
           exact per-voxel rasterisation against the pure pixel-set translation of the spot
           rasterised once at a reference position of the same triangle type.

Usage:
  uv run python scripts/speed_parity/s4/translation_correctness.py law [--n-random 100]
  uv run python scripts/speed_parity/s4/translation_correctness.py bitmap
  uv run python scripts/speed_parity/s4/translation_correctness.py cost --part 0 --parts 2
"""

import argparse
import copy
import json
import math
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C  # noqa: E402

from icenine import cost_functions as CF  # noqa: E402
from icenine import forward_simulation as FS  # noqa: E402
from icenine.image_data import ImageData  # noqa: E402

N_NEAR_VOXELS = 3
N_NEAR = 20
STEP_TARGETS_UM = (1.0, 3.0, 50.0, 150.0, 400.0, 1000.0)  # in voxel steps for the first two
EPS_INT = 1e-9


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------


class _Rec:
    """Image sink that records the vertices the simulator rasterises."""

    def __init__(self, log: List[Tuple[int, int, Tuple[float, ...]]], f: int, d: int) -> None:
        self.log, self.f, self.d = log, f, d

    def add_triangle_scanline(self, v0: Any, v1: Any, v2: Any, inten: Any = 1.0) -> None:
        self.log.append((self.f, self.d, (v0[0], v0[1], v1[0], v1[1], v2[0], v2[1])))


class _PixImg(ImageData):
    """Sparse ImageData whose pixel writes are collected: the simulator's raster, unchanged."""

    def __init__(self, rows: int, cols: int) -> None:
        super().__init__(rows, cols, mode="sparse")
        self.px: Set[Tuple[int, int]] = set()

    def add_to_pixel(self, j: Any, k: Any, value: Any) -> None:  # type: ignore[override]
        self.px.add((int(j), int(k)))


def raster(v: np.ndarray, rows: int, cols: int) -> Set[Tuple[int, int]]:
    """Pixels (col, row) ImageData.add_triangle_scanline lights for vertices v (3, 2)."""
    im = _PixImg(rows, cols)
    im.add_triangle_scanline(
        (float(v[0, 0]), float(v[0, 1])),
        (float(v[1, 0]), float(v[1, 1])),
        (float(v[2, 0]), float(v[2, 1])),
        1.0,
    )
    return im.px


def raster_iv(iv: np.ndarray, rows: int, cols: int) -> Set[Tuple[int, int]]:
    """The raster of integer vertices iv (3, 2) (already through the vertex rule). They are passed
    as iv + 0.5, which add_triangle_scanline truncates back to iv exactly (its float32 cast of a
    float64 coordinate could otherwise round a value just under an integer up to it)."""
    return raster(np.asarray(iv, np.float64) + 0.5, rows, cols)


FREE_OFFSET = 8192  # integer offset that puts any spot well inside a virtual detector


def raster_free(v: np.ndarray) -> Set[Tuple[int, int]]:
    """The same raster without the detector's clip: rasterised on a virtual 3 x 8192 px detector
    after an integer offset (exact by the rasteriser's integer translation property, which 1b(i)
    tests on the real detector), then shifted back."""
    big = 3 * FREE_OFFSET
    s = raster_iv(np.floor(v) + FREE_OFFSET, big, big)  # floor: no detector-edge rule
    return shift_set(s, np.array([-FREE_OFFSET, -FREE_OFFSET]))


def shift_set(s: Set[Tuple[int, int]], m: np.ndarray) -> Set[Tuple[int, int]]:
    a, b = int(m[0]), int(m[1])
    return {(c + a, r + b) for c, r in s}


def trunc_c(v: np.ndarray) -> np.ndarray:
    """The rasteriser's vertex rule: -1 for a negative coordinate, else floor."""
    return np.where(v < 0, -1.0, np.floor(v))


def int_dist(v: np.ndarray) -> float:
    """Smallest distance of any coordinate to an integer (pixel boundary)."""
    return float(np.min(np.abs(v - np.round(v))))


def flip_dist(v_got: np.ndarray, v_pred: np.ndarray) -> float:
    """Float64 distance to the nearest pixel boundary of the predicted coordinates whose integer
    vertex differs (the largest such distance; 0 if none differ)."""
    d = trunc_c(v_got) != trunc_c(v_pred)
    if not d.any():
        return 0.0
    return float(np.max(np.abs(v_pred[d] - np.round(v_pred[d]))))


def clipped(iv: np.ndarray, rows: int, cols: int) -> bool:
    """True if the Sutherland-Hodgman clip of the integer vertices changes the polygon (some vertex
    outside [0, cols-1) x [0, rows-1), the clipper's inside test)."""
    return bool(
        ((iv[:, 0] < 0) | (iv[:, 0] >= cols - 1) | (iv[:, 1] < 0) | (iv[:, 1] >= rows - 1)).any()
    )


def adjacent_ok(
    diff: Set[Tuple[int, int]], A: Set[Tuple[int, int]], B: Set[Tuple[int, int]]
) -> bool:
    """Every differing pixel has an 8-neighbour (or itself) in the other set: the two sets differ
    only in a one-pixel layer along their boundaries."""
    for q in diff:
        other = B if q in A else A
        c, r = q
        if not any((c + i, r + j) in other for i in (-1, 0, 1) for j in (-1, 0, 1)):
            return False
    return True


def boundary_ok(
    diff: Set[Tuple[int, int]], A: Set[Tuple[int, int]], B: Set[Tuple[int, int]]
) -> bool:
    """Every differing pixel lies on the boundary of the set it belongs to (a 4-neighbour outside)."""
    for q in diff:
        S = A if q in A else B
        c, r = q
        if all(n in S for n in ((c + 1, r), (c - 1, r), (c, r + 1), (c, r - 1))):
            return False
    return True


def f64_geometry(cf: Any, omegas: np.ndarray, normals: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Float64 stage B from float32 (omega, crystal-plane normal): full_rot (P,3,3), refl (P,3)."""
    B = cf.sample.sample_to_lab_matrix[:3, :3].double().numpy()
    b = cf.simulator.beam_direction.double().numpy()
    w = omegas.astype(np.float64)
    Rz = np.zeros((len(w), 3, 3))
    Rz[:, 0, 0], Rz[:, 0, 1], Rz[:, 1, 0], Rz[:, 1, 1], Rz[:, 2, 2] = (
        np.cos(w),
        -np.sin(w),
        np.sin(w),
        np.cos(w),
        1.0,
    )
    full = Rz @ B
    ln = np.einsum("pij,pj->pi", full, normals.astype(np.float64))
    refl = b[None, :] - 2.0 * (ln @ b)[:, None] * ln
    return full, refl


class Spy:
    """Captures the peaks and the stage-D arguments of VoxelCostFunction.evaluate."""

    def __init__(self) -> None:
        self.cap: Dict[str, Any] = {}
        self._calc = CF.calculate_diffraction_overlap_batched
        self._stage = CF._c_stage_d_overlap
        self.replace_verts: Optional[np.ndarray] = None

    def __enter__(self) -> "Spy":
        spy = self

        def calc(**kw: Any) -> Any:
            spy.cap["omegas"] = kw["peak_omegas"].numpy().copy()
            spy.cap["normals"] = kw["peak_normals"].numpy().copy()
            return spy._calc(**kw)

        def stage(*a: Any) -> Any:
            spy.cap["wedge"] = np.array(a[1]).copy()
            spy.cap["det_hit"] = np.array(a[2]).copy()
            spy.cap["verts"] = None if a[4] is None else np.array(a[4], np.float64).copy()
            spy.cap["centers"] = None if a[3] is None else np.array(a[3]).copy()
            spy.cap["args"] = a
            return spy._stage(*a)

        CF.calculate_diffraction_overlap_batched = calc  # type: ignore[assignment]
        CF._c_stage_d_overlap = stage  # type: ignore[assignment]
        return self

    def __exit__(self, *exc: Any) -> None:
        CF.calculate_diffraction_overlap_batched = self._calc  # type: ignore[assignment]
        CF._c_stage_d_overlap = self._stage  # type: ignore[assignment]

    def run(self, cf: Any, R: np.ndarray, verts: torch.Tensor, phase: int) -> Dict[str, Any]:
        self.cap = {}
        info = cf.evaluate(orientation=R, voxel_vertices=verts, phase_index=phase)
        out = dict(self.cap)
        out["info"] = info
        return out


def frame_valid_mask(cf: Any, omegas: np.ndarray) -> np.ndarray:
    idx = CF.frame_index_array(cf.range_map)
    fr = CF.frame_of_omegas(omegas, cf.range_map.low, cf.range_map.width, idx, cf.cpp_edge_rules)
    return fr


# ---------------------------------------------------------------------------------------------
# law: 1a and 1b
# ---------------------------------------------------------------------------------------------


def sim_setup() -> Dict[str, Any]:
    """The Phase D renderer's simulator (scripts/phase_d/render_full.py _worker)."""
    from icenine.config_file import ConfigFile
    from icenine.forward_simulation import ForwardSimulation
    from icenine.sample import Sample
    from icenine.simulation import Simulation

    cfg = ConfigFile.from_file(str(C.EX_D / "ConfigFiles" / "Example2.Simulation.config"))
    cfg.sample_filename = str(C.TRUTH_D)
    with C.quiet_stdout():
        sim = ForwardSimulation(cfg)
        sim.exp_setup.initialize_experiment()
        dets = sim.exp_setup.get_detector_list()
        rmap = sim.exp_setup.get_range_to_index_map()
        sim.simulator = Simulation(sim.exp_setup)
        sample = Sample()
        sim.exp_setup.initialize_sample(sample, dets[0])
    mic = sample.get_mic()
    return dict(sim=sim, dets=dets, rmap=rmap, sample=sample, mic=mic, template=mic.voxels[0])


def simulate_one(S: Dict[str, Any], R: np.ndarray, voxel: Any, pos: np.ndarray) -> List[Any]:
    """Recorded peaks of one voxel at ``pos``: list of (frame, verts (n_det, 3, 2) float64)."""
    v = copy.copy(S["template"])
    v.orientation = np.asarray(R, np.float64)
    v.position = np.asarray(pos, np.float64)
    v.side_length = float(voxel.side_length)
    v.points_up = bool(voxel.points_up)
    S["mic"].voxels = [v]
    log: List[Tuple[int, int, Tuple[float, ...]]] = []
    n_det = len(S["dets"])
    n_om = len(S["sim"].exp_setup.get_omega_range_list())
    images = [[_Rec(log, i, d) for d in range(n_det)] for i in range(n_om)]
    with C.quiet_stdout():
        S["sim"]._simulate_peaks_batched(images, S["dets"], S["sample"], S["rmap"], batch_size=None)
    peaks = []
    for k in range(0, len(log), n_det):
        grp = log[k : k + n_det]
        assert [g[1] for g in grp] == list(range(n_det)) and len({g[0] for g in grp}) == 1
        peaks.append((grp[0][0], np.array([g[2] for g in grp], np.float64).reshape(n_det, 3, 2)))
    return peaks


def pick_displacements(ds: C.DataSet, base: int, rng: np.random.Generator) -> List[int]:
    """Same-type grid voxels at about 1 and 3 voxel steps, 50, 150, 400 um and the farthest."""
    P = np.array([np.asarray(v.position, float) for v in ds.voxels])
    up = np.array([bool(v.points_up) for v in ds.voxels])
    same = np.nonzero(up == up[base])[0]
    same = same[same != base]
    d = np.linalg.norm(P[same] - P[base], axis=1)
    s = float(ds.voxels[base].side_length)
    out = []
    for i, tgt in enumerate(STEP_TARGETS_UM):
        t = tgt * s if i < 2 else tgt * 1e-3
        if i == len(STEP_TARGETS_UM) - 1:
            cand = same[d >= d.max() - 2 * s]
        else:
            cand = same[np.abs(d - t) <= max(0.5 * s, 0.02 * t)]
        if len(cand) == 0:
            cand = same[[int(np.argmin(np.abs(d - t)))]]
        out.append(int(rng.choice(cand)))
    return out


def match_sim_to_cf(sim_peaks: List[Any], frames: np.ndarray, U64: np.ndarray) -> np.ndarray:
    """Index of the cost-function peak each simulator peak is (same frame, every vertex within
    0.05 px on every detector), -1 if none. U64: (P, n_det, 3, 2) float64 cost-function vertices."""
    out = np.full(len(sim_peaks), -1, np.int64)
    for i, (f, v) in enumerate(sim_peaks):
        cand = np.nonzero(frames == f)[0]
        if len(cand) == 0:
            continue
        err = np.abs(U64[cand] - v[None]).reshape(len(cand), -1).max(axis=1)
        j = int(np.argmin(err))
        if err[j] < 0.05:
            out[i] = cand[j]
    return out


def compare_sets(
    got: Set[Tuple[int, int]],
    want: Set[Tuple[int, int]],
    iv_got: np.ndarray,
    iv_want: np.ndarray,
    src_iv: List[np.ndarray],
    rows: int,
    cols: int,
) -> str:
    """'equal', or the class of a difference: 'edge' (clipped at the detector), 'float' (an integer
    vertex differs: a float flip), 'fail' (same integer vertices, no clipping, different sets)."""
    if got == want:
        return "equal"
    if any(clipped(x, rows, cols) for x in [iv_got, iv_want] + src_iv):
        return "edge"
    if not np.array_equal(iv_got, iv_want):
        return "float"
    return "fail"


def cmd_law(a: argparse.Namespace) -> None:
    ds = C.load_set("phase_d")
    ref_ids = C.reference_voxels(ds)
    cap = C.capture_level0(ds, ds.voxels[ref_ids[0]])
    lcf = cap["local_cost_fn"]
    geoms = C.det_geoms(lcf)
    t_s = C.sample_translation(lcf)
    rows, cols = geoms[0].rows, geoms[0].cols
    S = sim_setup()
    rng = np.random.default_rng(C.REF_SEED + 7)
    R_all = C.all_orientations(cap)
    cases: List[Tuple[np.ndarray, int, str]] = []
    for k in rng.choice(len(R_all), a.n_random, replace=False):
        cases.append((R_all[int(k)], int(ref_ids[len(cases) % len(ref_ids)]), "random"))
    from scipy.spatial.transform import Rotation

    for vid in ref_ids[:N_NEAR_VOXELS]:
        for _ in range(a.n_near):
            ax = rng.normal(size=3)
            ax /= np.linalg.norm(ax)
            dR = Rotation.from_rotvec(np.radians(rng.uniform(0, 1.0)) * ax).as_matrix()
            cases.append((dR @ ds.truth[vid], vid, "near"))

    res: Dict[str, Any] = dict(
        n_orientations=len(cases),
        n_random=a.n_random,
        n_near=a.n_near * N_NEAR_VOXELS,
        dev64_max=0.0,
        dev32_max=0.0,
        base32_max=0.0,
        dev64_where=None,
        dev32_where=None,
        peakset_identical_cf=0,
        peakset_compared_cf=0,
        peakset_identical_sim=0,
        peakset_compared_sim=0,
        sim_matched=0,
        sim_unmatched=0,
        sim_vs_cf64_max=0.0,
    )
    cls_keys = ("equal", "edge", "edge_off", "float", "fail")
    tallies: Dict[str, Dict[str, int]] = {
        k: {c: 0 for c in cls_keys} for k in ("i32", "i64", "ii32", "ii64")
    }
    flips: Dict[str, List[float]] = {"i32": [], "ii32": [], "i64": [], "ii64": []}
    fails: List[Dict[str, Any]] = []
    iii_rel: List[float] = []
    iii_abs: List[int] = []
    iii_boundary_ok = 0
    iii_adjacent_ok = 0
    iii_interior_examples: List[Dict[str, Any]] = []
    iii_n = 0
    iii_equal = 0
    iii_edge = 0
    sim_peakset_differs = 0
    n_leave = 0
    n_enter = 0
    disp_um: List[float] = []
    t0 = time.time()
    with Spy() as spy:
        for ci, (R, base, kind) in enumerate(cases):
            vb = ds.voxels[base]
            Xb = C.voxel_vertices(vb)
            others = pick_displacements(ds, base, rng)
            pos_b = np.asarray(vb.position, float)
            c0 = spy.run(lcf, R, Xb, vb.phase)
            if "verts" not in c0 or c0["verts"] is None:
                continue
            om, nr = c0["omegas"], c0["normals"]
            fr_all = frame_valid_mask(lcf, om)
            keep = fr_all >= 0
            full, refl = f64_geometry(lcf, om[keep], nr[keep])
            A, cc, par, la, lb = C.affine_maps(full, refl, t_s, geoms)
            X64 = Xb.double().numpy()
            U_b, _ = C.project_direct64(full, refl, t_s, geoms, np.vstack([X64, X64.mean(0)]))
            # (P, nd, 4, 2): 3 vertices + centroid
            res["base32_max"] = max(
                res["base32_max"], float(np.abs(c0["verts"] - U_b[:, :, :3, :]).max())
            )
            frames_cf = fr_all[keep]
            assert np.array_equal(frames_cf, c0["wedge"])
            # simulator at the base
            sp_b = simulate_one(S, R, vb, pos_b)
            midx = match_sim_to_cf(sp_b, frames_cf, U_b[:, :, :3, :])
            res["sim_matched"] += int((midx >= 0).sum())
            res["sim_unmatched"] += int((midx < 0).sum())
            for i, (f, v) in enumerate(sp_b):
                if midx[i] >= 0:
                    res["sim_vs_cf64_max"] = max(
                        res["sim_vs_cf64_max"], float(np.abs(v - U_b[midx[i], :, :3, :]).max())
                    )
            S_b = {
                (i, d): raster(sp_b[i][1][d], rows, cols)
                for i in range(len(sp_b))
                if midx[i] >= 0
                for d in range(len(geoms))
            }
            S_b64 = {
                (i, d): raster_iv(trunc_c(U_b[midx[i], d, :3, :]), rows, cols)
                for i in range(len(sp_b))
                if midx[i] >= 0
                for d in range(len(geoms))
            }
            moves: List[Tuple[str, np.ndarray, Optional[Tuple[int, int, np.ndarray]]]] = []
            for o in others:
                moves.append(("grid", np.asarray(ds.voxels[o].position, float) - pos_b, None))
            for dz in (0.02, -0.05):  # out-of-layer (1a only)
                dxz = np.array([rng.uniform(-0.1, 0.1), rng.uniform(-0.1, 0.1), dz])
                moves.append(("z", dxz, None))
            # (i): integer shift of a chosen peak on a chosen detector
            matched = [i for i in range(len(sp_b)) if midx[i] >= 0]
            for i in rng.choice(matched, min(a.n_int, len(matched)), replace=False):
                d = int(rng.integers(len(geoms)))
                A2 = A[midx[i], d][:, :2]
                dx0 = rng.uniform(-0.3, 0.3, 2)
                m = np.round(A2 @ dx0)
                dxy = np.linalg.solve(A2, m)
                moves.append(("int", np.array([dxy[0], dxy[1], 0.0]), (int(i), d, m)))
            for kind_m, dx, chosen in moves:
                disp_um.append(float(np.linalg.norm(dx) * 1e3))
                Xn = Xb + torch.from_numpy(dx).float()  # the cost function's float32 vertices
                X64n = X64 + dx
                c1 = spy.run(lcf, R, Xn, vb.phase)
                res["peakset_compared_cf"] += 1
                same = (
                    np.array_equal(c1["omegas"], om)
                    and np.array_equal(c1["normals"], nr)
                    and np.array_equal(c1["wedge"], c0["wedge"])
                    and np.array_equal(c1["det_hit"], c0["det_hit"])
                )
                res["peakset_identical_cf"] += int(same)
                if not same:
                    continue
                U_n, _ = C.project_direct64(full, refl, t_s, geoms, np.vstack([X64n, X64n.mean(0)]))
                pred = U_b + np.einsum("pdaj,j->pda", A, dx)[:, :, None, :]
                e64 = np.abs(U_n - pred)
                e32 = np.abs(c1["verts"] - pred[:, :, :3, :])
                if e64.max() > res["dev64_max"]:
                    p, d, q, ax = np.unravel_index(int(np.argmax(e64)), e64.shape)
                    res["dev64_max"] = float(e64.max())
                    res["dev64_where"] = dict(
                        case=ci,
                        move=kind_m,
                        disp_um=disp_um[-1],
                        peak=int(p),
                        det=int(d),
                        coord=float(U_n[p, d, q, ax]),
                    )
                if e32.max() > res["dev32_max"]:
                    p, d, q, ax = np.unravel_index(int(np.argmax(e32)), e32.shape)
                    res["dev32_max"] = float(e32.max())
                    res["dev32_where"] = dict(
                        case=ci,
                        move=kind_m,
                        disp_um=disp_um[-1],
                        peak=int(p),
                        det=int(d),
                        coord=float(c1["verts"][p, d, q, ax]),
                    )
                if kind_m == "z":
                    continue
                # 1b: the simulator at the moved position
                sp_n = simulate_one(S, R, vb, pos_b + dx)
                res["peakset_compared_sim"] += 1
                if len(sp_n) != len(sp_b) or any(
                    sp_n[i][0] != sp_b[i][0] for i in range(len(sp_b))
                ):
                    sim_peakset_differs += 1
                    continue
                res["peakset_identical_sim"] += 1
                for i in matched:
                    p = int(midx[i])
                    for d in range(len(geoms)):
                        if chosen is not None and (i, d) != (chosen[0], chosen[1]):
                            continue
                        du = A[p, d] @ dx  # float64 continuous shift
                        v_b32, v_n32 = sp_b[i][1][d], sp_n[i][1][d]
                        v_b64, v_n64 = U_b[p, d, :3, :], U_n[p, d, :3, :]
                        got32 = raster(v_n32, rows, cols)
                        if chosen is None:
                            n_leave += int(bool(S_b[(i, d)]) and not got32)
                            n_enter += int(not S_b[(i, d)] and bool(got32))
                        got64 = raster_iv(trunc_c(v_n64), rows, cols)
                        if chosen is not None:
                            m = chosen[2]
                            for key, got, vb_, vn_, Sb in (
                                ("i32", got32, v_b32, v_n32, S_b[(i, d)]),
                                ("i64", got64, v_b64, v_n64, S_b64[(i, d)]),
                            ):
                                want = shift_set(Sb, m)
                                cl = compare_sets(
                                    got,
                                    want,
                                    trunc_c(vn_),
                                    trunc_c(vb_) + m,
                                    [trunc_c(vb_)],
                                    rows,
                                    cols,
                                )
                                tallies[key][cl] += 1
                                if cl == "edge" and not got:
                                    tallies[key][
                                        "edge_off"
                                    ] += 1  # (a subset of edge: no pixel on the detector at x')
                                if cl == "float":
                                    flips[key].append(flip_dist(vn_, vb_ + du))
                                if cl == "fail":
                                    fails.append(dict(test=key, case=ci, peak=i, det=d))
                            continue
                        m = np.floor(du)
                        phi = du - m
                        for key, got, vb_, vn_ in (
                            ("ii32", got32, v_b32, v_n32),
                            ("ii64", got64, v_b64, v_n64),
                        ):
                            want = shift_set(raster_free(vb_ + phi), m)
                            cl = compare_sets(
                                got, want, trunc_c(vn_), np.floor(vb_ + phi) + m, [], rows, cols
                            )
                            tallies[key][cl] += 1
                            if cl == "edge" and not got:
                                tallies[key][
                                    "edge_off"
                                ] += 1  # (a subset of edge: no pixel on the detector at x')
                            if cl == "float":
                                flips[key].append(flip_dist(vn_, vb_ + du))
                            if cl == "fail":
                                fails.append(
                                    dict(
                                        test=key,
                                        case=ci,
                                        peak=i,
                                        det=d,
                                        move=kind_m,
                                        v_base=vb_.tolist(),
                                        v_moved=vn_.tolist(),
                                        du=du.tolist(),
                                        got=sorted(got),
                                        want=sorted(want),
                                    )
                                )
                        # (iii) pure pixel-set translation by round(du)
                        T = shift_set(S_b[(i, d)], np.round(du))
                        if clipped(trunc_c(v_n32), rows, cols) or clipped(
                            trunc_c(v_b32), rows, cols
                        ):
                            iii_edge += 1
                            continue
                        diff = T ^ got32
                        iii_n += 1
                        iii_equal += int(not diff)
                        iii_abs.append(len(diff))
                        iii_rel.append(len(diff) / max(len(got32), 1))
                        ok_b = boundary_ok(diff, got32, T)
                        iii_boundary_ok += int(ok_b)
                        iii_adjacent_ok += int(adjacent_ok(diff, got32, T))
                        if not ok_b and len(iii_interior_examples) < 5:
                            iii_interior_examples.append(
                                dict(
                                    case=ci,
                                    peak=i,
                                    det=d,
                                    moved=sorted(got32),
                                    translated=sorted(T),
                                )
                            )
            if ci % 10 == 0:
                print(f"case {ci}/{len(cases)} {time.time() - t0:.0f}s", flush=True)
    rel = np.array(iii_rel)
    res.update(
        tallies=tallies,
        float_flip_dist_max={k: (max(v) if v else None) for k, v in flips.items()},
        float_flip_count={k: len(v) for k, v in flips.items()},
        fails=fails[:50],
        n_fails=len(fails),
        iii=dict(
            n_spots=iii_n,
            n_equal=iii_equal,
            n_edge_excluded=iii_edge,
            symdiff_px_median=float(np.median(iii_abs)) if iii_abs else None,
            symdiff_px_p95=float(np.percentile(iii_abs, 95)) if iii_abs else None,
            symdiff_px_max=int(max(iii_abs)) if iii_abs else None,
            rel_median=float(np.median(rel)) if len(rel) else None,
            rel_p95=float(np.percentile(rel, 95)) if len(rel) else None,
            rel_max=float(rel.max()) if len(rel) else None,
            rel_hist_edges=[0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 10.0],
            rel_hist=(
                np.histogram(rel, [0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 10.0])[0].tolist()
                if len(rel)
                else []
            ),
            all_on_boundary=iii_boundary_ok,
            all_adjacent_to_other_set=iii_adjacent_ok,
            interior_examples=iii_interior_examples,
        ),
        moves_sim_peakset_differs=sim_peakset_differs,
        spots_leaving_detector=n_leave,
        spots_entering_detector=n_enter,
        displacement_um_min=min(disp_um),
        displacement_um_max=max(disp_um),
        wall_s=time.time() - t0,
    )
    stop = (
        res["dev64_max"] > 1e-9
        or tallies["i32"]["fail"]
        or tallies["i64"]["fail"]
        or tallies["ii32"]["fail"]
        or tallies["ii64"]["fail"]
        or res["peakset_identical_cf"] != res["peakset_compared_cf"]
    )
    res["stop_rule_triggered"] = bool(stop)
    out = C.OUT / "step1_law.json"
    out.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps({k: v for k, v in res.items() if k != "fails"}, indent=1, default=float))


# ---------------------------------------------------------------------------------------------
# bitmap: 1e
# ---------------------------------------------------------------------------------------------


def cmd_bitmap(a: argparse.Namespace) -> None:
    ds = C.load_set("phase_d")
    cap = C.capture_level0(ds, ds.voxels[C.reference_voxels(ds)[0]])
    gcf = cap["global_cost_fn"]
    r = int(gcf.pixel_radius)
    band = r + 2
    frames = [0, 45, 90, 179]
    n_det = len(gcf.detector_list)
    res: Dict[str, Any] = dict(frames=frames, radius=r, band_outside=band, images=[])
    t0 = time.time()
    for f in frames:
        for d in range(n_det):
            img = gcf.exp_data.get_image(f, d).get_binary_numpy()
            rows, cols = img.shape
            dil = C.dilate(img, r)
            cys = np.arange(-band, rows + band)
            cxs = np.arange(-band, cols + band)
            n_calls = 0
            n_pix = 0
            bad = 0
            n_lit = 0
            n_bright = 0
            for cy in cys:
                inb_y = -r <= cy < rows + r
                inb = inb_y & (cxs >= -r) & (cxs < cols + r)
                br = np.zeros(len(cxs), bool)
                if inb_y:
                    sel = inb
                    br[sel] = dil[cy + r, cxs[sel] + r] != 0
                code = inb.astype(np.int8) * 2 + br.astype(np.int8)
                cuts = np.nonzero(np.diff(code))[0] + 1
                starts = np.concatenate([[0], cuts])
                ends = np.concatenate([cuts, [len(cxs)]])
                for s, e in zip(starts, ends):
                    M = int(e - s)
                    centers = np.zeros((M, 1, 2), np.int32)
                    centers[:, 0, 0] = cxs[s:e]
                    centers[:, 0, 1] = cy
                    out = CF._c_stage_d_overlap(
                        [img],
                        np.zeros(M, np.int32),
                        np.ones((M, 1), np.uint8),
                        centers,
                        None,
                        1,
                        r,
                        M,
                        rows,
                        cols,
                    )
                    want_lit = M * int(inb[s])
                    want_br = M * int(br[s])
                    bad += int(out[3] != want_lit) + int(out[2] != want_br)
                    n_calls += 1
                    n_pix += M
                n_lit += int(inb.sum())
                n_bright += int(br.sum())
            res["images"].append(
                dict(
                    frame=f,
                    det=d,
                    centres=n_pix,
                    runs=n_calls,
                    mismatched_runs=bad,
                    in_bounds=n_lit,
                    bright=n_bright,
                    raw_lit=int((img != 0).sum()),
                )
            )
            print(json.dumps(res["images"][-1]), flush=True)
    res["all_equal"] = all(im["mismatched_runs"] == 0 for im in res["images"])
    res["centres_total"] = sum(im["centres"] for im in res["images"])
    res["wall_s"] = time.time() - t0
    res["stop_rule_triggered"] = not res["all_equal"]
    (C.OUT / "step1_bitmap.json").write_text(json.dumps(res, indent=1) + "\n")
    print("all equal:", res["all_equal"])


# ---------------------------------------------------------------------------------------------
# cost: 1c
# ---------------------------------------------------------------------------------------------


def ref_position_voxel(ds: C.DataSet, up: bool) -> int:
    """The voxel of a triangle type nearest the sample centre (the spot shapes' reference)."""
    P = np.array([np.asarray(v.position, float) for v in ds.voxels])
    ups = np.array([bool(v.points_up) for v in ds.voxels])
    cand = np.nonzero(ups == up)[0]
    return int(cand[np.argmin(np.linalg.norm(P[cand, :2], axis=1))])


def cmd_cost(a: argparse.Namespace) -> None:
    from icenine.orientation_search import SearchCandidate, _spacing_filter

    ds = C.load_set("phase_d")
    ref_ids = C.reference_voxels(ds)
    # this part's voxels, grouped by triangle type (the reference-position evaluations are cached)
    by_type = sorted(ref_ids, key=lambda k: (bool(ds.voxels[k].points_up), k))
    vids = [int(x) for x in np.array_split(np.array(by_type), a.parts)[a.part]]
    if a.limit:
        vids = vids[: a.limit]
    cap = C.capture_level0(ds, ds.voxels[vids[0]])
    lcf = cap["local_cost_fn"]
    R_all = C.all_orientations(cap)
    n_local = len(cap["local_grid"])
    out_rows = []
    ref_cache: Dict[Tuple[bool, int], Dict[str, Any]] = {}
    t0 = time.time()
    with Spy() as spy:
        for vid in vids:
            v = ds.voxels[vid]
            z = np.load(C.CACHE / "ref" / f"phase_d_{vid}.npz")
            passers = z["passers"]
            recost_ref = z["recost"]
            fields = list(z["fields"])
            up = bool(v.points_up)
            vr = ds.voxels[ref_position_voxel(ds, up)]
            Xr, Xv = C.voxel_vertices(vr), C.voxel_vertices(v)
            ex = np.zeros((len(passers), len(fields)))
            tr = np.zeros((len(passers), len(fields)))
            n_edge = 0
            n_peaks_changed = 0
            n_peaks = 0
            for i, k in enumerate(passers):
                key = (up, int(k))
                if key not in ref_cache:
                    c_r = spy.run(lcf, R_all[k], Xr, vr.phase)
                    vr32 = c_r.get("verts")  # float32 values, kept as float32 (exact)
                    ref_cache[key] = dict(
                        verts=None if vr32 is None else vr32.astype(np.float32),
                        wedge=c_r.get("wedge"),
                    )
                cr = ref_cache[key]
                c_x = spy.run(lcf, R_all[k], Xv, v.phase)
                info = c_x["info"]
                ex[i] = [float(getattr(info, f)) for f in fields]
                if cr["verts"] is None or c_x.get("verts") is None:
                    tr[i] = ex[i]
                    continue
                assert np.array_equal(cr["wedge"], c_x["wedge"])
                vref = cr["verts"]  # (M, nd, 3, 2)
                shift = np.round((c_x["verts"] - vref).mean(axis=2))  # (M, nd, 2)
                vt = trunc_c(vref) + shift[:, :, None, :]
                n_edge += int((vt < 0).any(axis=(2, 3)).sum() + (vt >= 2047).any(axis=(2, 3)).sum())
                n_peaks += int(vt.shape[0] * vt.shape[1])
                n_peaks_changed += int(
                    (np.abs(trunc_c(c_x["verts"]) - vt).max(axis=(2, 3)) > 0).sum()
                )
                args = list(c_x["args"])
                args[4] = np.ascontiguousarray(vt.astype(np.float32))
                res = spy._stage(*args)
                tr[i] = [
                    float(res[0]),
                    float(res[1]),
                    float(res[2]),
                    float(res[3]),
                    float(res[4]),
                    float(res[5]),
                    float(res[6]),
                ]
            assert np.array_equal(ex, recost_ref), "exact re-cost differs from the Step 0 reference"
            fi = {f: j for j, f in enumerate(fields)}

            def conf(m: np.ndarray) -> np.ndarray:
                pon = m[:, fi["peak_on_detector"]]
                return np.where(pon > 0, m[:, fi["peak_overlap"]] / np.maximum(pon, 1), 0.0)

            rc_ex = np.where(ex[:, fi["peak_on_detector"]] > 0, 1.0 - conf(ex), 1.0)
            rc_tr = np.where(tr[:, fi["peak_on_detector"]] > 0, 1.0 - conf(tr), 1.0)
            q_ex, q_tr = ex[:, fi["quality"]], tr[:, fi["quality"]]
            # downstream: the spacing filter of each clique, exact vs translated re-cost
            sym = cap["symmetry_quats"]
            rad = float(cap["angular_radius"])
            clique = passers // n_local
            same_filter = 0
            n_cliques = 0
            kept_ex, kept_tr = [], []
            for cq in np.unique(clique):
                idx = np.nonzero(clique == cq)[0]
                ce = [
                    SearchCandidate(
                        orientation=R_all[passers[j]], cost=float(rc_ex[j]), clique=int(cq)
                    )
                    for j in idx
                ]
                ct = [
                    SearchCandidate(
                        orientation=R_all[passers[j]], cost=float(rc_tr[j]), clique=int(cq)
                    )
                    for j in idx
                ]
                fe = _spacing_filter(ce, rad, sym)
                ft = _spacing_filter(ct, rad, sym)
                ke = sorted(int(passers[idx[[id(c) for c in ce].index(id(x))]]) for x in fe)
                kt = sorted(int(passers[idx[[id(c) for c in ct].index(id(x))]]) for x in ft)
                same_filter += int(ke == kt)
                n_cliques += 1
                kept_ex += [(float(rc_ex[np.nonzero(passers == k)[0][0]]), k) for k in ke]
                kept_tr += [(float(rc_tr[np.nonzero(passers == k)[0][0]]), k) for k in kt]
            kept_ex.sort()
            kept_tr.sort()
            top = 30
            top_ex = [k for _, k in kept_ex[:top]]
            top_tr = [k for _, k in kept_tr[:top]]
            row = dict(
                vid=int(vid),
                points_up=up,
                ref_voxel=int(ref_position_voxel(ds, up)),
                disp_um=float(
                    np.linalg.norm(np.asarray(v.position) - np.asarray(vr.position)) * 1e3
                ),
                n_passers=int(len(passers)),
                peaks_compared=n_peaks,
                peaks_integer_vertices_changed=n_peaks_changed,
                peaks_near_edge=n_edge,
                passers_any_field_changed=int((np.abs(ex - tr).max(axis=1) > 0).sum()),
                passers_peak_overlap_changed=int(
                    (ex[:, fi["peak_overlap"]] != tr[:, fi["peak_overlap"]]).sum()
                ),
                passers_pixel_overlap_changed=int(
                    (ex[:, fi["pixel_overlap"]] != tr[:, fi["pixel_overlap"]]).sum()
                ),
                recost_abs_diff_max=float(np.abs(rc_ex - rc_tr).max()),
                recost_abs_diff_median_changed=(
                    float(np.median(np.abs(rc_ex - rc_tr)[rc_ex != rc_tr]))
                    if (rc_ex != rc_tr).any()
                    else 0.0
                ),
                passers_recost_changed=int((rc_ex != rc_tr).sum()),
                quality_abs_diff_max=float(np.abs(q_ex - q_tr).max()),
                quality_abs_diff_median=float(np.median(np.abs(q_ex - q_tr))),
                cliques=n_cliques,
                cliques_same_filtered_set=same_filter,
                kept_total_exact=len(kept_ex),
                kept_total_translated=len(kept_tr),
                best_same=bool(top_ex[:1] == top_tr[:1]),
                top30_overlap=len(set(top_ex) & set(top_tr)),
            )
            out_rows.append(row)
            print(json.dumps(row), flush=True)
    path = C.OUT / "cache" / f"step1_cost_part{a.part}.json"
    path.write_text(json.dumps(dict(rows=out_rows, wall_s=time.time() - t0), indent=1) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("law")
    p.add_argument("--n-random", type=int, default=100)
    p.add_argument("--n-near", type=int, default=N_NEAR)
    p.add_argument("--n-int", type=int, default=3, help="integer-shift peaks per orientation")
    p = sub.add_parser("bitmap")
    p = sub.add_parser("cost")
    p.add_argument("--part", type=int, default=0)
    p.add_argument("--parts", type=int, default=1)
    p.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    C.OUT.mkdir(parents=True, exist_ok=True)
    (C.OUT / "cache").mkdir(exist_ok=True)
    dict(law=cmd_law, bitmap=cmd_bitmap, cost=cmd_cost)[a.cmd](a)


if __name__ == "__main__":
    main()
