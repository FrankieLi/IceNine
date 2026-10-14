"""Shared helpers of the Task A step-0 parity probe (fixed-orientation C++ vs Python cost).

Python side: set up the Phase D clean full-sample reconstruction exactly as the no-start search
does, build the two cost functions the level-0 stage uses (screen: pixel radius 3, centre point,
Q-max 5; complete cost: pixel radius 0, triangle) with the C++ overlap rule, and evaluate one
orientation with a per-peak record. No library code is changed: the per-peak record is taken by
wrapping the module functions of ``icenine.cost_functions`` for one call.

C++ side: parser of the text written by the scratch driver (see
benchmarks/eval_reduction/A/step0/cpp_parity_patch.py.txt).
"""

import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
EX = ICE.parent / "Examples" / "Example2.ManyGrains"
CONFIG = ICE / "scripts" / "phase_d" / "configs" / "ReconstructPhaseD_mc_clean.config"
STACK = EX / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
TRUTH_MIC = EX / "SimInput" / "rand_500grains_1mm_neworient_s0.mic"
SIX_VOXELS = [15901, 17034, 18558, 19889, 22867, 5242]
PARITY_VOXELS = [15901, 22867, 5242]
SCREEN_QMAX = 5.0
SCREEN_RADIUS = 3

sys.path.insert(0, str(ICE / "scripts" / "common"))


def f32(x: Any) -> np.ndarray:
    return np.asarray(x, dtype=np.float32)


class PyContext:
    """Setup of the clean full-sample reconstruction (cwd is changed to the example directory)."""

    def __init__(self) -> None:
        os.chdir(EX)
        import torch

        torch.set_num_threads(1)
        from icenine import reconstructor as RC
        from icenine.config_file import ConfigFile
        from icenine.experimental_data import ExperimentalData

        self.RC = RC
        cfg = ConfigFile.from_file(str(CONFIG))
        exp = ExperimentalData.from_binary_memmap(STACK, 180, 2)
        self.setup = RC.setup_reconstruction(cfg, exp_data=exp)
        self.rec = RC.AdaptiveVoxelReconstructor(self.setup)
        self.eta_limit = self.setup.exp_setup.get_eta_limit()
        self.mic = self.setup.sample.get_mic()

    def verts(self, voxel: int) -> Any:
        return self.RC._get_voxel_vertices(self.mic.voxels[voxel])

    def truth(self, voxel: int) -> np.ndarray:
        from icenine.mic_file import MicFile

        if not hasattr(self, "_truth"):
            self._truth = MicFile.read(str(TRUTH_MIC))
        return np.asarray(self._truth.voxels[voxel].orientation, dtype=np.float64)

    def cost_fn(self, arm: str, cpp_edge_rules: bool = False) -> Any:
        """arm 'screen' (radius 3, centre, Q-max 5, C++ rule) or 'full' (radius 0, all
        reflections)."""
        from icenine.cost_functions import VoxelCostFunction

        s = self.setup
        kw: Dict[str, Any] = dict(
            simulator=s.simulator,
            detector_list=s.detector_list,
            range_map=s.range_map,
            exp_data=s.exp_data,
            sample=s.sample,
            structure_list=s.structure_list,
            mode="hard",
            eta_limit=self.eta_limit,
            cpp_overlap_rule=True,
            center_screen=True,
            cpp_edge_rules=cpp_edge_rules,
        )
        if arm == "screen":
            return VoxelCostFunction(pixel_radius=SCREEN_RADIUS, max_q=SCREEN_QMAX, **kw)
        return VoxelCostFunction(pixel_radius=0, **kw)

    def local_grid(self) -> np.ndarray:
        from icenine.sampling import generate_local_grid_multi_level

        return generate_local_grid_multi_level(self.setup.search_params.local_grid_radius, 0, 1)


def probe_eval(
    ctx: PyContext, cost_fn: Any, R: np.ndarray, verts: Any, phase: int
) -> Dict[str, Any]:
    """Evaluate ``R`` (3x3, cast to float32 inside) and return the aggregate counts plus one
    record per peak (including peaks the frame lookup drops, flagged ``nomatch``)."""
    from icenine import cost_functions as CF

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
        info = cost_fn.evaluate(np.asarray(R, dtype=np.float32), verts, phase)
    finally:
        CF.get_scattering_omegas_torch, CF.calculate_diffraction_overlap_batched = (
            o_bragg,
            o_batched,
        )
        CF._c_stage_d_overlap = o_stage
    out: Dict[str, Any] = dict(
        pixov=int(info.pixel_overlap),
        pixon=int(info.pixel_on_detector),
        pkov=int(info.peak_overlap),
        pkon=int(info.peak_on_detector),
        quality=float(info.quality),
        conf=float(info.confidence),
        peaks=[],
    )
    if "kw" not in cap:  # no observable peak at all
        return out
    kw = cap["kw"]
    omegas = kw["peak_omegas"].numpy()
    normals = kw["peak_normals"].numpy()
    br = cap["bragg"]
    obs = br.observable.numpy().astype(bool)
    k_obs = np.nonzero(obs)[0]
    K = len(k_obs)
    rows_omega = np.concatenate([br.omega1.numpy()[obs], br.omega2.numpy()[obs]])
    # the evaluate() eta filter keeps a subsequence of rows: match it in order
    row_of = np.full(len(omegas), -1, dtype=np.int64)
    j = 0
    for p in range(len(omegas)):
        while j < len(rows_omega) and rows_omega[j] != omegas[p]:
            j += 1
        row_of[p] = j
        j += 1
    rng = kw["range_map"]
    wedge = np.full(len(omegas), -1, dtype=np.int64)
    bins_f = (omegas - rng.low) / rng.width
    bins = (np.floor(bins_f) if getattr(cost_fn, "cpp_edge_rules", False) else bins_f).astype(int)
    for p, b in enumerate(bins):
        if 0 <= b < len(rng.index_list) and rng.index_list[b] is not None:
            wedge[p] = rng.index_list[b]
    valid = np.nonzero(wedge >= 0)[0]
    st = cap.get("stage")
    per_peak_res: Dict[int, Tuple[int, ...]] = {}
    pix: Dict[int, Any] = {}
    if st is not None:
        # args: image_list, local_wedge, det_hit, centers, verts, n_det, radius, M, rows, cols,
        # [rule]
        image_list, local_wedge, det_hit, centers, vverts = st[0], st[1], st[2], st[3], st[4]
        tail = st[5:]
        for m, p in enumerate(valid):
            r = o_stage(
                image_list,
                local_wedge[m : m + 1],
                det_hit[m : m + 1],
                None if centers is None else np.ascontiguousarray(centers[m : m + 1]),
                None if vverts is None else np.ascontiguousarray(vverts[m : m + 1]),
                tail[0],
                tail[1],
                1,
                *tail[3:],
            )
            per_peak_res[int(p)] = tuple(int(x) for x in r[:3]) + (int(r[3]), int(r[6]))
            pix[int(p)] = (
                det_hit[m].astype(bool).tolist(),
                None if centers is None else centers[m].tolist(),
                None if vverts is None else vverts[m].tolist(),
            )
    m_of = {int(p): m for m, p in enumerate(valid)}
    for p in range(len(omegas)):
        row = int(row_of[p])
        refl = int(k_obs[row % K])
        rec: Dict[str, Any] = dict(
            refl=refl,
            branch=int(row // K),
            omega=float(omegas[p]),
            dir=normals[p].tolist(),
            mag=float(cost_fn._phase_recip_vecs[phase][1][refl]),
            frame=int(wedge[p]),
            nomatch=bool(wedge[p] < 0),
        )
        if p in per_peak_res:
            pr = per_peak_res[p]
            rec.update(pixov=pr[0], pixon=pr[1], pkov=pr[2], pkon=pr[3], ndet=pr[4])
            rec["hit"], rec["centers"], rec["verts"] = pix[p]
        out["peaks"].append(rec)
    return out


def python_pixels(ctx: PyContext, R: np.ndarray, verts: Any, refl: int, branch: int) -> np.ndarray:
    """Float32 pixel coordinates [det][vertex] -> (col, row) of one peak (stage B/C of the
    batched overlap, repeated for diagnosis of mismatching pixel boundaries only)."""
    import torch
    from icenine.diffraction_core import get_scattering_omegas_torch

    s = ctx.setup
    cf = ctx.cost_fn("full")
    g_hkl, g_mag = next(iter(cf._phase_recip_vecs.values()))
    g = (torch.from_numpy(np.asarray(R, dtype=np.float32)).float() @ g_hkl.T).T
    br = get_scattering_omegas_torch(
        g, g_mag, s.simulator.beam_energy, s.simulator.beam_deflection_chi
    )
    om = (br.omega1 if branch == 0 else br.omega2)[refl]
    n = g[refl] / torch.norm(g[refl])
    c, sn = torch.cos(om), torch.sin(om)
    Rz = torch.zeros(3, 3)
    Rz[0, 0], Rz[0, 1], Rz[1, 0], Rz[1, 1], Rz[2, 2] = c, -sn, sn, c, 1.0
    base = s.sample.sample_to_lab_matrix
    full = Rz @ base[:3, :3]
    ln = full @ n
    bd = s.simulator.beam_direction
    refd = bd - 2.0 * (bd * ln).sum() * ln
    v4 = torch.cat([verts, torch.ones(3, 1)], dim=1)
    f4 = torch.zeros(4, 4)
    f4[:3, :3] = full
    f4[:3, 3] = base[:3, 3]
    f4[3, 3] = 1.0
    lab = torch.einsum("ij,vj->vi", f4, v4)[:, :3]
    res = []
    for det in s.detector_list:
        pl = det._detector_plane
        t = -((lab * pl.normal).sum(1) + pl.d) / (refd * pl.normal).sum()
        pts = lab + t.unsqueeze(1) * refd
        rel = pts - det._position
        pl_ = rel - det._lab_frame_coord_origin
        j = (pl_ * det._lab_frame_basis_j).sum(1)
        k = (pl_ * det._lab_frame_basis_k).sum(1)
        col = (j + det.pixel_half_width) / det.pixel_width
        row = (k + det.pixel_half_height) / det.pixel_height
        res.append(torch.stack([col, row], dim=1).numpy())
    return np.stack(res)  # (n_det, 3, 2)


# --------------------------------------------------------------------------------------------
# C++ output parser
# --------------------------------------------------------------------------------------------
def hx(tok: str) -> np.float32:
    return np.float32(float.fromhex(tok))


def parse_cpp(path: Path) -> Dict[str, Any]:
    """Records ``{id: {kind, prod, arms: {INT|EXS|EXF: {...}}}}`` plus header and voxel lines."""
    head: Dict[str, Any] = {}
    voxel: Dict[str, List[float]] = {}
    recs: Dict[int, Dict[str, Any]] = {}
    cur: Optional[Dict[str, Any]] = None
    arm: Optional[Dict[str, Any]] = None
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if not t:
                continue
            if t[0] == "H":
                head = {t[i]: int(t[i + 1]) for i in range(1, len(t), 2)}
            elif t[0] == "V":
                voxel[t[1]] = [float.fromhex(x) for x in t[2:5]]
            elif t[0] == "O":
                cur = dict(
                    id=int(t[1]),
                    kind=int(t[3]),
                    prod=np.array([hx(x) for x in t[5:14]], np.float32).reshape(3, 3),
                    arms={},
                )
                recs[cur["id"]] = cur
            elif t[0] == "A":
                arm = dict(peaks=[], obsrefl=None)
                assert cur is not None
                cur["arms"][t[1]] = arm
            elif t[0] == "S":
                assert arm is not None
                arm["nobs"] = int(t[2])
                arm["obsrefl"] = [int(x) for x in t[4:]]
            elif t[0] == "R":
                assert arm is not None
                arm.update(
                    pixov=int(t[2]),
                    pixon=int(t[4]),
                    pkov=int(t[6]),
                    pkon=int(t[8]),
                    quality=float(hx(t[10])),
                    conf=float(hx(t[12])),
                )
            elif t[0] == "P":
                assert arm is not None
                arm["peaks"].append(parse_peak(line))
    return dict(head=head, voxel=voxel, recs=recs)


def parse_peak(line: str) -> Dict[str, Any]:
    parts = [p.split() for p in line.split("|")]
    a = parts[0]
    d: Dict[str, Any] = dict(
        idx=int(a[1]), refl_code=int(a[3]), omega=float(hx(a[5])), frame=int(a[7])
    )
    d["refl"], d["branch"] = (
        (d["refl_code"] >> 1, d["refl_code"] & 1) if d["refl_code"] >= 0 else (-1, -1)
    )
    if len(a) > 8 and a[8] == "NOMATCH":
        d["nomatch"] = True
        d["dir"] = [float(hx(a[10])), float(hx(a[11])), float(hx(a[12]))]
        d["mag"] = float(hx(a[14]))
        return d
    d["nomatch"] = False
    d["dir"] = [float(hx(a[9])), float(hx(a[10])), float(hx(a[11]))]
    d["mag"] = float(hx(a[13]))
    d["dets"] = []
    for p in parts[1:-1]:
        h = int(p[1][1:])
        xy = [hx(x) for x in p[2 : 2 + 2 * h]]
        d["dets"].append(
            dict(
                hit=h,
                xy=np.array(xy, np.float32).reshape(-1, 2),
                lit=int(p[2 + 2 * h][3:]),
                ov=int(p[3 + 2 * h][2:]),
            )
        )
    z = parts[-1]
    d.update(pixov=int(z[1]), pixon=int(z[3]), pkov=int(z[5]), pkon=int(z[7]), ndet=int(z[9]))
    return d


def pix64(
    ctx: PyContext, cost_fn: Any, R: np.ndarray, verts: Any, refl: int, branch: int, phase: int
) -> np.ndarray:
    """Double-precision pixel coordinates (n_det, 3, 2) = (col, row) of one peak, from the same
    float32 inputs (orientation, reflection list, vertices, detector constants). The reference for
    deciding which side of a pixel boundary the true position lies on."""
    s = ctx.setup
    g_hkl = cost_fn._phase_recip_vecs[phase][0].numpy().astype(np.float64)
    g_mag = float(cost_fn._phase_recip_vecs[phase][1][refl])
    G = np.asarray(R, dtype=np.float64) @ g_hkl[refl]
    E = float(s.simulator.beam_energy)
    chi = float(s.simulator.beam_deflection_chi)
    k = 0.506773182 * E
    sin_t = g_mag / (2.0 * k)
    cos_c = G[2] / g_mag
    sin_c = np.sqrt(1.0 - cos_c * cos_c)
    num = sin_t + cos_c * np.sin(chi)
    den = sin_c * np.cos(chi)
    d0 = np.arctan2(G[0], G[1])
    b1 = np.arcsin(num / den)
    om = (b1 if branch == 0 else np.pi - b1) + d0
    om = (om + np.pi) % (2 * np.pi) - np.pi
    n = G / np.linalg.norm(G)
    Rz = np.array([[np.cos(om), -np.sin(om), 0], [np.sin(om), np.cos(om), 0], [0, 0, 1.0]])
    base = s.sample.sample_to_lab_matrix.numpy().astype(np.float64)
    full = Rz @ base[:3, :3]
    ln = full @ n
    bd = s.simulator.beam_direction.numpy().astype(np.float64)
    ref = bd - 2.0 * (bd @ ln) * ln
    lab = (verts.numpy().astype(np.float64) @ full.T) + base[:3, 3]
    out = []
    for det in s.detector_list:
        pl = det._detector_plane
        nrm = pl.normal.numpy().astype(np.float64)
        t = -((lab @ nrm) + float(pl.d)) / float(ref @ nrm)
        pts = lab + t[:, None] * ref
        rel = pts - det._position.numpy().astype(np.float64)
        pl_ = rel - det._lab_frame_coord_origin.numpy().astype(np.float64)
        j = pl_ @ det._lab_frame_basis_j.numpy().astype(np.float64)
        kk = pl_ @ det._lab_frame_basis_k.numpy().astype(np.float64)
        col = (j + float(det.pixel_half_width)) / float(det.pixel_width)
        row = (kk + float(det.pixel_half_height)) / float(det.pixel_height)
        out.append(np.stack([col, row], axis=1))
    return np.stack(out)
