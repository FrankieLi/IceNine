"""Shared set-up and kernels of the S4 study (translation correctness and loop ordering).

Plan: ``docs/todo_translation_screen_prototype.md``. Scripts only, no library change.

* ``load_set("phase_d" | "three")`` builds the reconstruction setup exactly as the reconstructor
  is built for that data set (Phase D clean: memmap stack of ``full/clean``, the
  ``ReconstructPhaseD_cma_clean`` config; ThreeVoxels: ``ReconstructQ8.config`` with the
  ``ScatteringData_Python`` images).
* ``capture_level0`` runs ``AdaptiveVoxelReconstructor.reconstruct_voxel`` up to the call of
  ``run_discrete_search_spaced`` at level 0 and returns the arguments it was given (global and
  local cost functions, FZ set, local grid, screen kwargs), so no setting is copied by hand.
* ``screen_reference`` is the level-0 screen exactly as ``run_discrete_search_spaced`` runs it
  with ``BatchedScreen`` (chunks of ``screen_chunk`` FZ points times the local grid).
* ``PeakTable`` holds, per orientation and peak, the frame and the float64 affine map
  ``u = A x + c`` of every detector (the peak set comes from the float32 torch operations of
  ``ScreenEvaluator._counts``, unchanged, so it matches the reference).
* ``dilate`` builds the padded, dilated binary frame that turns the pixel-radius test into one
  lookup; ``kernel_numpy`` is the reference implementation of the table-driven screen.
"""

import contextlib
import io
import math
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, Tuple

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[2]
ROOT = ICENINE_PY.parent
sys.path.insert(0, str(ICENINE_PY))
sys.path.insert(0, str(ICENINE_PY / "scripts" / "common"))

from icenine import cost_functions as CF  # noqa: E402
from icenine import reconstructor as REC  # noqa: E402
from icenine.diffraction_core import get_scattering_omegas_torch  # noqa: E402

EX_D = ROOT / "Examples" / "Example2.ManyGrains"
EX_3 = ROOT / "Examples" / "Example2.ThreeVoxels"
CFG_D = ICENINE_PY / "scripts" / "phase_d" / "configs" / "ReconstructPhaseD_cma_clean.config"
STACK_D = EX_D / "ScatteringData_PhaseD" / "stacks" / "clean.npy"
MIC_D = EX_D / "SimInput" / "rand_500grains_1mm_neworient_s0_grid.mic"
TRUTH_D = EX_D / "SimInput" / "rand_500grains_1mm_neworient_s0.mic"
OUT = ICENINE_PY / "benchmarks" / "s4_translation"
CACHE = OUT / "cache"
REF_SEED = 20261014  # seed of the reference voxel draw


# ---------------------------------------------------------------------------------------------
# data sets
# ---------------------------------------------------------------------------------------------


@dataclass
class DataSet:
    name: str
    rec: Any  # AdaptiveVoxelReconstructor
    setup: Any  # ReconstructionSetup
    voxels: List[Any]  # sample voxels (positions; orientations not used)
    truth: Optional[np.ndarray]  # (V, 3, 3) true orientations, or None
    example_dir: Path


def load_set(name: str, threads: int = 1) -> DataSet:
    """Reconstruction setup of a data set; cwd is the example directory afterwards."""
    torch.set_num_threads(threads)
    from icenine.config_file import ConfigFile
    from icenine.experimental_data import ExperimentalData
    from icenine.mic_file import MicFile
    from icenine.reconstructor import AdaptiveVoxelReconstructor, setup_reconstruction

    if name == "phase_d":
        os.chdir(EX_D)
        cfg = ConfigFile.from_file(str(CFG_D))
        exp = ExperimentalData.from_binary_memmap(STACK_D, 180, 2)
        with contextlib.redirect_stdout(io.StringIO()):
            setup = setup_reconstruction(cfg, exp_data=exp)
        voxels = list(setup.sample.get_mic().voxels)
        truth_mic = MicFile.read(str(TRUTH_D))
        truth = np.stack([np.asarray(v.orientation, np.float64) for v in truth_mic.voxels])
        assert len(truth) == len(voxels)
        return DataSet(name, AdaptiveVoxelReconstructor(setup), setup, voxels, truth, EX_D)
    if name == "three":
        from icenine.sampling import load_fundamental_zone_file

        os.chdir(EX_3)
        cfg = ConfigFile.from_file(str(EX_3 / "ConfigFiles" / "ReconstructQ8.config"))
        exp = ExperimentalData.from_image_directory(
            directory=str(EX_3 / "ScatteringData_Python"),
            basename="3Grains.sim",
            ext="d",
            serial_length=5,
            n_omega=180,
            n_detectors=2,
            num_rows=2048,
            num_cols=2048,
            mode="sparse",
        )
        exp.prepare_for_reconstruction()
        fz = load_fundamental_zone_file(cfg.fundamental_zone_filename)
        with contextlib.redirect_stdout(io.StringIO()):
            setup = setup_reconstruction(cfg, exp_data=exp, fz_orientations=fz)
        mic = MicFile.read(str(EX_3 / "SimInput" / "three_voxels.mic"))
        voxels = list(mic.voxels)
        truth = np.stack([np.asarray(v.orientation, np.float64) for v in voxels])
        return DataSet(name, AdaptiveVoxelReconstructor(setup), setup, voxels, truth, EX_3)
    raise ValueError(name)


def voxel_vertices(voxel: Any) -> torch.Tensor:
    return REC._get_voxel_vertices(voxel)


def centroid64(verts: torch.Tensor) -> np.ndarray:
    """Centroid as the centre screen computes it (float32 mean), widened to float64."""
    return verts.mean(dim=0).double().numpy()


def reference_voxels(ds: DataSet, n: int = 20) -> List[int]:
    """Seeded random Phase D voxels, half pointing up and half down (both triangle types)."""
    rng = np.random.default_rng(REF_SEED)
    up = np.array([bool(v.points_up) for v in ds.voxels])
    a = rng.choice(np.nonzero(up)[0], n // 2, replace=False)
    b = rng.choice(np.nonzero(~up)[0], n - n // 2, replace=False)
    return sorted(int(i) for i in np.concatenate([a, b]))


def timed_voxels(ds: DataSet, n: int = 10) -> List[int]:
    rng = np.random.default_rng(REF_SEED + 1)
    return sorted(int(i) for i in rng.choice(len(ds.voxels), n, replace=False))


# ---------------------------------------------------------------------------------------------
# level-0 capture
# ---------------------------------------------------------------------------------------------


class _Captured(Exception):
    pass


def capture_level0(ds: DataSet, voxel: Any) -> Dict[str, Any]:
    """The arguments ``reconstruct_voxel`` passes to ``run_discrete_search_spaced`` at level 0."""
    box: Dict[str, Any] = {}
    orig = REC.run_discrete_search_spaced

    def fake(**kw: Any) -> Any:
        box.update(kw)
        raise _Captured()

    REC.run_discrete_search_spaced = fake  # type: ignore[assignment]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            ds.rec.reconstruct_voxel(
                voxel_vertices=voxel_vertices(voxel),
                phase_index=voxel.phase,
                rng=np.random.default_rng(0),
            )
    except _Captured:
        pass
    finally:
        REC.run_discrete_search_spaced = orig  # type: ignore[assignment]
    if not box:
        raise RuntimeError("level-0 discrete search was not reached")
    return box


def chunk_orientations(fz: np.ndarray, local_grid: np.ndarray, chunk: int) -> Iterator[np.ndarray]:
    """(chunk * n_local, 3, 3) orientation blocks exactly as run_discrete_search_spaced builds."""
    for s in range(0, len(fz), chunk):
        c = np.asarray(fz[s : s + chunk])
        prod = np.matmul(local_grid[None, :, :, :], c[:, None, :, :])
        yield prod.reshape(-1, 3, 3)


def all_orientations(cap: Dict[str, Any]) -> np.ndarray:
    chunk = int(cap.get("screen_chunk", 32))
    return np.concatenate(
        list(chunk_orientations(cap["fz_orientations"], cap["local_grid"], chunk))
    )


def screen_reference(
    cap: Dict[str, Any], verts: torch.Tensor, phase: int
) -> Tuple[np.ndarray, np.ndarray]:
    """(peak_overlap, peak_on_detector), each (n_fz * n_local,), of the BatchedScreen screen."""
    screen = cap["screen"]
    assert not cap.get("screen_interp", False)
    chunk = int(cap.get("screen_chunk", 32))
    po, pon = [], []
    for block in chunk_orientations(cap["fz_orientations"], cap["local_grid"], chunk):
        a, b = screen.overlaps_exact(block, verts, phase)
        po.append(a)
        pon.append(b)
    return np.concatenate(po), np.concatenate(pon)


# ---------------------------------------------------------------------------------------------
# peaks of an orientation block (float32, as ScreenEvaluator) and their float64 affine maps
# ---------------------------------------------------------------------------------------------


@dataclass
class DetGeom:
    """Detector geometry in float64 (from the float32 tensors the cost function uses)."""

    n: np.ndarray
    d: float
    pos_org: np.ndarray  # _position + _lab_frame_coord_origin
    bj: np.ndarray
    bk: np.ndarray
    hw: float
    hh: float
    pw: float
    ph: float
    rows: int
    cols: int


def det_geoms(cf: Any) -> List[DetGeom]:
    out = []
    for det in cf.detector_list:
        pl = det._detector_plane
        f64 = lambda t: np.asarray(
            t.detach().double().numpy() if torch.is_tensor(t) else t, float
        )  # noqa: E731
        out.append(
            DetGeom(
                n=f64(pl.normal),
                d=float(pl.d),
                pos_org=f64(det._position) + f64(det._lab_frame_coord_origin),
                bj=f64(det._lab_frame_basis_j),
                bk=f64(det._lab_frame_basis_k),
                hw=float(det.pixel_half_width),
                hh=float(det.pixel_half_height),
                pw=float(det.pixel_width),
                ph=float(det.pixel_height),
                rows=int(det.num_rows) if hasattr(det, "num_rows") else 2048,
                cols=int(det.num_cols) if hasattr(det, "num_cols") else 2048,
            )
        )
    return out


@dataclass
class PeakTable:
    """Peaks of a set of orientations, sorted by orientation (segment).

    seg (P,) orientation index; frame (P,) global frame index; A (P, n_det, 2, 3) and
    c (P, n_det, 2) with u = A x + c (col, row) in float64; hit (P, n_det) the plane-hit flag
    (non-parallel and ray parameter > 0, checked constant over the sample); la, lb: the ray
    parameter lambda(x) = la . x + lb (P, n_det, 3) and (P, n_det) for the check; full_rot and
    refl (float32 inputs, widened) for re-use.
    """

    seg: np.ndarray
    frame: np.ndarray
    A: np.ndarray
    c: np.ndarray
    hit: np.ndarray
    la: np.ndarray
    lb: np.ndarray
    full_rot: np.ndarray
    refl: np.ndarray
    n_orient: int
    meta: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_peaks(self) -> int:
        return int(len(self.seg))


def affine_maps(
    full_rot: np.ndarray, refl: np.ndarray, t: np.ndarray, geoms: List[DetGeom]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Float64 A, c, parallel flag, lambda coefficients of peaks (full_rot (P,3,3), refl (P,3)).
    (macOS Accelerate matmul raises spurious floating-point warnings; inputs and outputs are
    asserted finite instead.)"""
    assert np.isfinite(full_rot).all() and np.isfinite(refl).all()
    with np.errstate(all="ignore"):
        out = _affine_maps(full_rot, refl, t, geoms)
    assert all(np.isfinite(a).all() for a in (out[0], out[1], out[3], out[4]))
    return out


def _affine_maps(
    full_rot: np.ndarray, refl: np.ndarray, t: np.ndarray, geoms: List[DetGeom]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    P = len(full_rot)
    nd = len(geoms)
    A = np.zeros((P, nd, 2, 3))
    c = np.zeros((P, nd, 2))
    par = np.zeros((P, nd), bool)
    la = np.zeros((P, nd, 3))
    lb = np.zeros((P, nd))
    for di, g in enumerate(geoms):
        denom = refl @ g.n  # (P,)
        par[:, di] = np.abs(denom) < 1e-8
        sd = np.where(par[:, di], 1.0, denom)
        # Q = I - k n^T / denom ; lab = M x + t ; pts = Q lab - k D / denom
        nM = np.einsum("i,pij->pj", g.n, full_rot)  # n^T M (P, 3)
        for ax, (b, half, size) in enumerate(((g.bj, g.hw, g.pw), (g.bk, g.hh, g.ph))):
            bk_ = refl @ b  # b . k (P,)
            bM = np.einsum("i,pij->pj", b, full_rot)  # b^T M
            # b^T Q M = b^T M - (b.k) (n^T M) / denom
            A[:, di, ax, :] = (bM - (bk_ / sd)[:, None] * nM) / size
            bQt = b @ t - bk_ * (g.n @ t) / sd
            c[:, di, ax] = (bQt - bk_ * g.d / sd - b @ g.pos_org + half) / size
        la[:, di, :] = -nM / sd[:, None]
        lb[:, di] = -(g.n @ t + g.d) / sd
    return A, c, par, la, lb


def project_direct64(
    full_rot: np.ndarray, refl: np.ndarray, t: np.ndarray, geoms: List[DetGeom], X: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Step-by-step float64 stages B and C at points X (Q, 3): (u (P, n_det, Q, 2), lam)."""
    with np.errstate(all="ignore"):
        return _project_direct64(full_rot, refl, t, geoms, X)


def _project_direct64(
    full_rot: np.ndarray, refl: np.ndarray, t: np.ndarray, geoms: List[DetGeom], X: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    lab = np.einsum("pij,qj->pqi", full_rot, X) + t  # (P, Q, 3)
    out = np.zeros((len(full_rot), len(geoms), len(X), 2))
    lam_out = np.zeros((len(full_rot), len(geoms), len(X)))
    for di, g in enumerate(geoms):
        denom = refl @ g.n
        numer = -(lab @ g.n + g.d)
        lam = numer / np.where(np.abs(denom) < 1e-8, 1.0, denom)[:, None]
        pts = lab + lam[:, :, None] * refl[:, None, :]
        loc = pts - g.pos_org
        out[:, di, :, 0] = (loc @ g.bj + g.hw) / g.pw
        out[:, di, :, 1] = (loc @ g.bk + g.hh) / g.ph
        lam_out[:, di] = lam
    return out, lam_out


def block_peaks_f32(cf: Any, orientations: np.ndarray, phase: int) -> Dict[str, Any]:
    """The float32 peaks of ScreenEvaluator.overlaps_exact/_counts up to the frame lookup, with
    the same operations (so the same bits for the same block shape): seg, frame, full_rot_v,
    reflected_v, the eta filter, and the stage-B/C float32 centre pixels (for the flip report).
    """
    K = len(orientations)
    g_hkl, g_mag = cf._phase_recip_vecs[phase]
    n = g_hkl.shape[0]
    R = torch.from_numpy(np.asarray(orientations, dtype=np.float64)).float()
    g = torch.matmul(R, g_hkl.T).transpose(1, 2).reshape(K * n, 3)
    br = get_scattering_omegas_torch(
        g, g_mag.repeat(K), cf.simulator.beam_energy, cf.simulator.beam_deflection_chi
    )
    obs = br.observable
    og = g[obs]
    normals = og / torch.norm(og, dim=1, keepdim=True)
    seg_one = torch.arange(K).repeat_interleave(n)[obs]
    omegas = torch.cat([br.omega1[obs], br.omega2[obs]])
    normals = torch.cat([normals, normals])
    seg = torch.cat([seg_one, seg_one])
    sim, sample = cf.simulator, cf.sample
    beam_dir = sim.beam_direction
    base_rot = sample.sample_to_lab_matrix[:3, :3]
    cos_w, sin_w = torch.cos(omegas), torch.sin(omegas)
    N = omegas.shape[0]
    Rz = torch.zeros(N, 3, 3)
    Rz[:, 0, 0], Rz[:, 0, 1], Rz[:, 1, 0], Rz[:, 1, 1], Rz[:, 2, 2] = (
        cos_w,
        -sin_w,
        sin_w,
        cos_w,
        1.0,
    )
    full_rot = Rz @ base_rot
    lab_normals = torch.bmm(full_rot, normals.unsqueeze(-1)).squeeze(-1)
    dot = (beam_dir.unsqueeze(0) * lab_normals).sum(dim=1, keepdim=True)
    reflected = beam_dir.unsqueeze(0) - 2.0 * dot * lab_normals
    rd_norms = torch.norm(reflected, dim=1)
    safe = torch.where(rd_norms > 0, rd_norms, torch.ones_like(rd_norms))
    eta = torch.atan2(torch.abs(reflected[:, 1]) / safe, torch.abs(reflected[:, 2]) / safe)
    valid = (eta < cf.eta_limit) & (rd_norms > 0)
    if cf.min_sin_eta > 0.0:
        valid = valid & (torch.sin(eta) >= cf.min_sin_eta)
    rng = cf.range_map
    index = np.array([-1 if i is None else int(i) for i in rng.index_list], dtype=np.int64)
    nb = len(index)
    bins = CF._omega_bins(omegas.numpy(), rng.low, rng.width, cf.cpp_edge_rules)
    ok_bin = (bins >= 0) & (bins < nb)
    wedge = np.where(ok_bin, index[np.clip(bins, 0, nb - 1)], -1)
    keep = np.nonzero(valid.numpy() & (wedge >= 0))[0]
    order = keep[np.argsort(seg.numpy()[keep], kind="stable")]
    o_idx = torch.from_numpy(order)
    return dict(
        seg=seg.numpy()[order].astype(np.int64),
        frame=wedge[order].astype(np.int64),
        full_rot=full_rot[o_idx],
        reflected=reflected[o_idx],
        omega=omegas[o_idx],
        normal=normals[o_idx],
        n_orient=K,
        n_before_filter=int(N),
    )


def centre_pixels_f32(cf: Any, pk: Dict[str, Any], verts: torch.Tensor) -> Tuple[np.ndarray, ...]:
    """Stage B/C of ScreenEvaluator._counts (float32), centre screen: (centres (M, n_det, 2)
    int32, det_hit (M, n_det) uint8, continuous float32 coordinates (M, n_det, 2))."""
    sample = cf.sample
    full_rot_v, reflected_v = pk["full_rot"], pk["reflected"]
    M = full_rot_v.shape[0]
    orig_matrix = sample.sample_to_lab_matrix
    full_4x4 = torch.zeros(M, 4, 4)
    full_4x4[:, :3, :3] = full_rot_v
    full_4x4[:, :3, 3] = orig_matrix[:3, 3]
    full_4x4[:, 3, 3] = 1.0
    if cf.center_screen:
        verts = verts.mean(dim=0, keepdim=True).repeat(3, 1)
    verts_4d = torch.cat([verts, torch.ones(3, 1)], dim=1)
    lab_verts = torch.einsum("mij,vj->mvi", full_4x4, verts_4d)[:, :, :3]
    det_hit, centers, conts = [], [], []
    for det in cf.detector_list:
        plane = det._detector_plane
        origins = lab_verts.reshape(M * 3, 3)
        dirs = reflected_v.unsqueeze(1).expand(M, 3, 3).reshape(M * 3, 3)
        denom = (dirs * plane.normal).sum(dim=1)
        numer = -((origins * plane.normal).sum(dim=1) + plane.d)
        parallel = torch.abs(denom) < 1e-8
        safe_d = torch.where(parallel, torch.ones_like(denom), denom)
        t = torch.where(parallel, torch.zeros_like(denom), numer / safe_d)
        hits = (~parallel) & (t > 0)
        pts = origins + t.unsqueeze(1) * dirs
        relative = pts - det._position
        pixel_loc = relative - det._lab_frame_coord_origin
        j = (pixel_loc * det._lab_frame_basis_j).sum(dim=1)
        k = (pixel_loc * det._lab_frame_basis_k).sum(dim=1)
        col = CF._pixel_coordinate(j, det.pixel_half_width, det.pixel_width, cf.cpp_edge_rules)
        row = CF._pixel_coordinate(k, det.pixel_half_height, det.pixel_height, cf.cpp_edge_rules)
        det_hit.append(hits.reshape(M, 3).all(dim=1).numpy().astype(np.uint8))
        pix = torch.stack([col.reshape(M, 3), row.reshape(M, 3)], dim=2)
        centers.append(pix[:, 0, :].numpy().astype(np.int32))
        conts.append(pix[:, 0, :].numpy())
    return (
        np.ascontiguousarray(np.stack(centers, axis=1)),
        np.ascontiguousarray(np.stack(det_hit, axis=1)),
        np.stack(conts, axis=1),
    )


def sample_translation(cf: Any) -> np.ndarray:
    return cf.sample.sample_to_lab_matrix[:3, 3].double().numpy()


def build_table(
    cf: Any,
    orientations_blocks: List[np.ndarray],
    phase: int,
    extent: Optional[np.ndarray] = None,
) -> PeakTable:
    """Peak table of several orientation blocks (each block runs the float32 ops once, so the
    block shape must equal the reference's for bit-identical peaks). ``extent`` (Q, 3) points
    (e.g. the corners of the sample's bounding box) at which the hit flag is checked constant."""
    geoms = det_geoms(cf)
    t = sample_translation(cf)
    segs, frames, As, cs, hits, las, lbs, frs, rfs = [], [], [], [], [], [], [], [], []
    off = 0
    for blk in orientations_blocks:
        pk = block_peaks_f32(cf, blk, phase)
        fr = pk["full_rot"].double().numpy()
        rf = pk["reflected"].double().numpy()
        A, c, par, la, lb = affine_maps(fr, rf, t, geoms)
        if extent is not None:
            lam = np.einsum("pdj,qj->pdq", la, extent) + lb[:, :, None]
            pos = lam > 0
            assert (pos.all(axis=2) | (~pos).all(axis=2)).all(), "ray sign varies over sample"
            hit = (~par) & pos[:, :, 0]
        else:
            hit = ~par
        segs.append(pk["seg"] + off)
        frames.append(pk["frame"])
        As.append(A)
        cs.append(c)
        hits.append(hit)
        las.append(la)
        lbs.append(lb)
        frs.append(fr)
        rfs.append(rf)
        off += len(blk)
    return PeakTable(
        seg=np.concatenate(segs),
        frame=np.concatenate(frames),
        A=np.concatenate(As),
        c=np.concatenate(cs),
        hit=np.concatenate(hits),
        la=np.concatenate(las),
        lb=np.concatenate(lbs),
        full_rot=np.concatenate(frs),
        refl=np.concatenate(rfs),
        n_orient=off,
    )


def sample_extent(voxels: List[Any]) -> np.ndarray:
    """Corners of the bounding box of every voxel vertex (the ray parameter is affine in x)."""
    P = np.array([np.asarray(v.position, float) for v in voxels])
    s = max(float(v.side_length) for v in voxels)
    lo, hi = P.min(axis=0) - s, P.max(axis=0) + s
    return np.array(
        [[x, y, z] for x in (lo[0], hi[0]) for y in (lo[1], hi[1]) for z in (lo[2], hi[2])]
    )


# ---------------------------------------------------------------------------------------------
# dilated frames and the numpy reference kernel
# ---------------------------------------------------------------------------------------------


def dilate(binary: np.ndarray, r: int) -> np.ndarray:
    """(rows + 2r, cols + 2r) uint8: entry [y + r, x + r] = any lit pixel of the image in the
    (2r+1)^2 box centred on pixel (x, y) (pixels outside the image count as unlit)."""
    from scipy.ndimage import maximum_filter1d

    rows, cols = binary.shape
    pad = np.zeros((rows + 4 * r, cols + 4 * r), np.uint8)
    pad[2 * r : 2 * r + rows, 2 * r : 2 * r + cols] = binary != 0
    m = maximum_filter1d(pad, 2 * r + 1, axis=0, mode="constant", cval=0)
    m = maximum_filter1d(m, 2 * r + 1, axis=1, mode="constant", cval=0)
    return np.ascontiguousarray(m[r : r + rows + 2 * r, r : r + cols + 2 * r])


def dilate_sparse(binary: np.ndarray, r: int) -> np.ndarray:
    """``dilate`` built by stamping the box around each lit pixel (fast for sparse images)."""
    rows, cols = binary.shape
    out = np.zeros((rows + 2 * r, cols + 2 * r), np.uint8)
    ys, xs = np.nonzero(binary)
    for dy in range(-r, r + 1):
        for dx in range(-r, r + 1):
            # lit pixel (x, y) makes centre (x - dx, y - dy) bright: index (y - dy + r, x - dx + r)
            out[ys - dy + r, xs - dx + r] = 1
    return out


def trunc_index(u: np.ndarray) -> np.ndarray:
    """Pixel index of a continuous coordinate as the screen takes it: float -> int32 cast, i.e.
    truncation toward zero (``.astype(np.int32)``); clipped far outside to stay in range."""
    return np.trunc(np.clip(u, -1e7, 1e7)).astype(np.int64)


def qualified(lit: np.ndarray, bright: np.ndarray, cpp_rule: bool) -> Tuple[np.ndarray, np.ndarray]:
    """count_qualified_peaks_c on (..., n_det) bool arrays: (peak_on_det, peak_overlap)."""
    nd = lit.shape[-1]
    if cpp_rule:
        valid_sim = lit[..., 0].copy()
        on_prev = lit[..., 0].copy()
        for i in range(1, nd):
            li = lit[..., i]
            new_valid = valid_sim & (li | on_prev)
            on_prev = np.where(valid_sim & ~li & on_prev, False, on_prev)
            valid_sim = new_valid
        prev_ov = bright[..., 0].copy()
        valid_ov = bright[..., 0].copy()
        for i in range(1, nd):
            bi, li = bright[..., i], lit[..., i]
            new_valid = valid_ov & (bi | (prev_ov & ~li))
            prev_ov = np.where(valid_ov & ~bi & prev_ov & ~li, False, prev_ov)
            valid_ov = new_valid
        return valid_sim, valid_ov
    valid_sim = lit[..., 0].copy()
    gap = np.zeros_like(valid_sim)
    for i in range(1, nd):
        valid_sim &= ~(gap & lit[..., i])
        gap |= ~lit[..., i]
    valid_ov = bright[..., 0].copy()
    gap = np.zeros_like(valid_ov)
    for i in range(1, nd):
        valid_ov &= ~(gap & lit[..., i] & bright[..., i])
        gap |= ~bright[..., i] & lit[..., i]
    return valid_sim, valid_sim & valid_ov


class Bitmaps:
    """Bit-packed dilated frames, [frame, det] -> (rows + 2r, ceil((cols + 2r) / 8)) uint8."""

    def __init__(self, packed: np.ndarray, r: int, rows: int, cols: int) -> None:
        self.packed = packed  # (n_frames, n_det, R2, C8)
        self.r = r
        self.rows = rows
        self.cols = cols

    def lookup(self, f: np.ndarray, d: int, cx: np.ndarray, cy: np.ndarray) -> np.ndarray:
        """(in_bounds, bright) of centre pixels (cx, cy) on frame f (arrays broadcast)."""
        r = self.r
        inb = (cx >= -r) & (cx < self.cols + r) & (cy >= -r) & (cy < self.rows + r)
        X = np.where(inb, cx + r, 0)
        Y = np.where(inb, cy + r, 0)
        byte = self.packed[f, d, Y, X >> 3]
        bit = (byte >> (7 - (X & 7))) & 1
        return inb, inb & (bit == 1)


def build_bitmaps(
    exp_data: Any,
    n_frames: int,
    n_det: int,
    r: int,
    frames: Optional[np.ndarray] = None,
    sparse: bool = False,
) -> Bitmaps:
    """Packed dilated bitmaps of the frames (``sparse``: built with ``dilate_sparse``)."""
    first = exp_data.get_image(0, 0).get_binary_numpy()
    rows, cols = first.shape
    C8 = (cols + 2 * r + 7) // 8
    packed = np.zeros((n_frames, n_det, rows + 2 * r, C8), np.uint8)
    use = range(n_frames) if frames is None else [int(f) for f in frames]
    for f in use:
        for d in range(n_det):
            img = exp_data.get_image(f, d).get_binary_numpy()
            dil = dilate_sparse(img, r) if sparse else dilate(img, r)
            packed[f, d] = np.packbits(dil, axis=1)
    return Bitmaps(packed, r, rows, cols)


def kernel_numpy(
    tab: PeakTable,
    X: np.ndarray,
    bm: Bitmaps,
    cpp_rule: bool,
    peak_slice: Optional[slice] = None,
    return_coords: bool = False,
) -> Dict[str, np.ndarray]:
    """Reference table-driven screen: peak_overlap and peak_on_det (n_orient, V) for voxel
    centroids X (V, 3) float64."""
    sl = peak_slice or slice(0, tab.n_peaks)
    A, c, hit, f, seg = tab.A[sl], tab.c[sl], tab.hit[sl], tab.frame[sl], tab.seg[sl]
    P, nd = hit.shape
    V = len(X)
    U = np.einsum("pdaj,vj->pdva", A, X) + c[:, :, None, :]  # (P, nd, V, 2)
    lit = np.zeros((P, V, nd), bool)
    bright = np.zeros((P, V, nd), bool)
    for d in range(nd):
        cx = trunc_index(U[:, d, :, 0])
        cy = trunc_index(U[:, d, :, 1])
        inb, br = bm.lookup(f[:, None], d, cx, cy)
        lit[:, :, d] = inb & hit[:, d][:, None]
        bright[:, :, d] = br & hit[:, d][:, None]
    on_det, ovl = qualified(lit, bright, cpp_rule)
    po = np.zeros((tab.n_orient, V), np.int64)
    pon = np.zeros((tab.n_orient, V), np.int64)
    np.add.at(po, seg, ovl.astype(np.int64))
    np.add.at(pon, seg, on_det.astype(np.int64))
    out = dict(po=po, pon=pon)
    if return_coords:
        out["U"] = U
    return out


def preflight_json(path: Path) -> Dict[str, Any]:
    """Preflight saved as JSON, with "quiet" from preflight.require_quiet (not raised: the machine
    has been contended, owner decision D-20261014-I16; the label says so)."""
    import preflight

    info = preflight.preflight()
    try:
        preflight.require_quiet(info=info)
        info["quiet"] = True
    except preflight.MachineBusyError as e:
        info["quiet"] = False
        info["not_quiet_reason"] = str(e)
    path.parent.mkdir(parents=True, exist_ok=True)
    import json

    path.write_text(json.dumps(info, indent=2) + "\n")
    return info


def label_of(info: Dict[str, Any]) -> str:
    return "quiet" if info.get("quiet") else "contended"


@contextlib.contextmanager
def quiet_stdout() -> Iterator[None]:
    with contextlib.redirect_stdout(io.StringIO()):
        yield


def now() -> float:
    return time.perf_counter()
