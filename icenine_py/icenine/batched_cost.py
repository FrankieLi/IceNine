"""Batched local cost: many orientations per call, equal to ``VoxelCostFunction.evaluate``.

``VoxelCostFunction.evaluate`` scores one orientation per call. Of its ~330 us at Q-max 8 only
about 10 us is physics and pixel work; the rest is per-call torch and Python overhead (S3,
MIGRATION_HISTORY "Python speed parity"). ``BatchedCostEvaluator`` scores B orientations of the
*local* cost (hard mode, ``pixel_radius`` 0, no ray rule) in one call, used by the lockstep quick
MC (``batched_quick_mc.py``, config key ``BatchedQuickMC``, opt-in).

Two forward variants, both reproducing ``evaluate`` (stages A-C in torch, stage D the existing C
call once per orientation):

* **compact** (CPU): rows of the observable reflections only, boolean-mask compaction after the eta
  filter and after the frame lookup, as ``ScreenEvaluator``. Data-dependent shapes (not suited to a
  GPU). On the CPU it equals ``evaluate`` in every ``OverlapInfo`` field. That needs (found in S3
  Step 0): g computed with one 3x3 gemm per orientation against the ``g_hkl.T`` view, and kept
  component-major (3, B K). A batched matmul rounds 45% of the g rows differently.
* **fixed** (GPU, ``mps`` or ``cuda``; also runnable on the CPU, where it equals ``evaluate``
  too): all 2K omega rows of every orientation, validity as a mask, no data-dependent shape until
  the single packed transfer to the host. On a GPU the transcendental functions differ from the
  CPU's in the last bits, so a GPU result is *not* equal to ``evaluate`` (about 11% of the
  orientations differ in a count field on MPS at K = 112; compare wrong counts, not fields).

The GPU variant is used only for batches of at least ``GPU_MIN_BATCH`` orientations (the measured
crossover against the CPU compact variant, S3 Step 0: 56-62 us at B = 64 against 2.3 ms at B = 1);
smaller batches use the compact variant on the CPU.

Supported: the default rule, ``CppOverlapRule`` and ``CppEdgeRules`` (S1 Step 4: the frame lookup
floors, a detector coordinate below 0 is pixel -1, the vertex-lit rule is in the C call's rule bits;
needs a build with ``HAS_VERTEX_LIT`` >= 2). Not supported (``make_batched_cost_evaluator`` returns
None, the caller falls back to ``evaluate``): the ray rule, soft mode, ``pixel_radius`` > 0, and a
build without the ``_rasterize`` extension.

The binary frames for the C call are the cost function's own ``FrameTable``
(``VoxelCostFunction.frame_table``), shared with ``evaluate``; the omega-to-frame lookup is
``cost_functions.frame_of_omegas``.
"""

import time
import warnings
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import torch

from . import cost_functions as CF
from .cost_functions import OverlapInfo, VoxelCostFunction
from .diffraction_core import get_scattering_omegas_torch

# Smallest batch sent to the GPU variant (S3 Step 0: MPS beats the CPU compact variant from about
# a batch of 64; Step 4 re-checks it).
GPU_MIN_BATCH = 64
DEVICES = ("cpu", "mps", "cuda")

_F32 = torch.float32


def device_available(device: str) -> bool:
    """Whether torch can run on ``device`` ("cpu", "mps" or "cuda") in this process."""
    if device == "cpu":
        return True
    if device == "mps":
        return bool(getattr(torch.backends, "mps", None) and torch.backends.mps.is_available())
    if device == "cuda":
        return bool(torch.cuda.is_available())
    return False


class _DeviceConstants:
    """The cost function's constants as float32 tensors on one device."""

    def __init__(self, cf: VoxelCostFunction, device: torch.device) -> None:
        kw: Dict[str, Any] = dict(dtype=_F32, device=device)
        sim = cf.simulator
        self.beam = sim.beam_direction.to(**kw)
        m = cf.sample.sample_to_lab_matrix
        self.matrix = m.clone()  # the values the constants were taken from (see ``current``)
        self.base_rot = m[:3, :3].to(**kw)
        self.base_t = m[:3, 3].to(**kw)
        self.dets: List[Dict[str, Any]] = []
        for d in cf.detector_list:
            pl = d._detector_plane
            assert pl is not None  # set by the detector's constructor
            self.dets.append(
                dict(
                    n=pl.normal.to(**kw),
                    d=torch.as_tensor(pl.d, **kw),
                    pos=d._position.to(**kw),
                    org=d._lab_frame_coord_origin.to(**kw),
                    bj=d._lab_frame_basis_j.to(**kw),
                    bk=d._lab_frame_basis_k.to(**kw),
                    hw=float(d.pixel_half_width),
                    pw=float(d.pixel_width),
                    hh=float(d.pixel_half_height),
                    ph=float(d.pixel_height),
                )
            )
        # per phase: (K, g_hkl^T, |g|); on the CPU the transposed *view* evaluate multiplies with
        self.phase: Dict[int, Tuple[int, torch.Tensor, torch.Tensor]] = {}
        for p, (g_hkl, g_mag) in cf._phase_recip_vecs.items():
            gT = g_hkl.T if device.type == "cpu" else g_hkl.T.contiguous().to(device)
            self.phase[p] = (g_hkl.shape[0], gT, g_mag.to(device))

    def current(self, cf: VoxelCostFunction) -> bool:
        """False when the sample's matrix changed (also in place) since these were taken."""
        return torch.equal(cf.sample.sample_to_lab_matrix, self.matrix)


class BatchedCostEvaluator:
    """Batched local cost of one ``VoxelCostFunction`` (hard mode, ``pixel_radius`` 0).

    ``evaluate_many(R, voxel_vertices, phase_index)`` returns one ``OverlapInfo`` per orientation
    and adds B to ``cost_fn.eval_count``. ``timing`` accumulates seconds in "torch" (stages A-C and
    the host transfer), "host_prep" and "stage_d" (the C calls); ``counts`` the calls and
    orientations per variant.
    """

    def __init__(self, cost_fn: VoxelCostFunction, device: str = "cpu") -> None:
        reason = unsupported_reason(cost_fn)
        if reason is not None:
            raise NotImplementedError(f"batched local cost: {reason}")
        if device not in DEVICES:
            raise ValueError(f"device must be one of {DEVICES}, got {device!r}")
        if not device_available(device):
            raise RuntimeError(f"SearchDevice {device!r} is not available in this torch build")
        self.cf = cost_fn
        self.device = device
        self._consts: Dict[str, _DeviceConstants] = {}
        rng = cost_fn.range_map
        self._low, self._width = rng.low, rng.width
        self._index = CF.frame_index_array(rng)
        self._n_det = len(cost_fn.detector_list)
        # the binary frames for stage D are the cost function's FrameTable (cf.frame_table())
        self._rule_bits = (1 if cost_fn.cpp_overlap_rule else 0) | (
            2 if cost_fn.cpp_edge_rules else 0
        )
        self.timing: Dict[str, float] = {"torch": 0.0, "host_prep": 0.0, "stage_d": 0.0}
        self.counts: Dict[str, int] = {
            "calls_compact": 0,
            "calls_fixed": 0,
            "orient_compact": 0,
            "orient_fixed": 0,
        }

    # ------------------------------------------------------------------ public

    def evaluate_many(
        self,
        orientations: np.ndarray,
        voxel_vertices: torch.Tensor,
        phase_index: int = 0,
        variant: Optional[str] = None,
    ) -> List[OverlapInfo]:
        """``OverlapInfo`` of each of B orientations (B, 3, 3), as ``cost_fn.evaluate`` returns.

        ``variant``: None picks the GPU variant for batches of at least ``GPU_MIN_BATCH`` on a
        GPU device and the CPU compact variant otherwise; "compact" or "fixed" forces one (the
        fixed one runs on ``self.device``, which is "cpu" for a CPU evaluator).
        """
        cf = self.cf
        R_np = np.asarray(orientations, dtype=np.float64)
        B = len(R_np)
        if B == 0:
            return []
        if phase_index not in cf._phase_recip_vecs:
            return [cf.evaluate(R_np[b], voxel_vertices, phase_index) for b in range(B)]
        if variant is None:
            use_fixed = self.device != "cpu" and B >= GPU_MIN_BATCH
        elif variant in ("compact", "fixed"):
            use_fixed = variant == "fixed"
        else:
            raise ValueError(f"variant must be None, 'compact' or 'fixed', got {variant!r}")
        dev_name = self.device if use_fixed else "cpu"
        consts = self._constants(dev_name)
        K, g_hkl_T, g_mag = consts.phase[phase_index]
        dev = torch.device(dev_name)
        n_det = self._n_det

        t0 = time.perf_counter()
        R = torch.from_numpy(R_np).float().to(dev)
        verts = voxel_vertices.to(dev)
        if use_fixed:
            arr = self._forward_fixed(consts, R, verts, K, g_hkl_T, g_mag).cpu().numpy()
            t1 = time.perf_counter()
            omega = arr[..., 6 * n_det + n_det]
            wedge = self._wedge(omega)
            keep = (arr[..., -1] > 0.5) & (wedge >= 0)
            counts = keep.sum(axis=1)
            sel = arr[keep]
            pix = np.ascontiguousarray(sel[:, : 6 * n_det].reshape(-1, n_det, 3, 2))
            hit = np.ascontiguousarray(sel[:, 6 * n_det : 7 * n_det] > 0.5).astype(np.uint8)
            wedge_k = np.ascontiguousarray(wedge[keep].astype(np.int32))
            self.counts["calls_fixed"] += 1
            self.counts["orient_fixed"] += B
        else:
            got = self._forward_compact(consts, R, verts, B, K, g_hkl_T, g_mag)
            self.counts["calls_compact"] += 1
            self.counts["orient_compact"] += B
            if got is None:  # no observable reflection or peak in the whole batch
                cf.eval_count += B
                return [OverlapInfo() for _ in range(B)]
            pix, hit, wedge_k, counts = got
            t1 = time.perf_counter()
        t2 = time.perf_counter()
        out = self._stage_d(counts, wedge_k, hit, pix)
        t3 = time.perf_counter()
        cf.eval_count += B
        self.timing["torch"] += t1 - t0
        self.timing["host_prep"] += t2 - t1
        self.timing["stage_d"] += t3 - t2
        return out

    # ----------------------------------------------------------------- helpers

    def _constants(self, dev_name: str) -> _DeviceConstants:
        c = self._consts.get(dev_name)
        if c is None or not c.current(self.cf):
            c = _DeviceConstants(self.cf, torch.device(dev_name))
            self._consts[dev_name] = c
        return c

    def _wedge(self, omega: np.ndarray) -> np.ndarray:
        """Frame index of each omega exactly as calculate_diffraction_overlap_batched (-1: none)."""
        return CF.frame_of_omegas(
            omega, self._low, self._width, self._index, self.cf.cpp_edge_rules
        )

    def _stage_d(
        self, counts: np.ndarray, wedge_k: np.ndarray, hit: np.ndarray, pix: np.ndarray
    ) -> List[OverlapInfo]:
        """One stage D call per orientation (the existing C function), in peak order."""
        B = len(counts)
        table = self.cf.frame_table()
        if len(wedge_k):
            table.ensure(wedge_k)
        offs = np.concatenate([[0], np.cumsum(counts)])
        rb = (self._rule_bits,) if self._rule_bits else ()
        out: List[OverlapInfo] = []
        for b in range(B):
            s, e = int(offs[b]), int(offs[b + 1])
            if e == s:
                out.append(OverlapInfo())
                continue
            r = CF._c_stage_d_overlap(
                table.images,
                wedge_k[s:e],
                hit[s:e],
                None,
                pix[s:e],
                self._n_det,
                0,
                e - s,
                table.rows,
                table.cols,
                *rb,
            )
            out.append(
                OverlapInfo(
                    pixel_overlap=r[0],
                    pixel_on_detector=r[1],
                    peak_overlap=r[2],
                    peak_on_detector=r[3],
                    quality=r[4],
                    n_quality_points=r[5],
                    detectors_overlap=r[6],
                )
            )
        return out

    # --------------------------------------------------------- torch stages A-C

    def _bragg(
        self,
        cst: _DeviceConstants,
        R: torch.Tensor,
        K: int,
        g_hkl_T: torch.Tensor,
        g_mag: torch.Tensor,
    ) -> Tuple[torch.Tensor, ...]:
        """g = R g_hkl^T for every orientation, then the Bragg solve. On the CPU this reproduces
        ``evaluate`` bit for bit: one 3x3 @ 3xK gemm per orientation against the transposed view
        of g_hkl (a batched matmul, or the same gemm on a contiguous copy of g_hkl^T, rounds some
        g rows differently), and g kept component-major (3, B K) so the elementwise kernels see
        contiguous columns as in ``evaluate``. On a GPU: one batched matmul."""
        B = R.shape[0]
        if R.device.type == "cpu":
            gT = torch.cat([R[b] @ g_hkl_T for b in range(B)], dim=1)  # (3, B K)
            g = gT.T
        else:
            g = torch.matmul(R, g_hkl_T).transpose(1, 2).reshape(B * K, 3)
        sim = self.cf.simulator
        br = get_scattering_omegas_torch(
            g, g_mag.repeat(B), sim.beam_energy, sim.beam_deflection_chi
        )
        return g.contiguous(), br.observable, br.omega1, br.omega2

    def _rows_eta(
        self, cst: _DeviceConstants, omegas: torch.Tensor, normals: torch.Tensor
    ) -> Tuple[torch.Tensor, ...]:
        """Rz(omega) @ base, reflected beam and the eta mask; any leading shape."""
        cf = self.cf
        c, s = torch.cos(omegas), torch.sin(omegas)
        z, o = torch.zeros_like(c), torch.ones_like(c)
        Rz = torch.stack([c, -s, z, s, c, z, z, z, o], dim=-1).reshape(*c.shape, 3, 3)
        full_rot = Rz @ cst.base_rot
        lab_n = torch.matmul(full_rot, normals.unsqueeze(-1)).squeeze(-1)
        dot = (cst.beam * lab_n).sum(dim=-1, keepdim=True)
        refl = cst.beam - 2.0 * dot * lab_n
        rd = torch.norm(refl, dim=-1)
        safe = torch.where(rd > 0, rd, torch.ones_like(rd))
        eta = torch.atan2(torch.abs(refl[..., 1]) / safe, torch.abs(refl[..., 2]) / safe)
        valid = (eta < cf.eta_limit) & (rd > 0)
        if cf.min_sin_eta > 0.0:
            valid = valid & (torch.sin(eta) >= cf.min_sin_eta)
        return full_rot, refl, valid

    def _stage_bc(
        self,
        cst: _DeviceConstants,
        full_rot: torch.Tensor,
        refl: torch.Tensor,
        verts: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """Stages B and C: lab vertices and detector pixel coordinates; any leading shape.
        The pixel coordinate is ``CF._pixel_coordinate``, the function ``evaluate`` uses (with
        ``CppEdgeRules`` a detector coordinate below 0 is pixel -1). Returns (all three vertices
        hit the detector plane (..., n_det), pixel [col, row] of the three vertices
        (..., n_det, 3, 2))."""
        lead = full_rot.shape[:-2]
        dev = full_rot.device
        edge = self.cf.cpp_edge_rules  # a detector coordinate below 0 is pixel -1 (C++ ToColPixel)
        f4 = torch.zeros(*lead, 4, 4, dtype=_F32, device=dev)
        f4[..., :3, :3] = full_rot
        f4[..., :3, 3] = cst.base_t
        f4[..., 3, 3] = 1.0
        v4 = torch.cat([verts, torch.ones(3, 1, device=dev)], dim=1)  # (3, 4)
        lab = torch.einsum("...ij,vj->...vi", f4, v4)[..., :3]  # (..., 3, 3)
        dirs = refl.unsqueeze(-2).expand(*lead, 3, 3)
        hit_all, pix = [], []
        for d in cst.dets:
            denom = (dirs * d["n"]).sum(dim=-1)
            numer = -((lab * d["n"]).sum(dim=-1) + d["d"])
            parallel = torch.abs(denom) < 1e-8
            safe_d = torch.where(parallel, torch.ones_like(denom), denom)
            t = torch.where(parallel, torch.zeros_like(denom), numer / safe_d)
            hits = (~parallel) & (t > 0)
            rel = lab + t.unsqueeze(-1) * dirs - d["pos"]
            loc = rel - d["org"]
            j = (loc * d["bj"]).sum(dim=-1)
            k = (loc * d["bk"]).sum(dim=-1)
            col = CF._pixel_coordinate(j, d["hw"], d["pw"], edge)
            row = CF._pixel_coordinate(k, d["hh"], d["ph"], edge)
            hit_all.append(hits.all(dim=-1))
            pix.append(torch.stack([col, row], dim=-1))
        return torch.stack(hit_all, dim=-1), torch.stack(pix, dim=-3)

    # ------------------------------------------------------------ the variants

    def _forward_fixed(
        self,
        cst: _DeviceConstants,
        R: torch.Tensor,
        verts: torch.Tensor,
        K: int,
        g_hkl_T: torch.Tensor,
        g_mag: torch.Tensor,
    ) -> torch.Tensor:
        """Fixed-shape forward: (B, 2K, 6 n_det + n_det + 2) float32 with the pixel coordinates,
        the detector hits, the omega and the validity (observable and eta) of every omega row."""
        B = R.shape[0]
        g, obs, om1, om2 = self._bragg(cst, R, K, g_hkl_T, g_mag)
        normals = (g / torch.norm(g, dim=1, keepdim=True)).reshape(B, K, 3)
        omegas = torch.cat([om1.reshape(B, K), om2.reshape(B, K)], dim=1)  # (B, 2K)
        normals = torch.cat([normals, normals], dim=1)
        obs2 = torch.cat([obs.reshape(B, K)] * 2, dim=1)
        full_rot, refl, valid = self._rows_eta(cst, omegas, normals)
        valid = valid & obs2
        det_hit, pix = self._stage_bc(cst, full_rot, refl, verts)
        return torch.cat(
            [
                pix.reshape(B, 2 * K, -1),
                det_hit.float(),
                omegas.unsqueeze(-1),
                valid.float().unsqueeze(-1),
            ],
            dim=-1,
        )

    def _forward_compact(
        self,
        cst: _DeviceConstants,
        R: torch.Tensor,
        verts: torch.Tensor,
        B: int,
        K: int,
        g_hkl_T: torch.Tensor,
        g_mag: torch.Tensor,
    ) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
        """CPU compact forward: (pixels, detector hits, wedge, peaks per orientation) of the
        surviving peaks ordered by orientation, or None when no orientation has a peak."""
        g, obs, om1, om2 = self._bragg(cst, R, K, g_hkl_T, g_mag)
        if not bool(obs.any()):
            return None
        og = g[obs]
        normals = og / torch.norm(og, dim=1, keepdim=True)
        seg1 = torch.arange(B).repeat_interleave(K)[obs]
        omegas = torch.cat([om1[obs], om2[obs]])
        normals = torch.cat([normals, normals])
        seg = torch.cat([seg1, seg1])
        full_rot, refl, valid = self._rows_eta(cst, omegas, normals)
        vidx = torch.nonzero(valid).squeeze(1)
        om_np = omegas[vidx].numpy()
        wedge_all = self._wedge(om_np)
        ok = wedge_all >= 0
        if not ok.any():
            return None
        sel_i = vidx[torch.from_numpy(np.nonzero(ok)[0])]
        seg_np = seg[sel_i].numpy()
        order = np.argsort(seg_np, kind="stable")
        sel_i = sel_i[torch.from_numpy(order)]
        seg_np = seg_np[order]
        wedge_k = np.ascontiguousarray(wedge_all[ok][order].astype(np.int32))
        det_hit_t, pix_t = self._stage_bc(cst, full_rot[sel_i], refl[sel_i], verts)
        pix = np.ascontiguousarray(pix_t.numpy())
        hit = np.ascontiguousarray(det_hit_t.numpy().astype(np.uint8))
        counts = np.bincount(seg_np, minlength=B)
        return pix, hit, wedge_k, counts


def unsupported_reason(cost_fn: VoxelCostFunction) -> Optional[str]:
    """Why ``BatchedCostEvaluator`` cannot score ``cost_fn``'s costs (None: it can)."""
    if cost_fn.pixel_radius != 0:
        return "pixel_radius > 0 (the screen cost; use ScreenEvaluator)"
    if cost_fn.mode != "hard":
        return "soft mode"
    if cost_fn.ray_rule is not None:
        return "the ray rule"
    if not CF._HAS_C_RASTERIZE:
        return "the _rasterize C extension is missing"
    if cost_fn.cpp_edge_rules and not CF._C_HAS_VERTEX_LIT:
        return "the _rasterize extension predates the exact CppEdgeRules (HAS_VERTEX_LIT)"
    if cost_fn.cpp_overlap_rule and not CF._C_HAS_CPP_RULE:
        return "the _rasterize extension predates CppOverlapRule"
    if cost_fn.eta_limit is None:
        return "no eta limit"
    return None


def make_batched_cost_evaluator(
    cost_fn: VoxelCostFunction, device: str = "cpu", warn: bool = True
) -> Optional[BatchedCostEvaluator]:
    """A ``BatchedCostEvaluator``, or None when the cost function is not supported (the caller
    falls back to ``evaluate``; a RuntimeWarning names the reason when ``warn``). A device that is
    not available raises RuntimeError (silently using the CPU would change the results)."""
    reason = unsupported_reason(cost_fn)
    if reason is not None:
        if warn:
            warnings.warn(
                f"BatchedQuickMC is not supported for this cost function ({reason}): "
                "using the sequential quick MC",
                RuntimeWarning,
                stacklevel=2,
            )
        return None
    return BatchedCostEvaluator(cost_fn, device)
