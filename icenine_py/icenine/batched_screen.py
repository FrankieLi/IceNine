"""Batched level-0 screen: many orientations per call, exact or with C++'s omega interpolation.

The coarse screen of the discrete search scores every (FZ point, local-grid point) pair with the
pixel-radius cost (``VoxelCostFunction`` with ``pixel_radius`` > 0) and keeps the pair when
``peak_overlap`` > 0. ``VoxelCostFunction.evaluate`` does this one orientation at a time; about
0.5 ms of each call is Python and torch overhead (the C overlap stage is 5% of it).
``ScreenEvaluator`` scores many orientations in one call (config key ``BatchedScreen``, default on
since 2026-10-13), and optionally computes the omegas of the local-grid orientations the way C++
does (opt-in config key ``CppOmegaInterp``):

* exact (``overlaps_exact``): the same arithmetic as ``evaluate``, with the per-orientation torch
  operations stacked into one batch. Pass sets can differ from the unbatched screen only where a
  float32 operation gives a different last bit for another batch shape (counted in
  MIGRATION_HISTORY, "Task A").
* interpolated (``overlaps_interp``): C++ ``LocalAngularInterpProcess``. For each FZ point the
  observable reflections and their two omegas are solved once at the FZ centre
  (``HEDM::PrecomputeScatteringVector``, DiffractionCore.h:262); for each local-grid rotation L the
  reflection direction is rotated exactly (``L @ dir``) but the omega is first-order interpolated,
  ``omega' = omega + (Rz(omega)[0] . G' - x_lab) / (Rz(omega)[1] . G')`` with ``G' = L @ dir`` and
  ``x_lab = (Rz(omega) @ dir)_x`` (``HEDM::InterpolateScatteringOmegas`` and
  ``InterpolateObservablePeaks``, DiffractionCore.h:187-251 and 322-348). A peak whose
  ``|omega' - omega|`` exceeds 1 radian is dropped, and the observable set is frozen at the centre.

C++ reference: Src/DiscreteSearch.h:145-230, Src/DiffractionCore.h:187-348,
Src/CostFunctions.h:197-228 (InterpolatedVoxelCostFunction).
"""

import math
from typing import Any, Optional, Tuple

import numpy as np
import torch

from . import cost_functions as CF
from .diffraction_core import get_scattering_omegas_torch


class ScreenEvaluator:
    """Batched pixel-radius screen for one ``VoxelCostFunction`` (a level's global cost function).

    Supports ``mode='hard'`` without the ray rule and needs the ``_rasterize`` extension; the
    reconstructor falls back to the per-orientation path otherwise.
    """

    def __init__(self, cost_fn: "CF.VoxelCostFunction") -> None:
        if cost_fn.pixel_radius <= 0:
            raise ValueError("ScreenEvaluator needs a pixel-radius (screen) cost function")
        if cost_fn.ray_rule is not None or cost_fn.mode != "hard" or not CF._HAS_C_RASTERIZE:
            raise NotImplementedError("batched screen: hard mode, no ray rule, C extension")
        self.cf = cost_fn
        rng = cost_fn.range_map
        n_bins = len(rng.index_list)
        self._index = np.array(
            [-1 if i is None else int(i) for i in rng.index_list], dtype=np.int64
        )
        self._n_bins = n_bins
        self._rule_bits = (1 if cost_fn.cpp_overlap_rule else 0) | (
            2 if cost_fn.cpp_edge_rules else 0
        )

    # ------------------------------------------------------------------ public

    def overlaps_exact(
        self, orientations: np.ndarray, voxel_vertices: torch.Tensor, phase_index: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """(peak_overlap, peak_on_detector), each (K,), of K orientations (K, 3, 3)."""
        cf = self.cf
        K = len(orientations)
        if phase_index not in cf._phase_recip_vecs:
            return np.zeros(K, np.int64), np.zeros(K, np.int64)
        g_hkl, g_mag = cf._phase_recip_vecs[phase_index]
        n = g_hkl.shape[0]
        R = torch.from_numpy(np.asarray(orientations, dtype=np.float64)).float()  # (K, 3, 3)
        g = (
            torch.matmul(R, g_hkl.T).transpose(1, 2).reshape(K * n, 3)
        )  # (K*n, 3), orientation-major
        br = get_scattering_omegas_torch(
            g, g_mag.repeat(K), cf.simulator.beam_energy, cf.simulator.beam_deflection_chi
        )
        obs = br.observable
        if not obs.any():
            return np.zeros(K, np.int64), np.zeros(K, np.int64)
        og = g[obs]
        normals = og / torch.norm(og, dim=1, keepdim=True)
        seg_one = torch.arange(K).repeat_interleave(n)[obs]
        omegas = torch.cat([br.omega1[obs], br.omega2[obs]])
        normals = torch.cat([normals, normals])
        seg = torch.cat([seg_one, seg_one])
        return self._counts(omegas, normals, seg, K, voxel_vertices)

    def interp_peaks(
        self, centres: np.ndarray, local_grid: np.ndarray, phase_index: int
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """The peaks of every (centre, local rotation) pair as C++ builds them: (omegas (M,),
        normals (M, 3), segment (M,)) with segment = centre index * J + local index. The omegas are
        interpolated from the centre, the normals are the exactly rotated directions."""
        cf = self.cf
        Kc, J = len(centres), len(local_grid)
        empty = (torch.zeros(0), torch.zeros(0, 3), torch.zeros(0, dtype=torch.long))
        if phase_index not in cf._phase_recip_vecs:
            return empty
        g_hkl, g_mag = cf._phase_recip_vecs[phase_index]
        n = g_hkl.shape[0]
        C = torch.from_numpy(np.asarray(centres, dtype=np.float64)).float()
        L = torch.from_numpy(np.asarray(local_grid, dtype=np.float64)).float()  # (J, 3, 3)
        g = torch.matmul(C, g_hkl.T).transpose(1, 2).reshape(Kc * n, 3)
        br = get_scattering_omegas_torch(
            g, g_mag.repeat(Kc), cf.simulator.beam_energy, cf.simulator.beam_deflection_chi
        )
        obs = br.observable
        if not obs.any():
            return empty
        og = g[obs]  # (R, 3)
        d0 = og / torch.norm(og, dim=1, keepdim=True)  # unit direction at the centre (oRecipVecDir)
        cidx = torch.arange(Kc).repeat_interleave(n)[obs]  # centre of each observed reflection
        # directions after the local rotation: (R, J, 3)
        dn = torch.einsum("jab,rb->rja", L, d0)
        peaks_o, peaks_n, peaks_s = [], [], []
        for omega0 in (br.omega1[obs], br.omega2[obs]):
            c0, s0 = torch.cos(omega0), torch.sin(omega0)  # BuildActiveEulerMatrix(omega, 0, 0)
            x_lab = c0 * d0[:, 0] + (-s0) * d0[:, 1] + 0.0 * d0[:, 2]  # (Rz @ dir)_x
            gx, gy, gz = dn[..., 0], dn[..., 1], dn[..., 2]
            # InterpolateScatteringOmegas: fTop = m00 Gx - x_lab + m01 Gy + m02 Gz,
            # fBottom = row 1 . G
            top = ((c0[:, None] * gx - x_lab[:, None]) + (-s0)[:, None] * gy) + 0.0 * gz
            bottom = (s0[:, None] * gx + c0[:, None] * gy) + 0.0 * gz
            d_om = top / bottom
            keep = torch.abs(d_om) <= 1.0  # also drops NaN and inf
            new_om = omega0[:, None] + d_om
            ridx, jidx = torch.nonzero(keep, as_tuple=True)
            peaks_o.append(new_om[ridx, jidx])
            peaks_n.append(dn[ridx, jidx])
            peaks_s.append(cidx[ridx] * J + jidx)
        return torch.cat(peaks_o), torch.cat(peaks_n), torch.cat(peaks_s)

    def overlaps_interp(
        self,
        centres: np.ndarray,
        local_grid: np.ndarray,
        voxel_vertices: torch.Tensor,
        phase_index: int,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """(peak_overlap, peak_on_detector), each (Kc, J): FZ centres (Kc, 3, 3), local rotations
        (J, 3, 3), omegas interpolated from the centre as C++ does."""
        Kc, J = len(centres), len(local_grid)
        omegas, normals, seg = self.interp_peaks(centres, local_grid, phase_index)
        po, pon = self._counts(omegas, normals, seg, Kc * J, voxel_vertices)
        return po.reshape(Kc, J), pon.reshape(Kc, J)

    # ----------------------------------------------------------------- shared

    def _counts(
        self,
        omegas: torch.Tensor,
        normals: torch.Tensor,
        seg: torch.Tensor,
        n_seg: int,
        voxel_vertices: torch.Tensor,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Peak counts per segment from a flat peak list (omega, unit normal, segment index)."""
        cf = self.cf
        po = np.zeros(n_seg, np.int64)
        pon = np.zeros(n_seg, np.int64)
        if omegas.numel() == 0:
            return po, pon
        sim, sample = cf.simulator, cf.sample
        beam_dir = sim.beam_direction
        base_rot = sample.sample_to_lab_matrix[:3, :3]

        # eta filter, as VoxelCostFunction.evaluate (same operations on the stacked peaks)
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

        # frame lookup (Stage A of calculate_diffraction_overlap_batched)
        bins = CF._omega_bins(
            omegas.numpy(), cf.range_map.low, cf.range_map.width, cf.cpp_edge_rules
        )
        ok_bin = (bins >= 0) & (bins < self._n_bins)
        wedge = np.where(ok_bin, self._index[np.clip(bins, 0, self._n_bins - 1)], -1)
        keep = np.nonzero(valid.numpy() & (wedge >= 0))[0]
        if keep.size == 0:
            return po, pon
        # peaks sorted by segment so each segment is a contiguous slice
        order = keep[np.argsort(seg.numpy()[keep], kind="stable")]
        o_idx = torch.from_numpy(order)
        v_wedge = wedge[order]
        v_seg = seg.numpy()[order]
        M = len(order)
        full_rot_v = full_rot[o_idx]
        reflected_v = reflected[o_idx]

        # Stage B/C of calculate_diffraction_overlap_batched (same operations)
        orig_matrix = sample.sample_to_lab_matrix
        full_4x4 = torch.zeros(M, 4, 4)
        full_4x4[:, :3, :3] = full_rot_v
        full_4x4[:, :3, 3] = orig_matrix[:3, 3]
        full_4x4[:, 3, 3] = 1.0
        if cf.center_screen:
            voxel_vertices = voxel_vertices.mean(dim=0, keepdim=True).repeat(3, 1)
        verts_4d = torch.cat([voxel_vertices, torch.ones(3, 1)], dim=1)
        lab_verts = torch.einsum("mij,vj->mvi", full_4x4, verts_4d)[:, :, :3]
        n_det = len(cf.detector_list)
        det_hit, centers = [], []
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
            row = CF._pixel_coordinate(
                k, det.pixel_half_height, det.pixel_height, cf.cpp_edge_rules
            )
            det_hit.append(hits.reshape(M, 3).all(dim=1).numpy().astype(np.uint8))
            pix = torch.stack([col.reshape(M, 3), row.reshape(M, 3)], dim=2)  # (M, 3, 2)
            centers.append(pix[:, 0, :].numpy().astype(np.int32))
        det_hit_a = np.ascontiguousarray(np.stack(det_hit, axis=1))
        centers_a = np.ascontiguousarray(np.stack(centers, axis=1))  # (M, n_det, 2)

        # images of the frames that occur, then one C call per segment
        unique = np.unique(v_wedge)
        local_of = np.full(int(v_wedge.max()) + 1, -1, dtype=np.int32)
        local_of[unique] = np.arange(len(unique), dtype=np.int32)
        local_wedge = local_of[v_wedge].astype(np.int32)
        image_list = [
            cf.exp_data.get_image(int(w), d).get_binary_numpy()
            for w in unique
            for d in range(n_det)
        ]
        first = cf.exp_data.get_image(int(unique[0]), 0)
        rows, cols = first.num_rows, first.num_cols
        bounds = np.nonzero(np.diff(v_seg))[0] + 1
        starts = np.concatenate([[0], bounds])
        ends = np.concatenate([bounds, [M]])
        for s, e in zip(starts, ends):
            res = CF._c_stage_d_overlap(
                image_list,
                local_wedge[s:e],
                det_hit_a[s:e],
                centers_a[s:e],
                None,
                n_det,
                cf.pixel_radius,
                int(e - s),
                rows,
                cols,
                *((self._rule_bits,) if self._rule_bits else ()),
            )
            po[v_seg[s]] = res[2]
            pon[v_seg[s]] = res[3]
        return po, pon


def make_screen_evaluator(cost_fn: "CF.VoxelCostFunction") -> Optional[ScreenEvaluator]:
    """A ScreenEvaluator, or None when the cost function is not supported (caller falls back)."""
    try:
        return ScreenEvaluator(cost_fn)
    except NotImplementedError:
        return None
