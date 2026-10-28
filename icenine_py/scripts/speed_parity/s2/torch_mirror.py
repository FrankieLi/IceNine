"""The torch operations of ``VoxelCostFunction.evaluate`` (hard mode, no ray rule) statement by
statement, with every intermediate kept, for the S2 Step 0 per-intermediate comparison.

``mirror(cf, R, verts, phase)`` runs the same tensor operations on tensors of the same shapes as
``cost_functions.py`` (``evaluate`` 1151-1260 and ``calculate_diffraction_overlap_batched``
624-886 on the C fast path) and returns the ``OverlapInfo`` tuple and the intermediates, keyed by
row ``s * K + k`` (solution s of reflection k). ``s2_parity.py`` checks that the mirror's result
equals the library's on every orientation, so the intermediates are those of the library.
"""

import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))

from icenine import cost_functions as CF  # noqa: E402
from icenine.diffraction_core import BraggConstants, scattering_omegas_fixed_beam  # noqa: E402


def reference64(cf: Any, orientation: np.ndarray, voxel_vertices: torch.Tensor, phase: int) -> Dict:
    """The same quantities in float64 (numpy) from the same float32 inputs (R rounded to float32,
    g_hkl, the Bragg constants, the geometry), for every row ``s * K + k``: a reference that tells
    which of two float32 results is nearer the exact value. No filtering (all rows computed).
    With ``cpp_edge_rules`` a negative detector coordinate gives pixel -1, as in the library."""
    f = np.float64
    if cf.center_screen and cf.pixel_radius > 0:
        voxel_vertices = voxel_vertices.mean(dim=0, keepdim=True).repeat(3, 1)
    V = voxel_vertices.numpy().astype(f)
    g_hkl_t, g_mag_t = cf._phase_recip_vecs[phase]
    consts = BraggConstants(g_mag_t, cf.simulator.beam_energy, cf.simulator.beam_deflection_chi)
    R = torch.from_numpy(orientation).float().numpy().astype(f)
    g = g_hkl_t.numpy().astype(f) @ R.T  # (K, 3)
    K = len(g)
    eps = float(np.float32(consts.epsilon))
    cos_chi = np.clip(g[:, 2] / consts.g_mag_eps.numpy().astype(f), -1, 1)
    sin_chi = np.sqrt(1 - cos_chi * cos_chi + eps)
    num = consts.sin_theta.numpy().astype(f) + cos_chi * float(consts.sin_chi_laue.reshape(-1)[0])
    den = sin_chi * float(consts.cos_chi_laue.reshape(-1)[0]) + eps
    d0 = np.arctan2(g[:, 0], g[:, 1])
    b1 = np.arcsin(np.clip(num / den, -1, 1))
    o1, o2 = b1 + d0, np.pi - b1 + d0
    o1 = np.where(o1 > np.pi, o1 - 2 * np.pi, o1)
    o1 = np.where(o1 < -np.pi, o1 + 2 * np.pi, o1)
    o2 = np.where(o2 > np.pi, o2 - 2 * np.pi, o2)
    o2 = np.where(o2 < -np.pi, o2 + 2 * np.pi, o2)
    om = np.concatenate([o1, o2])
    n = g / np.linalg.norm(g, axis=1, keepdims=True)
    n = np.concatenate([n, n])
    c, s = np.cos(om), np.sin(om)
    Rz = np.zeros((2 * K, 3, 3))
    Rz[:, 0, 0], Rz[:, 0, 1], Rz[:, 1, 0], Rz[:, 1, 1], Rz[:, 2, 2] = c, -s, s, c, 1
    s2l = cf.sample.sample_to_lab_matrix.numpy().astype(f)
    F = Rz @ s2l[:3, :3]
    ln = np.einsum("nij,nj->ni", F, n)
    beam = cf.simulator.beam_direction.numpy().astype(f)
    refl = beam[None] - 2 * (ln @ beam)[:, None] * ln
    rn = np.linalg.norm(refl, axis=1)
    eta = np.arctan2(np.abs(refl[:, 1]) / rn, np.abs(refl[:, 2]) / rn)
    LV = np.einsum("nij,vj->nvi", F, V) + s2l[:3, 3][None, None]
    n_det = len(cf.detector_list)
    pix = np.zeros((2 * K, n_det, 3, 2))
    t_all = np.zeros((2 * K, n_det, 3))
    for d, det in enumerate(cf.detector_list):
        pl = det._detector_plane
        pn = pl.normal.numpy().astype(f)
        denom = refl @ pn
        numer = -(LV @ pn + float(pl.d))
        t = numer / denom[:, None]
        pts = LV + t[:, :, None] * refl[:, None, :]
        loc = pts - det._position.numpy().astype(f) - det._lab_frame_coord_origin.numpy().astype(f)
        j = loc @ det._lab_frame_basis_j.numpy().astype(f)
        k = loc @ det._lab_frame_basis_k.numpy().astype(f)
        pix[:, d, :, 0] = (j + float(np.float32(det.pixel_half_width))) / float(
            np.float32(det.pixel_width)
        )
        pix[:, d, :, 1] = (k + float(np.float32(det.pixel_half_height))) / float(
            np.float32(det.pixel_height)
        )
        if cf.cpp_edge_rules:
            # CppEdgeRules -1 rule, as cost_functions._pixel_coordinate (C++ ToColPixel /
            # ToRowPixel, Detector.cpp 528-548): a negative detector coordinate is pixel -1
            pix[:, d, :, 0] = np.where(j < 0, -1.0, pix[:, d, :, 0])
            pix[:, d, :, 1] = np.where(k < 0, -1.0, pix[:, d, :, 1])
        t_all[:, d] = t
    return dict(K=K, omega=om, eta=eta, pix=pix, t=t_all, obs_margin=np.abs(num) - np.abs(den))


def mirror(cf: Any, orientation: np.ndarray, voxel_vertices: torch.Tensor, phase: int) -> Dict:
    assert cf.mode == "hard" and cf.ray_rule is None
    if cf.center_screen and cf.pixel_radius > 0:
        voxel_vertices = voxel_vertices.mean(dim=0, keepdim=True).repeat(3, 1)
    g_hkl_batch, g_mag_batch = cf._phase_recip_vecs[phase]
    K = g_hkl_batch.shape[0]
    n_det = len(cf.detector_list)
    out: Dict[str, Any] = dict(K=K)
    orientation_t = torch.from_numpy(orientation).float()
    g_lab_batch = (orientation_t @ g_hkl_batch.T).T
    out["g_lab"] = g_lab_batch.numpy().copy()
    consts = BraggConstants(g_mag_batch, cf.simulator.beam_energy, cf.simulator.beam_deflection_chi)
    br = scattering_omegas_fixed_beam(g_lab_batch, consts)
    obs_mask = br.observable
    out["obs"] = obs_mask.numpy().astype(np.uint8)
    out["omega1"] = br.omega1.numpy().copy()
    out["omega2"] = br.omega2.numpy().copy()
    # Bragg margins: |numerator| vs |denominator| (recomputed with the same operations)
    eps = consts.epsilon
    cos_chi = torch.clamp(g_lab_batch[:, 2] / consts.g_mag_eps, -1.0, 1.0)
    sin_chi = torch.sqrt(1.0 - cos_chi * cos_chi + eps)
    out["bragg_num"] = (consts.sin_theta + cos_chi * consts.sin_chi_laue).abs().numpy()
    out["bragg_den"] = (sin_chi * consts.cos_chi_laue + eps).abs().numpy()
    rows = 2 * K
    out["eta"] = np.full(rows, np.nan, np.float32)
    out["valid"] = np.zeros(rows, np.uint8)
    out["frame"] = np.full(rows, -1, np.int32)
    out["omega"] = np.full(rows, np.nan, np.float32)
    out["pix"] = np.full((rows, n_det, 3, 2), np.nan, np.float32)
    out["hit"] = np.zeros((rows, n_det), np.uint8)
    out["t"] = np.full((rows, n_det, 3), np.nan, np.float32)
    out["denom"] = np.full((rows, n_det), np.nan, np.float32)
    out["jk"] = np.full((rows, n_det, 3, 2), np.nan, np.float32)  # detector coordinates (mm)
    out["result"] = (0, 0, 0, 0, 0.0, 0, 0)
    if not obs_mask.any():
        return out
    obs_k = torch.nonzero(obs_mask).reshape(-1).numpy()
    beam_dir = cf.simulator.beam_direction
    base_rot = cf.sample.sample_to_lab_matrix[:3, :3]
    obs_g = g_lab_batch[obs_mask]
    obs_mag = torch.norm(obs_g, dim=1, keepdim=True)
    obs_normals = obs_g / obs_mag
    all_omegas = torch.cat([br.omega1[obs_mask], br.omega2[obs_mask]])
    all_normals = torch.cat([obs_normals, obs_normals])
    row_of = np.concatenate([obs_k, K + obs_k])  # row key of each of the 2K_obs rows
    N = all_omegas.shape[0]
    cos_w = torch.cos(all_omegas)
    sin_w = torch.sin(all_omegas)
    Rz = torch.zeros(N, 3, 3)
    Rz[:, 0, 0] = cos_w
    Rz[:, 0, 1] = -sin_w
    Rz[:, 1, 0] = sin_w
    Rz[:, 1, 1] = cos_w
    Rz[:, 2, 2] = 1.0
    full_rot = Rz @ base_rot
    lab_normals = torch.bmm(full_rot, all_normals.unsqueeze(-1)).squeeze(-1)
    dot = (beam_dir.unsqueeze(0) * lab_normals).sum(dim=1, keepdim=True)
    reflected = beam_dir.unsqueeze(0) - 2.0 * dot * lab_normals
    rd_norms = torch.norm(reflected, dim=1)
    safe_norms = torch.where(rd_norms > 0, rd_norms, torch.ones_like(rd_norms))
    ry = torch.abs(reflected[:, 1]) / safe_norms
    rz_val = torch.abs(reflected[:, 2]) / safe_norms
    eta = torch.atan2(ry, rz_val)
    valid = (eta < cf.eta_limit) & (rd_norms > 0)
    sin_eta = torch.sin(eta)
    if cf.min_sin_eta > 0.0:
        valid = valid & (sin_eta >= cf.min_sin_eta)
    out["omega"][row_of] = all_omegas.numpy()
    out["eta"][row_of] = eta.numpy()
    out["sin_eta"] = np.full(rows, np.nan, np.float32)
    out["sin_eta"][row_of] = sin_eta.numpy()
    out["valid"][row_of] = valid.numpy().astype(np.uint8)
    peak_omegas = all_omegas[valid]
    peak_normals = all_normals[valid]
    vrow = row_of[valid.numpy()]
    if peak_omegas.shape[0] == 0:
        return out

    # calculate_diffraction_overlap_batched, fast path
    table = cf.frame_table()
    frame_index = cf._frame_index
    rmap = cf.range_map
    wedge_indices = CF.frame_of_omegas(
        peak_omegas.numpy(), rmap.low, rmap.width, frame_index, cf.cpp_edge_rules
    )
    out["frame_f"] = np.full(rows, np.nan, np.float64)
    out["frame_f"][vrow] = (peak_omegas.numpy() - rmap.low) / rmap.width
    valid_mask = wedge_indices >= 0
    out["frame"][vrow] = np.where(valid_mask, wedge_indices, -1)
    if not np.any(valid_mask):
        return out
    valid_indices = np.where(valid_mask)[0]
    v_omegas = peak_omegas[valid_indices]
    v_normals = peak_normals[valid_indices]
    v_wedge = wedge_indices[valid_indices]
    mrow = vrow[valid_indices]
    M = len(valid_indices)
    orig_matrix = cf.sample.sample_to_lab_matrix
    cos_w = torch.cos(v_omegas)
    sin_w = torch.sin(v_omegas)
    Rz = torch.zeros(M, 3, 3)
    Rz[:, 0, 0] = cos_w
    Rz[:, 0, 1] = -sin_w
    Rz[:, 1, 0] = sin_w
    Rz[:, 1, 1] = cos_w
    Rz[:, 2, 2] = 1.0
    base_rot = orig_matrix[:3, :3]
    full_rot = Rz @ base_rot
    lab_normals = torch.bmm(full_rot, v_normals.unsqueeze(-1)).squeeze(-1)
    dot = (beam_dir.unsqueeze(0) * lab_normals).sum(dim=1, keepdim=True)
    reflected = beam_dir.unsqueeze(0) - 2.0 * dot * lab_normals
    full_4x4 = torch.zeros(M, 4, 4)
    full_4x4[:, :3, :3] = full_rot
    full_4x4[:, :3, 3] = orig_matrix[:3, 3]
    full_4x4[:, 3, 3] = 1.0
    verts_4d = torch.cat([voxel_vertices, torch.ones(3, 1)], dim=1)
    lab_verts = torch.einsum("mij,vj->mvi", full_4x4, verts_4d)[:, :, :3]
    per_det_all_hit = []
    per_det_pixels = []
    for det_idx in range(n_det):
        detector = cf.detector_list[det_idx]
        plane = detector._detector_plane
        plane_n = plane.normal
        plane_d = plane.d
        origins = lab_verts.reshape(M * 3, 3)
        dirs = reflected.unsqueeze(1).expand(M, 3, 3).reshape(M * 3, 3)
        denom = (dirs * plane_n).sum(dim=1)
        numer = -((origins * plane_n).sum(dim=1) + plane_d)
        parallel = torch.abs(denom) < 1e-8
        safe_denom = torch.where(parallel, torch.ones_like(denom), denom)
        t = torch.where(parallel, torch.zeros_like(denom), numer / safe_denom)
        hits = (~parallel) & (t > 0)
        pts = origins + t.unsqueeze(1) * dirs
        relative = pts - detector._position
        pixel_loc = relative - detector._lab_frame_coord_origin
        j = (pixel_loc * detector._lab_frame_basis_j).sum(dim=1)
        k = (pixel_loc * detector._lab_frame_basis_k).sum(dim=1)
        col = CF._pixel_coordinate(
            j, detector.pixel_half_width, detector.pixel_width, cf.cpp_edge_rules
        )
        row = CF._pixel_coordinate(
            k, detector.pixel_half_height, detector.pixel_height, cf.cpp_edge_rules
        )
        hits_mv = hits.reshape(M, 3)
        all_hit = hits_mv.all(dim=1)
        pixels = torch.stack([col.reshape(M, 3), row.reshape(M, 3)], dim=2)
        per_det_all_hit.append(all_hit)
        per_det_pixels.append(pixels)
        out["pix"][mrow, det_idx] = pixels.numpy()
        out["hit"][mrow, det_idx] = all_hit.numpy().astype(np.uint8)
        out["t"][mrow, det_idx] = t.reshape(M, 3).numpy()
        out["denom"][mrow, det_idx] = denom.reshape(M, 3)[:, 0].numpy()
        out["jk"][mrow, det_idx] = torch.stack([j.reshape(M, 3), k.reshape(M, 3)], dim=2).numpy()
    v_wedge_np = np.asarray(v_wedge, dtype=np.int32)
    table.ensure(v_wedge_np)
    det_hit = np.stack([per_det_all_hit[d].numpy().astype(np.uint8) for d in range(n_det)], axis=1)
    rule_bits = (1 if cf.cpp_overlap_rule else 0) | (2 if cf.cpp_edge_rules else 0)
    if cf.pixel_radius > 0:
        centers = np.ascontiguousarray(
            np.stack(
                [per_det_pixels[d][:, 0, :].numpy().astype(np.int32) for d in range(n_det)],
                axis=1,
            )
        )
        res = CF._c_stage_d_overlap(
            table.images, v_wedge_np, det_hit, centers, None, n_det, cf.pixel_radius, M,
            table.rows, table.cols, *((rule_bits,) if rule_bits else ()),
        )  # fmt: skip
    else:
        verts = np.ascontiguousarray(
            np.stack([per_det_pixels[d].numpy().astype(np.float32) for d in range(n_det)], axis=1)
        )
        res = CF._c_stage_d_overlap(
            table.images, v_wedge_np, det_hit, None, verts, n_det, cf.pixel_radius, M,
            table.rows, table.cols, *((rule_bits,) if rule_bits else ()),
        )  # fmt: skip
    out["result"] = tuple(res)
    out["stage_d_rows"] = mrow
    return out
