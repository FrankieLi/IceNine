"""Check the per-peak sensitivities behind the angular-resolution estimate, on Example2 voxel 0.

Frame:  grad omega*  = -m / (z.m)                       (verified earlier)
Pixels: grad s (px/rad) = (L tan2theta / (a cos theta)) [t + (z.t) grad omega*],  s = arc length along the ring
with g-hat, m = (g x k)/|g x k|, t = g-hat x m in the lab frame at omega*.
Pixel gradient is checked by finite differences of the simulator's projected spot centroid.
"""

import sys
from pathlib import Path

import numpy as np
import torch

ICENINE_PY = Path("/Users/sfli/Research/IceNine/icenine_py")
sys.path.insert(0, str(ICENINE_PY / "scripts"))
from generate_toy_orientation_dataset import setup_example  # noqa: E402

from icenine.constants import KEV_OVER_HBAR_C_IN_ANG  # noqa: E402
from icenine.diffraction_core import get_scattering_omegas_torch  # noqa: E402
from icenine.orientation_nn import _project_peak_on_detector, _restore_and_rotate, define_roi_set  # noqa: E402
from icenine.peak_filters import TrivialAcceptFn  # noqa: E402

DEG = np.pi / 180
FRAME = 1.0 * DEG


def rz(w):
    c, s = np.cos(w), np.sin(w)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def expm(v):
    th = np.linalg.norm(v)
    if th == 0:
        return np.eye(3)
    k = v / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


example_dir = ICENINE_PY.parent / "Examples" / "Example2.ThreeVoxels"
mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, get_vertices = setup_example(example_dir)
E = float(exp_setup.beam_energy)
k_vec = np.array([KEV_OVER_HBAR_C_IN_ANG * E, 0.0, 0.0])
zhat = np.array([0.0, 0.0, 1.0])

voxel = mic.voxels[0]
R_nom = voxel.orientation.astype(np.float64)
verts = get_vertices(voxel)
roi = define_roi_set(torch.from_numpy(R_nom).float(), verts, sample, detector_list, range_map,
                     exp_setup, structure_list, simulator, phase_index=voxel.phase)
base = sample.sample_to_lab_matrix[:3, :3].clone()
accept = TrivialAcceptFn()
a_pix = detector_list[0].pixel_width
L = [float(d._position[0]) for d in detector_list]
r_perp = float(np.hypot(*voxel.position[:2]))
print(f"voxel 0: {len(roi)} ROI peaks; pixel a = {a_pix} mm; L = {L} mm; r_perp = {r_perp} mm")
print(f"detector usage: {np.bincount([p.detector_index for p in roi])}")


def centroid(p, delta):
    """Simulator spot centroid (col, row) and omega* for peak p at offset delta; None if lost."""
    g_s = expm(delta) @ R_nom @ (R_nom.T @ (p.g_hkl.double().numpy()))  # g_hkl in ROIPeak is crystal frame
    g = torch.from_numpy(g_s)[None]
    r = get_scattering_omegas_torch(g, torch.norm(g, dim=1), E, 0.0, epsilon=0.0)
    if not bool(r.observable[0]):
        return None
    w = float((r.omega1 if p.omega_branch == 1 else r.omega2)[0])
    _restore_and_rotate(sample, base, w)
    res = _project_peak_on_detector(simulator, sample, detector_list[p.detector_index], verts,
                                    torch.from_numpy(g_s / np.linalg.norm(g_s)).float(), accept)
    _restore_and_rotate(sample, base, 0.0)
    if res is None:
        return None
    row0, col0, _, _ = res
    return np.array([col0, row0]), w


H = 1e-3  # rad; float32 projection needs a step well above ~1e-4 px noise
rows = []
for p in roi:
    # ROIPeak.g_hkl is the crystal-frame vector; sample frame at nominal:
    g_s0 = R_nom @ p.g_hkl.double().numpy()
    w0 = p.nominal_omega
    g = rz(w0) @ g_s0
    gm = np.linalg.norm(g)
    sin_t = gm / (2 * np.linalg.norm(k_vec))
    theta = np.arcsin(sin_t)
    m = np.cross(g, k_vec)
    m /= np.linalg.norm(m)
    ghat = g / gm
    t = np.cross(ghat, m)
    grad_w_lab = -m / (zhat @ m)
    grad_w = rz(w0).T @ grad_w_lab  # sample-frame gradient (delta is applied in the sample frame)
    Ld = L[p.detector_index]
    pred_mag = Ld * np.tan(2 * theta) / (a_pix * np.cos(theta)) * np.linalg.norm(t + (zhat @ t) * grad_w_lab)

    J_u = np.zeros((2, 3))
    ok = True
    for i in range(3):
        e = np.zeros(3)
        e[i] = H
        cp, cm = centroid(p, e), centroid(p, -e)
        if cp is None or cm is None:
            ok = False
            break
        J_u[:, i] = (cp[0] - cm[0]) / (2 * H)
    if not ok:
        continue
    sv = np.linalg.svd(J_u, compute_uv=False)
    rows.append(dict(sin_eta=abs(zhat @ m), theta=theta, det=p.detector_index, grad_w=grad_w, J_u=J_u,
                     fd_mag=sv[0], sv2=sv[1], pred_mag=pred_mag, z_sens=np.linalg.norm(J_u @ zhat)))

P = len(rows)
fd = np.array([r["fd_mag"] for r in rows])
pr = np.array([r["pred_mag"] for r in rows])
sv2 = np.array([r["sv2"] for r in rows])
zs = np.array([r["z_sens"] for r in rows])
th = np.array([r["theta"] for r in rows])
print(f"\n{P} peaks with finite-difference pixel gradients")
print("A. Spot moves along a line (rank-1 pixel gradient): second / first singular value "
      f"median {np.median(sv2 / fd):.2e}, max {np.max(sv2 / fd):.2e}")
rel = np.abs(fd - pr) / pr
print(f"B. |grad s| vs L tan2theta/(a cos theta)|t + (z.t) grad w|: rel err median {np.median(rel):.2e}, "
      f"p90 {np.percentile(rel, 90):.2e}, max {rel.max():.2e}")
print(f"   typical pixel sensitivity: median {np.median(fd):.0f} px/rad = {np.median(fd) * DEG:.2f} px/deg")
print(f"C. Pixel sensitivity to rotation about z: median {np.median(zs):.3f} px/rad "
      f"(parallax scale r_perp/a = {r_perp / a_pix:.1f} px/rad); vs {np.median(fd):.0f} px/rad perpendicular")

# Information matrix: uniform quantisation variance Delta^2/12 per measurement, independent errors.
J_frame = sum(np.outer(r["grad_w"], r["grad_w"]) for r in rows) * 12 / FRAME ** 2
J_pix = sum(r["J_u"].T @ r["J_u"] for r in rows) * 12  # pixel units, variance 1/12 per coordinate
for name, J in (("frames only", J_frame), ("pixels only", J_pix), ("frames + pixels", J_frame + J_pix)):
    try:
        C = np.linalg.inv(J)
        sig = np.sqrt(np.diag(C)) / DEG
        print(f"D. {name:16s} sigma_x,y,z = {sig[0]:.2e}, {sig[1]:.2e}, {sig[2]:.2e} deg   "
              f"RMS misorientation sqrt(tr) = {np.sqrt(np.trace(C)) / DEG:.2e} deg")
    except np.linalg.LinAlgError:
        ev = np.linalg.eigvalsh(J)
        print(f"D. {name:16s} singular: eigenvalues {ev}")
evp = np.linalg.eigh(J_pix)
print(f"   pixel-info eigenvalues {evp[0]}; weakest direction {np.round(evp[1][:, 0], 4)}")

# Closed forms
two_theta_rms = np.sqrt(np.mean((np.tan(2 * th) / np.cos(th)) ** 2))
Lp = np.array([L[r["det"]] for r in rows])
L_rms = np.sqrt(np.mean(Lp ** 2))
sig_z_cf = FRAME / np.sqrt(12 * P)
sig_perp_cf = (a_pix / L_rms) / (two_theta_rms * np.sqrt(6 * P))
print(f"\nE. Closed forms: sigma_z = Delta_w/sqrt(12P) = {sig_z_cf / DEG:.2e} deg; "
      f"sigma_perp = (a/L)/(kappa sqrt(6P)) = {sig_perp_cf / DEG:.2e} deg  "
      f"(kappa_rms = {two_theta_rms:.3f}, L_rms = {L_rms:.3f} mm)")
