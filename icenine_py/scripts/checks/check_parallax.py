"""How does voxel distance from the rotation axis (r_perp) change the spot-motion and resolution results?

Re-runs the check_resolution analysis for Example2 voxel 0 moved to several r_perp values.
Parallax prediction: a rotation beta about z shifts omega* by -beta, which moves the voxel's lab position by
beta * (z x x_v) with |z x x_v| = r_perp. The spot then shifts by that displacement's component perpendicular
to the diffracted ray, projected onto the detector: roughly r_perp / a px per rad.
Crossover: parallax pixel information about the z component beats frame information when r_perp > a / Delta_omega_f.
"""

import copy
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
H = 1e-3


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
base = sample.sample_to_lab_matrix[:3, :3].clone()
accept = TrivialAcceptFn()
a_pix = detector_list[0].pixel_width
zhat = np.array([0.0, 0.0, 1.0])

rng = np.random.default_rng(1)
MAX_PEAKS = 250  # subsample for speed

for r_perp in (0.012, 0.05, 0.1, 0.25, 0.5):
    voxel = copy.deepcopy(mic.voxels[0])
    phi = np.arctan2(voxel.position[1], voxel.position[0])
    voxel.position = np.array([r_perp * np.cos(phi), r_perp * np.sin(phi), voxel.position[2]])
    R_nom = voxel.orientation.astype(np.float64)
    verts = get_vertices(voxel)
    roi = define_roi_set(torch.from_numpy(R_nom).float(), verts, sample, detector_list, range_map,
                         exp_setup, structure_list, simulator, phase_index=voxel.phase)
    if len(roi) > MAX_PEAKS:
        roi = [roi[i] for i in rng.choice(len(roi), MAX_PEAKS, replace=False)]

    def centroid(p, delta):
        g_s = expm(delta) @ R_nom @ p.g_hkl.double().numpy()
        g = torch.from_numpy(g_s)[None]
        r = get_scattering_omegas_torch(g, torch.norm(g, dim=1), E, 0.0, epsilon=0.0)
        if not bool(r.observable[0]):
            return None
        w = float((r.omega1 if p.omega_branch == 1 else r.omega2)[0])
        _restore_and_rotate(sample, base, w)
        res = _project_peak_on_detector(simulator, sample, detector_list[p.detector_index], verts,
                                        torch.from_numpy(g_s / np.linalg.norm(g_s)).float(), accept)
        _restore_and_rotate(sample, base, 0.0)
        return None if res is None else np.array([res[1], res[0]])

    ratios, zsens, perp = [], [], []
    J_frame = np.zeros((3, 3))
    J_pix = np.zeros((3, 3))
    for p in roi:
        Ju = np.zeros((2, 3))
        ok = True
        for i in range(3):
            e = np.zeros(3)
            e[i] = H
            cp, cm = centroid(p, e), centroid(p, -e)
            if cp is None or cm is None:
                ok = False
                break
            Ju[:, i] = (cp - cm) / (2 * H)
        if not ok:
            continue
        sv = np.linalg.svd(Ju, compute_uv=False)
        ratios.append(sv[1] / sv[0])
        zsens.append(np.linalg.norm(Ju @ zhat))
        perp.append(sv[0])
        # frame gradient in the sample frame, by finite differences of omega*
        gw = np.zeros(3)
        for i in range(3):
            e = np.zeros(3)
            e[i] = 1e-7
            ws = []
            for sgn in (1, -1):
                g_s = expm(sgn * e) @ R_nom @ p.g_hkl.double().numpy()
                g = torch.from_numpy(g_s)[None]
                r = get_scattering_omegas_torch(g, torch.norm(g, dim=1), E, 0.0, epsilon=0.0)
                ws.append(float((r.omega1 if p.omega_branch == 1 else r.omega2)[0]))
            gw[i] = (ws[0] - ws[1]) / 2e-7
        J_frame += np.outer(gw, gw) * 12 / FRAME ** 2
        J_pix += Ju.T @ Ju * 12
    P = len(ratios)
    C = np.linalg.inv(J_frame + J_pix)
    Cz_frames = np.linalg.inv(J_frame)[2, 2]
    sig = np.sqrt(np.diag(C)) / DEG
    print(f"r_perp = {r_perp * 1000:5.0f} um  (a/Delta_w = {a_pix / FRAME * 1000:.0f} um)  P = {P}")
    print(f"   rank-one ratio sv2/sv1: median {np.median(ratios):.3f}, max {np.max(ratios):.3f}")
    print(f"   pixel sensitivity to z-rotation: median {np.median(zsens):.1f} px/rad "
          f"(r_perp/a = {r_perp / a_pix:.1f}); perpendicular: median {np.median(perp):.0f} px/rad")
    print(f"   sigma_x,y,z = {sig[0]:.1e}, {sig[1]:.1e}, {sig[2]:.1e} deg;  frames-only sigma_z = "
          f"{np.sqrt(Cz_frames) / DEG:.1e} deg;  Delta_w/sqrt(12P) = {FRAME / np.sqrt(12 * P) / DEG:.1e} deg")
