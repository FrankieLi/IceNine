"""Check the closed-form spot-motion Jacobian against simulator finite differences.

Gamma_L = (1/a) B Q [ -(L/|k|) [g]x (I + z grad_w^T) + (z x x_v) grad_w^T ],   Gamma = Gamma_L Rz(w*)
B: rows = detector column/row axes (lab), Q = I - khat' n^T/(n.khat'), L = ray length voxel->detector,
grad_w = -c/(z.c), c = (g x k)/|g x k|, x_v = voxel lab position at w*.
"""

import copy, sys
from pathlib import Path
import numpy as np, torch

ICENINE_PY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ICENINE_PY / "scripts"))
from generate_toy_orientation_dataset import setup_example
from icenine.constants import KEV_OVER_HBAR_C_IN_ANG
from icenine.diffraction_core import get_scattering_omegas_torch
from icenine.orientation_nn import _project_peak_on_detector, _restore_and_rotate, define_roi_set
from icenine.peak_filters import TrivialAcceptFn


def expm(v):
    th = np.linalg.norm(v)
    if th == 0:
        return np.eye(3)
    k = v / th
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


def rz(w):
    c, s = np.cos(w), np.sin(w)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def skew(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


mic, sample, dl, rm, es, sim, sl, gv = setup_example(
    ICENINE_PY.parent / "Examples" / "Example2.ThreeVoxels"
)
E = float(es.beam_energy)
base = sample.sample_to_lab_matrix[:3, :3].clone()
acc = TrivialAcceptFn()
print("sample translation:", sample.sample_to_lab_matrix[:3, 3].numpy())
k = np.array([KEV_OVER_HBAR_C_IN_ANG * E, 0, 0])
kmag = np.linalg.norm(k)
zhat = np.array([0, 0, 1.0])
target = (362, 2)
rows = []
for r_perp in (0.012, 0.1, 0.25, 0.5):
    vox = copy.deepcopy(mic.voxels[0])
    phi = np.arctan2(vox.position[1], vox.position[0])
    vox.position = np.array([r_perp * np.cos(phi), r_perp * np.sin(phi), vox.position[2]])
    R = vox.orientation.astype(np.float64)
    verts = gv(vox)
    roi = define_roi_set(
        torch.from_numpy(R).float(), verts, sample, dl, rm, es, sl, sim, phase_index=vox.phase
    )
    for p in roi:
        det = dl[p.detector_index]

        def spot(delta):
            g_s = expm(delta) @ R @ p.g_hkl.double().numpy()
            g = torch.from_numpy(g_s)[None]
            r = get_scattering_omegas_torch(g, torch.norm(g, dim=1), E, 0.0, epsilon=0.0)
            if not bool(r.observable[0]):
                return None
            w = float((r.omega1 if p.omega_branch == 1 else r.omega2)[0])
            _restore_and_rotate(sample, base, w)
            res = _project_peak_on_detector(
                sim, sample, det, verts, torch.from_numpy(g_s / np.linalg.norm(g_s)).float(), acc
            )
            _restore_and_rotate(sample, base, 0.0)
            return None if res is None else np.array([res[1], res[0]])

        H = 1e-3
        G_fd = np.zeros((2, 3))
        ok = True
        for i in range(3):
            e = np.zeros(3)
            e[i] = H
            up, um = spot(e), spot(-e)
            if up is None or um is None:
                ok = False
                break
            G_fd[:, i] = (up - um) / (2 * H)
        if not ok:
            continue
        # closed form
        w = p.nominal_omega
        g = rz(w) @ R @ p.g_hkl.double().numpy()
        c = np.cross(g, k)
        c /= np.linalg.norm(c)
        gw = -c / (zhat @ c)
        kp = (k + g) / kmag
        x_v = rz(w) @ verts.double().numpy().mean(axis=0)
        n = det.detector_plane.normal.double().numpy()
        dplane = float(det.detector_plane.d)
        L = -(n @ x_v + dplane) / (n @ kp)
        Q = np.eye(3) - np.outer(kp, n) / (n @ kp)
        Bm = np.stack(
            [
                det._lab_frame_basis_j.double().numpy() / det.pixel_width,
                det._lab_frame_basis_k.double().numpy() / det.pixel_height,
            ]
        )
        ring = -(L / kmag) * skew(g) @ (np.eye(3) + np.outer(zhat, gw))
        par = np.outer(np.cross(zhat, x_v), gw)
        G_cf = Bm @ Q @ (ring + par) @ rz(w)
        G_ring = Bm @ Q @ ring @ rz(w)
        G_par = Bm @ Q @ par @ rz(w)
        err = np.linalg.norm(G_cf - G_fd) / np.linalg.norm(G_fd)
        rows.append((r_perp, err))
        if (p.reflection_index, p.omega_branch) == target:
            np.set_printoptions(precision=1, suppress=True)
            print(
                f"\nr_perp={r_perp*1000:.0f} um, example peak: FD\n{G_fd}\nclosed form\n{G_cf}\n  ring part\n{G_ring}\n  parallax part\n{G_par}\n  rel err {err:.3e}; blind axis check |Gamma g_s|/|Gamma| = {np.linalg.norm(G_cf @ (R @ p.g_hkl.double().numpy())) / np.linalg.norm(G_cf) / np.linalg.norm(p.g_hkl.double().numpy()):.1e}"
            )
rows = np.array(rows)
for r_perp in np.unique(rows[:, 0]):
    e = rows[rows[:, 0] == r_perp, 1]
    print(
        f"r_perp={r_perp*1000:4.0f} um: {len(e)} peaks, rel Frobenius error median {np.median(e):.2e}, p90 {np.percentile(e,90):.2e}, max {e.max():.2e}"
    )
