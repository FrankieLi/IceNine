"""Check the 'frame boundaries are parallel planes in delta-space' claim (nn_inverse_problem_formulation.md §3.1).

delta is a rotation vector applied in the sample frame: R(delta) = Exp([delta]x) R_nom, so g_s(delta) = Exp(delta) g_s0.
"""

import sys
from pathlib import Path

import numpy as np
import torch

ICENINE_PY = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ICENINE_PY / "scripts"))
from generate_toy_orientation_dataset import setup_example  # noqa: E402

from icenine.constants import KEV_OVER_HBAR_C_IN_ANG  # noqa: E402
from icenine.diffraction_core import get_scattering_omegas_torch  # noqa: E402

DEG = np.pi / 180
FRAME = 1.0 * DEG


def rz(w):
    c, s = np.cos(w), np.sin(w)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def expm(v):
    """Rodrigues: rotation matrices for rotation vectors v (N,3)."""
    v = np.atleast_2d(v)
    th = np.linalg.norm(v, axis=1)
    out = np.tile(np.eye(3), (len(v), 1, 1))
    nz = th > 0
    k = v[nz] / th[nz, None]
    K = np.zeros((nz.sum(), 3, 3))
    K[:, 0, 1], K[:, 0, 2] = -k[:, 2], k[:, 1]
    K[:, 1, 0], K[:, 1, 2] = k[:, 2], -k[:, 0]
    K[:, 2, 0], K[:, 2, 1] = -k[:, 1], k[:, 0]
    s, c = np.sin(th[nz])[:, None, None], np.cos(th[nz])[:, None, None]
    out[nz] = np.eye(3) + s * K + (1 - c) * K @ K
    return out


def omega_star(g_batch, E, branch):
    g = torch.from_numpy(g_batch)
    r = get_scattering_omegas_torch(g, torch.norm(g, dim=1), float(E), 0.0, epsilon=0.0)
    w = (r.omega1 if branch == 1 else r.omega2).numpy()
    return w, r.observable.numpy()


def wrap(d):
    return (d + np.pi) % (2 * np.pi) - np.pi


example_dir = ICENINE_PY.parent / "Examples" / "Example2.ThreeVoxels"
mic, sample, detector_list, range_map, exp_setup, simulator, structure_list, get_vertices = (
    setup_example(example_dir)
)
E = float(exp_setup.beam_energy)
k_vec = np.array([KEV_OVER_HBAR_C_IN_ANG * E, 0.0, 0.0])
refls = structure_list[mic.voxels[0].phase].get_reflection_vectors()
g_hkl = np.stack([r.q_vec for r in refls]).astype(np.float64)

rng = np.random.default_rng(0)
N_DIR = 200
dirs = rng.normal(size=(N_DIR, 3))
dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
RADII_DEG = [0.1, 0.5, 1.0, 2.0]
N_GRAD_DIR = 12  # directions at which to re-estimate the gradient (plane-normal tilt)
FD = 1e-7

records = []
z_shift_err = []
for voxel in mic.voxels:
    g_s0 = (voxel.orientation.astype(np.float64) @ g_hkl.T).T
    w1, w2 = [omega_star(g_s0, E, b)[0] for b in (1, 2)]
    obs0 = omega_star(g_s0, E, 1)[1]
    for i in np.where(obs0)[0]:
        g = g_s0[i]
        gm = np.linalg.norm(g)
        sin_t = gm / (2 * np.linalg.norm(k_vec))
        chi = np.arccos(np.clip(g[2] / gm, -1, 1))
        sin_eta = np.sqrt(max(np.sin(chi) ** 2 - sin_t**2, 0.0)) / np.sqrt(1 - sin_t**2)
        for branch, w0 in ((1, w1[i]), (2, w2[i])):
            g_lab = rz(w0) @ g
            fprime = k_vec @ np.cross([0.0, 0.0, 1.0], g_lab)
            grad = (
                -(rz(w0).T @ np.cross(g_lab, k_vec)) / fprime
            )  # analytic d omega*/d delta at delta=0

            # Exact z property: rotation about the sample rotation axis shifts omega* by exactly -beta.
            for beta in (0.01 * DEG, 1.0 * DEG, 10.0 * DEG):
                wz, _ = omega_star((expm(np.array([[0, 0, beta]]))[0] @ g)[None], E, branch)
                z_shift_err.append(abs(wrap(wz[0] - w0) + beta))

            rec = dict(
                sin_eta=sin_eta,
                tan_theta=sin_t / np.sqrt(1 - sin_t**2),
                grad_norm=np.linalg.norm(grad),
                grad_z=grad[2],
            )
            for r_deg in RADII_DEG:
                d = dirs * (r_deg * DEG)
                gp = np.einsum("nij,j->ni", expm(d), g)
                wp, obs = omega_star(gp, E, branch)
                lin = w0 + d @ grad
                err = np.abs(wrap(wp - lin))[obs]
                rec[f"lost_{r_deg}"] = (~obs).mean()
                rec[f"linerr_{r_deg}"] = err.max() / FRAME if err.size else np.nan

                # Plane-normal tilt: gradient at perturbed points vs at delta=0 (central differences).
                tilts, norm_ratio = [], []
                for dd in d[:N_GRAD_DIR]:
                    gr = np.zeros(3)
                    ok = True
                    for a in range(3):
                        e = np.zeros(3)
                        e[a] = FD
                        wpp, o1 = omega_star((expm((dd + e)[None])[0] @ g)[None], E, branch)
                        wmm, o2 = omega_star((expm((dd - e)[None])[0] @ g)[None], E, branch)
                        ok &= bool(o1[0] and o2[0])
                        gr[a] = wrap(wpp[0] - wmm[0]) / (2 * FD)
                    if ok:
                        cosang = gr @ grad / (np.linalg.norm(gr) * np.linalg.norm(grad))
                        tilts.append(np.degrees(np.arccos(np.clip(cosang, -1, 1))))
                        norm_ratio.append(np.linalg.norm(gr) / np.linalg.norm(grad))
                rec[f"tilt_{r_deg}"] = max(tilts) if tilts else np.nan
                rec[f"spacing_{r_deg}"] = max(abs(np.log(norm_ratio))) if norm_ratio else np.nan
            records.append(rec)

R = {k: np.array([r[k] for r in records]) for k in records[0]}
print(f"{len(records)} (reflection, branch) peaks over {len(mic.voxels)} voxels")

print("\n1. Analytic gradient checks")
print(f"   | |grad| * |sin eta| - 1 |   max {np.abs(R['grad_norm'] * R['sin_eta'] - 1).max():.2e}")
print(
    f"   | grad_z + 1 |               max {np.abs(R['grad_z'] + 1).max():.2e}   (every peak: d omega*/d delta_z = -1)"
)
print(
    f"   exact z-shift: | omega*(beta z) - omega*(0) + beta |, beta up to 10 deg: max {max(z_shift_err):.2e} rad"
)

bins = [(0.0, 0.1), (0.1, 0.3), (0.3, 1.01)]
print(
    "\n2. Linearisation error: max over 200 directions on the sphere |delta| = r, in frames (1 frame = 1 deg)."
)
print(
    "   Per |sin eta| bin: median / max over peaks.  'lost' = fraction of directions where the peak stops diffracting."
)
hdr = "   r (deg) | " + " | ".join(
    f"|sin eta| in [{a:.1f},{min(b,1):.1f}) n={((R['sin_eta']>=a)&(R['sin_eta']<b)).sum():4d}"
    for a, b in bins
)
print(hdr)
for r_deg in RADII_DEG:
    cells = []
    for a, b in bins:
        m = (R["sin_eta"] >= a) & (R["sin_eta"] < b)
        e = R[f"linerr_{r_deg}"][m]
        lost = R[f"lost_{r_deg}"][m].mean()
        cells.append(f"{np.nanmedian(e):8.2e} / {np.nanmax(e):8.2e}  lost {lost:5.1%}")
    print(f"   {r_deg:7.1f} | " + " | ".join(cells))

print(
    "\n3. Parallelism: max tilt of the plane normal (deg) and max |log(|grad| ratio)| (spacing change),"
)
print(
    "   comparing the gradient at 12 points on |delta| = r with the gradient at delta = 0. median / max over peaks."
)
for r_deg in RADII_DEG:
    cells = []
    for a, b in bins:
        m = (R["sin_eta"] >= a) & (R["sin_eta"] < b)
        t, s = R[f"tilt_{r_deg}"][m], R[f"spacing_{r_deg}"][m]
        cells.append(
            f"tilt {np.nanmedian(t):7.3f}/{np.nanmax(t):7.3f}  dspacing {np.nanmedian(s):6.3f}/{np.nanmax(s):6.3f}"
        )
    print(f"   {r_deg:7.1f} | " + " | ".join(cells))

# Heuristic scale for where linearisation fails: curvature/gradient ~ tan(theta)/sin^2(eta)
print(
    "\n4. Heuristic check: linearisation error (frames) at r = 1 deg vs  r^2 * tan(theta) / |sin eta|^3  (log-log fit)"
)
x = (1.0 * DEG) ** 2 * R["tan_theta"] / R["sin_eta"] ** 3
y = R["linerr_1.0"] * FRAME
ok = np.isfinite(y) & (y > 0)
slope, icpt = np.polyfit(np.log(x[ok]), np.log(y[ok]), 1)
print(f"   slope {slope:.3f}, prefactor {np.exp(icpt):.3f}  (slope 1 => the scaling holds)")
