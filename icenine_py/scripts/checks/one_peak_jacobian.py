import copy, sys
from pathlib import Path
import numpy as np, torch
ICENINE_PY = Path("/Users/sfli/Research/IceNine/icenine_py")
sys.path.insert(0, str(ICENINE_PY / "scripts"))
from generate_toy_orientation_dataset import setup_example
from icenine.constants import KEV_OVER_HBAR_C_IN_ANG
from icenine.diffraction_core import get_scattering_omegas_torch
from icenine.orientation_nn import _project_peak_on_detector, _restore_and_rotate, define_roi_set
from icenine.peak_filters import TrivialAcceptFn
DEG = np.pi / 180
def expm(v):
    th = np.linalg.norm(v)
    if th == 0: return np.eye(3)
    k = v / th; K = np.array([[0,-k[2],k[1]],[k[2],0,-k[0]],[-k[1],k[0],0]])
    return np.eye(3) + np.sin(th)*K + (1-np.cos(th))*K@K
mic, sample, dl, rm, es, sim, sl, gv = setup_example(ICENINE_PY.parent / "Examples" / "Example2.ThreeVoxels")
E = float(es.beam_energy); base = sample.sample_to_lab_matrix[:3,:3].clone(); acc = TrivialAcceptFn()
k = np.array([KEV_OVER_HBAR_C_IN_ANG*E,0,0])
for r_perp in (0.012, 0.5):
    vox = copy.deepcopy(mic.voxels[0]); phi = np.arctan2(vox.position[1], vox.position[0])
    vox.position = np.array([r_perp*np.cos(phi), r_perp*np.sin(phi), vox.position[2]])
    R = vox.orientation.astype(np.float64); verts = gv(vox)
    roi = define_roi_set(torch.from_numpy(R).float(), verts, sample, dl, rm, es, sl, sim, phase_index=vox.phase)
    # pick the same reflection/branch both times: moderate eta, detector 0
    def key(p): return (p.reflection_index, p.omega_branch)
    if r_perp == 0.012:
        cands = []
        for p in roi:
            g = R @ p.g_hkl.double().numpy(); gl = np.array([[np.cos(p.nominal_omega),-np.sin(p.nominal_omega),0],[np.sin(p.nominal_omega),np.cos(p.nominal_omega),0],[0,0,1]]) @ g
            m = np.cross(gl, k); m /= np.linalg.norm(m)
            cands.append((abs(abs(m[2]) - 0.5), p))
        target = key(min(cands, key=lambda c: c[0])[1])
    p = next(q for q in roi if key(q) == target)
    def spot(delta):
        g_s = expm(delta) @ R @ p.g_hkl.double().numpy(); g = torch.from_numpy(g_s)[None]
        r = get_scattering_omegas_torch(g, torch.norm(g, dim=1), E, 0.0, epsilon=0.0)
        w = float((r.omega1 if p.omega_branch == 1 else r.omega2)[0])
        _restore_and_rotate(sample, base, w)
        res = _project_peak_on_detector(sim, sample, dl[p.detector_index], verts, torch.from_numpy(g_s/np.linalg.norm(g_s)).float(), acc)
        _restore_and_rotate(sample, base, 0.0)
        return np.array([res[1], res[0]]), w
    H = 1e-3; Ju = np.zeros((2,3)); dw = np.zeros(3)
    for i in range(3):
        e = np.zeros(3); e[i] = H
        (up, wp), (um, wm) = spot(e), spot(-e)
        Ju[:, i] = (up - um)/(2*H); dw[i] = (wp - wm)/(2*H)
    (u0, w0) = spot(np.zeros(3))
    U, S, Vt = np.linalg.svd(Ju)
    phiv = np.degrees(np.arctan2(*(np.array([[np.cos(w0),-np.sin(w0)],[np.sin(w0),np.cos(w0)]]) @ vox.position[:2])[::-1]))
    print(f"\nr_perp = {r_perp*1000:.0f} um | reflection {p.reflection_index} branch {p.omega_branch} | omega* = {np.degrees(w0):.2f} deg | spot (col,row) = ({u0[0]:.1f}, {u0[1]:.1f}) | voxel lab azimuth {phiv:.0f} deg")
    print("  Jacobian (px per rad), rows = (col, row), columns = rotation about sample x, y, z:")
    for row, name in zip(Ju, ("col", "row")): print(f"    {name}: " + "  ".join(f"{v:9.1f}" for v in row))
    print(f"  per degree:  col " + "  ".join(f"{v*DEG:7.2f}" for v in Ju[0]) + "   row " + "  ".join(f"{v*DEG:7.2f}" for v in Ju[1]))
    print(f"  d omega*/d delta = {np.round(dw, 3)}   (z-component should be -1)")
    print(f"  singular values {np.round(S,1)} px/rad; ratio {S[1]/S[0]:.3f}")
    print(f"  strongest detector direction {np.round(U[:,0],3)} <- from rotation axis {np.round(Vt[0],3)}")
    print(f"  weakest  detector direction {np.round(U[:,1],3)} <- from rotation axis {np.round(Vt[1],3)}")
    print(f"  null (blind) rotation axis {np.round(Vt[2],3)}; g_s direction {np.round((R@p.g_hkl.double().numpy())/np.linalg.norm(R@p.g_hkl.double().numpy()),3)}")
