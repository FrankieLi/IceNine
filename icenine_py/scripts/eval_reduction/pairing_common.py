"""Shared loaders of the Task D (spot pairing) scripts: geometry, provenance, truth links.

`load_geometry()` builds the Example2.Simulation geometry once (detectors, omega frames, voxel
centroids and grains, hull points, Cu ring angles) through the same classes as render_full.py and
caches it as a small npz. It checks the frame assumption of ``icenine.spot_pairing`` (lab =
Rz(omega) sample: identity sample-to-lab transform, beam along +x, no beam deflection) and refuses
any other set-up. `load_provenance()` reads the provenance.npz written by
scripts/phase_d/render_provenance.py. `truth_links()` returns the true (detector-0, detector-1)
component pairs of each frame: two components are a true pair when they share a
(voxel, reflection) source in that frame.
"""

import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EX = ROOT / "Examples" / "Example2.ManyGrains"
PHASE_D = EX / "ScatteringData_PhaseD"
BENCH = HERE.parents[1] / "benchmarks" / "eval_reduction" / "D"
CACHE = BENCH / "cache"
sys.path.insert(0, str(HERE.parents[1]))


def load_geometry(seed: int = 0, refresh: bool = False) -> Dict[str, Any]:
    """Geometry dict (numpy arrays and floats); cached in cache/geometry_s<seed>.npz."""
    path = CACHE / f"geometry_s{seed}.npz"
    if path.exists() and not refresh:
        z = np.load(path, allow_pickle=True)
        if "frame_checked" in z.files and bool(z["frame_checked"]):
            return {k: z[k] for k in z.files}
        # an older cache written before the frame check: rebuild it (and check)
    cwd = os.getcwd()
    os.chdir(EX)
    try:
        from icenine.config_file import ConfigFile
        from icenine.constants import KEV_OVER_HBAR_C_IN_ANG
        from icenine.forward_simulation import ForwardSimulation
        from icenine.sample import Sample
        from icenine.simulation import Simulation

        cfg = ConfigFile.from_file(str(EX / "ConfigFiles" / "Example2.Simulation.config"))
        cfg.sample_filename = f"SimInput/rand_500grains_1mm_neworient_s{seed}.mic"
        sim = ForwardSimulation(cfg)
        sim.exp_setup.initialize_experiment()
        oms = sim.exp_setup.get_omega_range_list()
        dets = sim.exp_setup.get_detector_list()
        sim.simulator = Simulation(sim.exp_setup)
        sample = Sample()
        sim.exp_setup.initialize_sample(sample, dets[0])
        check_frame_assumption(cfg, sample)
        vox = sample.get_mic().voxels
        side = float(vox[0].side_length)
        cen = np.zeros((len(vox), 2))
        tri = np.zeros((len(vox), 3, 2))
        ori = np.zeros((len(vox), 9))
        s3 = np.sqrt(3.0)
        for i, v in enumerate(vox):
            x, y = float(v.position[0]), float(v.position[1])
            s = float(v.side_length)
            if v.points_up:
                tv = [[x, y], [x + s, y], [x + s / 2, y + s / 2 * s3]]
            else:
                tv = [[x, y], [x + s / 2, y - s / 2 * s3], [x + s, y]]
            tri[i] = tv
            cen[i] = np.mean(tv, axis=0)
            ori[i] = np.asarray(v.orientation, dtype=np.float64).reshape(-1)
        _, grain = np.unique(np.round(ori, 7), axis=0, return_inverse=True)
        energy = sim.exp_setup.beam_energy
        k = KEV_OVER_HBAR_C_IN_ANG * energy
        # structure list: an empty placeholder phase first, then the crystal (largest list)
        st = max(sample.get_structure_list(), key=lambda t: len(t.get_reflection_vectors()))
        qs = np.array([rv.q_mag for rv in st.get_reflection_vectors()])
        refl_two_theta = 2.0 * np.arcsin(qs / (2.0 * k))
        rings = np.unique(np.round(refl_two_theta, 9))
        from icenine.spot_pairing import geometry_from_detectors

        gs = geometry_from_detectors(dets)
        out: Dict[str, Any] = {
            "omega_lo": np.array([o.low for o in oms]),
            "omega_hi": np.array([o.high for o in oms]),
            "voxel_centroid": cen,
            "voxel_tri": tri,
            "voxel_grain": grain.astype(np.int64),
            "voxel_side": side,
            "refl_two_theta": refl_two_theta,
            "rings_two_theta": rings,
            "registration_mm": float(cfg.detector_dist_deviation),
            "beam_energy_kev": float(energy),
            "num_rows": float(dets[0].num_rows),
            "num_cols": float(dets[0].num_cols),
            "frame_checked": True,
        }
        for di, g in enumerate(gs):
            for f in ("pos", "org", "bj", "bk", "n"):
                out[f"det{di}_{f}"] = getattr(g, f)
            for f in ("pw", "ph", "hw", "hh", "d"):
                out[f"det{di}_{f}"] = float(getattr(g, f))
    finally:
        os.chdir(cwd)
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez(path, **out)
    return out


def check_frame_assumption(cfg: Any, sample: Any) -> None:
    """Raise unless lab = Rz(omega) sample holds: the sample-to-lab transform is the identity
    (zero SampleLocation and SampleOrientation), the beam runs along +x and is not deflected.
    ``icenine.spot_pairing`` and ``pairing_attribution.predicted_peaks`` assume this."""
    m = np.asarray(sample.sample_to_lab_matrix.double().numpy())
    if not np.allclose(m, np.eye(4), rtol=0.0, atol=1e-12):
        raise ValueError(f"sample-to-lab transform is not the identity:\n{m}")
    beam = np.asarray(cfg.beam_direction, dtype=np.float64)
    if not np.allclose(beam / np.linalg.norm(beam), [1.0, 0.0, 0.0], rtol=0.0, atol=1e-12):
        raise ValueError(f"beam direction {beam} is not +x")
    if float(cfg.beam_deflection_chi_laue) != 0.0:
        raise ValueError(f"beam deflection {cfg.beam_deflection_chi_laue} is not 0")


def detector_geoms(geo: Dict[str, Any]) -> List[Any]:
    from icenine.spot_pairing import DetectorGeometry

    return [
        DetectorGeometry(
            pos=geo[f"det{d}_pos"],
            org=geo[f"det{d}_org"],
            bj=geo[f"det{d}_bj"],
            bk=geo[f"det{d}_bk"],
            pw=float(geo[f"det{d}_pw"]),
            ph=float(geo[f"det{d}_ph"]),
            hw=float(geo[f"det{d}_hw"]),
            hh=float(geo[f"det{d}_hh"]),
            n=geo[f"det{d}_n"],
            d=float(geo[f"det{d}_d"]),
        )
        for d in (0, 1)
    ]


def load_provenance(tag: str = "full") -> Dict[str, np.ndarray]:
    z = np.load(PHASE_D / tag / "provenance" / "provenance.npz")
    return {k: z[k] for k in z.files}


def truth_links(prov: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """True pairs: unique (comp0, comp1) over (frame, voxel, refl) peaks seen on both detectors.

    Returns (comp0, comp1, n_peaks): global component ids and the number of (voxel, reflection)
    peaks the pair shares. A peak drawn on one detector only (clipped off the other) gives no link.
    """
    fr, det = prov["tri_frame"], prov["tri_det"]
    key = (
        fr.astype(np.int64) * (1 << 45)
        + prov["tri_voxel"].astype(np.int64) * (1 << 12)
        + prov["tri_refl"].astype(np.int64)
    )
    comp = prov["tri_comp"]
    m0, m1 = det == 0, det == 1
    # a (frame, voxel, refl) peak can have several triangles per detector only if both omegas
    # fall into one frame; keep the first of each
    k0, i0 = np.unique(key[m0], return_index=True)
    k1, i1 = np.unique(key[m1], return_index=True)
    common, a, b = np.intersect1d(k0, k1, return_indices=True)
    c0 = comp[m0][i0[a]]
    c1 = comp[m1][i1[b]]
    pair, n_peaks = np.unique(np.stack([c0, c1], axis=1), axis=0, return_counts=True)
    return pair[:, 0], pair[:, 1], n_peaks


def component_table(prov: Dict[str, np.ndarray], geo: Dict[str, Any]) -> Dict[str, np.ndarray]:
    """Per-component truth: distinct voxels, grains, sample-frame truth centroid (mean of the
    distinct source voxel centroids), mean exact omega, and the detector-border flag."""
    ncomp = len(prov["comp_npix"])
    comp = prov["tri_comp"].astype(np.int64)
    vox = prov["tri_voxel"].astype(np.int64)
    uk = np.unique(comp * (1 << 32) + vox)
    uc, uv = uk >> 32, uk & ((1 << 32) - 1)
    n_vox = np.bincount(uc, minlength=ncomp)
    cen = geo["voxel_centroid"]
    tx = np.bincount(uc, weights=cen[uv, 0], minlength=ncomp) / np.maximum(n_vox, 1)
    ty = np.bincount(uc, weights=cen[uv, 1], minlength=ncomp) / np.maximum(n_vox, 1)
    gk = np.unique(uc * (1 << 32) + geo["voxel_grain"][uv])
    n_grain = np.bincount(gk >> 32, minlength=ncomp)
    om = np.bincount(comp, weights=prov["tri_omega"], minlength=ncomp) / np.maximum(
        np.bincount(comp, minlength=ncomp), 1
    )
    n_rows, n_cols = int(geo["num_rows"]), int(geo["num_cols"])
    edge = (
        (prov["comp_cmin"] == 0)
        | (prov["comp_rmin"] == 0)
        | (prov["comp_cmax"] == n_cols - 1)
        | (prov["comp_rmax"] == n_rows - 1)
    )
    return {
        "n_src": prov["comp_n_src"],
        "n_vox": n_vox,
        "n_grain": n_grain,
        "truth_xy": np.stack([tx, ty], axis=1),
        "omega_mean": om,
        "edge": edge,
        "frame": prov["comp_frame"],
        "det": prov["comp_det"],
        "npix": prov["comp_npix"],
        "c_col": prov["comp_c_col"],
        "c_row": prov["comp_c_row"],
        "i_col": prov["comp_i_col"],
        "i_row": prov["comp_i_row"],
        "cmin": prov["comp_cmin"],
        "cmax": prov["comp_cmax"],
        "rmin": prov["comp_rmin"],
        "rmax": prov["comp_rmax"],
    }


CLASS_NAMES = ["single-source", "one grain, several voxels", "several grains"]


def component_class(ct: Dict[str, np.ndarray]) -> np.ndarray:
    """0 = one (voxel, reflection) source, 1 = several voxels of one grain, 2 = several grains."""
    return np.where(ct["n_src"] == 1, 0, np.where(ct["n_grain"] == 1, 1, 2))


def pair_keys(c0: np.ndarray, c1: np.ndarray, ncomp: int) -> np.ndarray:
    return c0.astype(np.int64) * ncomp + c1.astype(np.int64)


def both_detector_peaks(prov: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """The (frame, voxel, reflection) peaks drawn on both detectors, with their two components
    and the exact omega (the voxel's true, observable peaks)."""
    fr, det = prov["tri_frame"], prov["tri_det"]
    key = (
        fr.astype(np.int64) * (1 << 45)
        + prov["tri_voxel"].astype(np.int64) * (1 << 12)
        + prov["tri_refl"].astype(np.int64)
    )
    m0, m1 = det == 0, det == 1
    k0, i0 = np.unique(key[m0], return_index=True)
    k1, i1 = np.unique(key[m1], return_index=True)
    common, a, b = np.intersect1d(k0, k1, return_indices=True)
    sel0 = np.nonzero(m0)[0][i0[a]]
    sel1 = np.nonzero(m1)[0][i1[b]]
    return {
        "frame": prov["tri_frame"][sel0],
        "voxel": prov["tri_voxel"][sel0],
        "refl": prov["tri_refl"][sel0],
        "omega": prov["tri_omega"][sel0],
        "c0": prov["tri_comp"][sel0],
        "c1": prov["tri_comp"][sel1],
    }
