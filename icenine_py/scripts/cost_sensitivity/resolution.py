#!/usr/bin/env python3
"""
Phase A2 (finisher/MC study): the centroid-quantisation scale of a voxel's orientation, from the
experiment geometry. This is a scale, NOT a lower bound (see the end).

Model. A voxel with orientation R has M recorded Bragg reflections. Reflection m is read as a
point (frame omega_m, detector column, detector row) = (w_m, u_m, v_m); each coordinate is
quantised to one bin, which we treat as an error with variance (bin)^2 / 12, with bins
  omega: the frame width dw (rad), pixel: 1 px (the detector pixel, 1.48 um), both axes.
A small rotation delta (rotation vector, rad, applied on the left of R) moves reflection m by
  d w_m = g_m . delta,  d (u_m, v_m) = G_m delta,
g_m (3,) and G_m (2,3) from central differences of the batched ray tracer (BatchedObserver) at
+-h. Treating the 3M coordinates as independent Gaussian-equivalent measurements with those
variances, the information matrix on delta is
  J = sum_m [ 12/dw^2 g_m g_m^T + 12 G_m^T G_m ]          (rad^-2)
and the covariance J^-1. Reported: per-axis sigma_i = sqrt((J^-1)_ii) and the 3-D RMS angle
sqrt(tr J^-1) (deg), comparable with a misorientation error (angle of R_hat R^T); the pixel part
(G) and frame part (g) alone too. Per-peak drivers: |sin eta| (frame term), sin theta and the
detector distance (pixel term: the spot moves L * d(angle)).

Why it is a scale and not a bound. Uniform quantisation noise violates the regularity conditions of
the Cramer-Rao bound (the likelihood is not differentiable in the parameter), so an estimator can
beat sqrt(tr J^-1): the exact pixel sets carry sub-pixel edge information, and the finisher's
continuation reaches 0.006 deg, below this scale. Also, even an efficient estimator with
3-D Gaussian errors is below its RMS scale in only about 61% of cases (chi-square, 3 dof),
so individual cases below or above the scale are expected. The equality with
ExactBayes.information_matrix
(asserted) only checks the same algorithm computed two ways; it is not an independent validation.

Usage (from icenine_py/):
  uv run python scripts/cost_sensitivity/resolution.py
"""

import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

HERE = Path(__file__).resolve().parent
ICENINE_PY = HERE.parents[1]

H_DEG = 0.005
OUT_DIR = ICENINE_PY / "benchmarks" / "cost_sensitivity"


def fisher_quantisation(
    grad_omega: np.ndarray, grad_uv: np.ndarray, frame_width_rad: float,
    omega_bin: float = 1.0, pixel_bin: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:  # fmt: skip
    """(J_frame, J_pixel), each (3, 3) in rad^-2, for reflections with frame gradient
    grad_omega (M, 3) [rad per rad] and pixel-coordinate gradient grad_uv (M, 2, 3) [px per rad].
    Uniform quantisation error: variance (bin)^2 / 12 with bin = omega_bin * frame_width_rad for
    the frame and pixel_bin pixels for each detector coordinate."""
    jf = (12.0 / (omega_bin * frame_width_rad) ** 2) * grad_omega.T @ grad_omega
    jp = (12.0 / pixel_bin**2) * np.einsum("mai,maj->ij", grad_uv, grad_uv)
    return jf, jp


def crb(J: np.ndarray) -> Dict[str, Any]:
    """Cramer-Rao summary of an information matrix (rad^-2): per-axis sigma (deg), the 3-D RMS
    angle sqrt(tr J^-1) (deg) and the principal sigmas (deg, ascending). Singular J -> inf."""
    w = np.linalg.eigvalsh(0.5 * (J + J.T))
    if w.min() <= 1e-12 * max(w.max(), 1e-300):
        return dict(sigma_axes_deg=[math.inf] * 3, rms3_deg=math.inf, principal_deg=[math.inf] * 3)
    C = np.linalg.inv(J)
    return dict(
        sigma_axes_deg=[math.degrees(math.sqrt(C[i, i])) for i in range(3)],
        rms3_deg=math.degrees(math.sqrt(float(np.trace(C)))),
        principal_deg=[math.degrees(1.0 / math.sqrt(x)) for x in w[::-1]],
    )


def gradients(obs: Any, h_deg: float = H_DEG) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Central-difference (grad_omega (Mp,3) rad/rad, grad_uv (Mp,2,3) px/rad, mask (M,)) at the
    observer's nominal orientation, over the peaks present at all probe offsets."""
    import torch

    pts = [np.zeros(3)]
    for i in range(3):
        e = np.zeros(3)
        e[i] = h_deg
        pts += [e, -e]
    o = obs.observe(torch.as_tensor(np.array(pts), dtype=obs.dtype))
    present = o.present.all(dim=0)
    om = o.omega[:, present].numpy()
    cent = o.verts[:, present].mean(dim=2).numpy()
    h = math.radians(h_deg)
    gw = np.stack([(om[1 + 2 * i] - om[2 + 2 * i]) / (2 * h) for i in range(3)], axis=-1)
    gu = np.stack([(cent[1 + 2 * i] - cent[2 + 2 * i]) / (2 * h) for i in range(3)], axis=-1)
    return gw, gu, present.numpy()


def voxel_resolution(ctx: Any, vidx: int, min_sin_eta: Optional[float]) -> Dict[str, Any]:
    """Resolution bound of voxel vidx at its true orientation, using the reflection set of
    build_problem(min_sin_eta=...) (None: the net's/sweep's filter ctx.args.min_sin_eta)."""
    import torch
    from generate_toy_orientation_dataset import build_problem
    from icenine.constants import KEV_OVER_HBAR_C_IN_ANG
    from icenine.orientation_eval import BatchedObserver, ExactBayes

    a = ctx.args
    mse = a.min_sin_eta if min_sin_eta is None else min_sin_eta
    pr = build_problem(
        ctx.example_dir, vidx, max_q=a.max_q, detectors=a.detectors, min_sin_eta=mse,
        setup=ctx.setup,
    )  # fmt: skip
    obs = BatchedObserver(
        pr["R_nom"], pr["vertices"], pr["sample"], pr["detector_list"], pr["range_map"],
        pr["exp_setup"], pr["roi_list"],
    )  # fmt: skip
    gw, gu, mask = gradients(obs)
    jf, jp = fisher_quantisation(gw, gu, obs.frame_width_rad)
    # cross-check against the library's linearised information matrix (same formula)
    j_lib = ExactBayes(obs, 1.0, use_pixels=True).information_matrix(np.zeros(3), H_DEG)
    assert np.allclose(jf + jp, j_lib, rtol=1e-6), "A2 information matrix != ExactBayes"
    sin_eta = obs.sin_eta(torch.zeros(1, 3, dtype=obs.dtype))[0].numpy()[mask]
    sin_th = (obs.g_hkl.norm(dim=-1) / (2.0 * KEV_OVER_HBAR_C_IN_ANG * obs.energy)).numpy()[mask]
    det = obs.det_idx.numpy()[mask]
    # lab-frame distance of each detector plane from the origin, per detector index
    dist = np.abs(obs.d_plane.numpy())[mask]
    out = dict(
        vidx=int(vidx), n_peaks=int(mask.sum()), frame_width_deg=math.degrees(obs.frame_width_rad),
        pixel_size=float(obs.d_pw[0]),
        median_abs_sin_eta=float(np.median(sin_eta)) if mask.any() else math.nan,
        median_two_theta_deg=(
            float(2 * np.degrees(np.arcsin(np.median(sin_th)))) if mask.any() else math.nan
        ),
        det_counts=[int((det == k).sum()) for k in range(2)],
        det_dist=[
            float(np.median(dist[det == k])) if (det == k).any() else math.nan for k in range(2)
        ],
        both=crb(jf + jp), frame_only=crb(jf), pixel_only=crb(jp),
    )  # fmt: skip
    return out


def main() -> None:
    sys.path.insert(0, str(HERE.parent / "finisher_diagnosis"))
    sys.path.insert(0, str(HERE.parent / "nn_hybrid"))
    sys.path.insert(0, str(HERE.parent / "common"))
    sys.path.insert(0, str(HERE.parent))
    import findoptimal_sweep as fs
    import run as nnrun
    import diagnose as D

    wargs, _ = nnrun.worker_args(["realistic_s0"])
    fs.init_worker(wargs)
    ctx = fs._W.ctx
    items, _ = D.build_items(HERE / "cache" / "tmp_res", 0)
    vox = sorted({it[0] for it in items})
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    res = []
    for v in vox:
        r_net = voxel_resolution(ctx, v, None)
        r_all = voxel_resolution(ctx, v, 0.0)
        res.append(dict(net_filter=r_net, all_peaks=r_all))
        print(
            v,
            r_net["n_peaks"],
            r_net["both"]["rms3_deg"],
            r_all["n_peaks"],
            r_all["both"]["rms3_deg"],
            flush=True,
        )
    (HERE / "cache" / "resolution.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
