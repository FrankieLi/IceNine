"""GO5q (I59): the ``strain`` sample plus an in-plane linear strain gradient inside every grain.

  uv run python scripts/grain_opt/make_strain_gradient_sample.py --tag straing

Same orientations, grain map, levels and grain-mean strains as ``make_strained_sample.py`` (same
seed and draw order); then per grain, eps(v) = eps_mean + D d_v with d_v = (xy_v - c) / r (c the
grain centroid, r the RMS voxel distance from it) and D a 6 x 2 Gaussian matrix (off-diagonal
rows with half the variance) scaled so that the gradient part has tensor RMS over the grain's
voxels equal to the grain's level. Only the in-plane gradient is observable in a one-voxel-thick
layer. Writes ``benchmarks/grain_opt/samples/<tag>/`` like ``make_strained_sample.py`` plus
``strain_vox.npy`` (per-voxel [e11 e22 e33 e23 e13 e12]) and per grain ``grad`` (D, 6 x 2),
``centroid`` and ``radius`` in ``sample.json``. Then render with
  uv run python scripts/phase_d/render_full.py --tag <tag> --workers 4 --variants clean --strain \\
      --mic <icenine_py>/benchmarks/grain_opt/samples/<tag>/truth.mic
and build the stores with build_store.py --tag <tag>.
"""

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ICE = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar3"))
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar1"))

import ar3_common as C  # noqa: E402,F401  (sets threads)
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402
from make_strained_sample import LEVELS, draw_strain  # noqa: E402


def tensor_rms(e: np.ndarray) -> np.ndarray:
    e = np.asarray(e, float).reshape(-1, 6)
    return np.sqrt(((e[:, :3] ** 2).sum(1) + 2 * (e[:, 3:] ** 2).sum(1)) / 3.0)


def main() -> None:
    import misses as M
    import torch

    from icenine import grain_fit as GF
    from icenine.mic_file import MicFile

    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="straing")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--mode", default="linear", choices=["linear", "nonlinear"])
    ap.add_argument("--corr", type=float, default=3.0, help="GP correlation length, voxel spacings")
    a = ap.parse_args()
    out_dir = ICE / "benchmarks" / "grain_opt" / "samples" / a.tag
    out_dir.mkdir(parents=True, exist_ok=True)
    mic = MicFile.read(str(G.TRUTH_MIC))
    gmap = np.load(M.GRAINMAP)
    n = len(mic.voxels)
    base = json.loads((G.SAMPLE / "sample.json").read_text())
    gids = sorted(int(g) for g in np.unique(gmap))
    rng = np.random.default_rng(a.seed)
    lv_idx = rng.integers(0, len(LEVELS), size=len(gids))  # the strain sample's draws, in order
    eps = {g: draw_strain(rng, LEVELS[i]) for g, i in zip(gids, lv_idx)}
    level = {g: float(LEVELS[i]) for g, i in zip(gids, lv_idx)}
    xy = np.array([v.position[:2] for v in mic.voxels], float)
    eps_vox = np.zeros((n, 6))
    grad, cen, rad, lin_share = {}, {}, {}, {}
    from scipy.spatial import cKDTree

    uxy = np.unique(np.round(xy.astype(np.float64), 6), axis=0)  # triangle pairs share positions
    spacing = float(np.median(cKDTree(uxy).query(uxy, k=2)[0][:, 1]))

    def scaled(part: np.ndarray, target: float) -> np.ndarray:
        s_ = float(np.sqrt((tensor_rms(part) ** 2).mean())) if len(part) > 1 else 0.0
        return part * (target / s_) if s_ > 0 else part * 0.0
    for g in gids:
        m = np.nonzero(gmap == g)[0]
        c = xy[m].mean(0)
        r = float(np.sqrt(((xy[m] - c) ** 2).sum(1).mean())) or 1.0
        d = (xy[m] - c) / r
        D = rng.standard_normal((6, 2))
        D[3:] *= np.sqrt(0.5)
        part = d @ D.T  # (k, 6)
        s = float(np.sqrt((tensor_rms(part) ** 2).mean())) if len(m) > 1 else 1.0
        if a.mode == "linear":
            D *= level[g] / s if s > 0 else 0.0
            intra = d @ D.T
        else:  # GO5s: linear + quadratic + smooth Gaussian random field, equal RMS level / sqrt(3)
            tgt = level[g] / np.sqrt(3.0)
            D *= tgt / s if s > 0 else 0.0
            lin = d @ D.T
            f2 = np.stack([d[:, 0] ** 2, d[:, 0] * d[:, 1], d[:, 1] ** 2], 1)
            f2 = f2 - f2.mean(0)
            Q = rng.standard_normal((6, 3))
            Q[3:] *= np.sqrt(0.5)
            A1 = np.c_[np.ones(len(m)), d]  # amendment: parts orthogonal to the linear basis

            def orth(x_: np.ndarray) -> np.ndarray:
                c_, *_ = np.linalg.lstsq(A1, x_, rcond=None)
                return x_ - A1 @ c_

            quad = scaled(orth(f2 @ Q.T), tgt)
            dist2 = ((xy[m][:, None, :] - xy[m][None, :, :]) ** 2).sum(-1)
            K = np.exp(-0.5 * dist2 / (a.corr * spacing) ** 2)
            if not np.isfinite(K).all():
                raise SystemExit(f"non-finite kernel: grain {g}, {len(m)} voxels, spacing {spacing}")
            ev, U = np.linalg.eigh(K)  # PSD square root (the kernel is near-singular)
            Lk = U * np.sqrt(np.clip(ev, 0.0, None))
            z = Lk @ rng.standard_normal((len(m), 6))
            z[:, 3:] *= np.sqrt(0.5)
            z = z - z.mean(0)
            gp = scaled(orth(z), tgt)
            intra = lin + quad + gp
        eps_vox[m] = eps[g] + intra
        if len(m) > 3:  # oracle: least-squares linear projection of the intragranular field
            A = np.c_[np.ones(len(m)), d]
            coef, *_ = np.linalg.lstsq(A, intra, rcond=None)
            proj = A @ coef - (A @ coef).mean(0)
            tot = float((tensor_rms(intra - intra.mean(0)) ** 2).sum())
            lin_share[g] = float((tensor_rms(proj) ** 2).sum()) / tot if tot > 0 else 1.0
        grad[g], cen[g], rad[g] = D, c, r
    Mv = GF.strain_to_m(torch.as_tensor(eps_vox, dtype=torch.float64)).numpy()
    for v, mm in zip(mic.voxels, Mv):
        v.deformation = mm
    mic.write(str(out_dir / "truth.mic"), deformation_digits=9)
    shutil.copyfile(G.GRID_MIC, out_dir / "grid.mic")
    back = MicFile.read(str(out_dir / "truth.mic"))
    eps_back = GF.m_to_strain(np.stack([v.deformation for v in back.voxels]).astype(np.float64))
    sha = hashlib.sha256((out_dir / "truth.mic").read_bytes()).hexdigest()
    info = dict(base)
    info.update(
        tag=a.tag, strain_seed=a.seed, strain_levels=list(LEVELS), truth_sha256=sha,
        strain_level_meaning="grain-mean tensor RMS; gradient part has the same tensor RMS over "
        "the grain's voxels; eps(v) = eps_mean + D (xy_v - c) / r",
        max_eps_roundtrip_error=float(np.abs(eps_back - eps_vox).max()),
        mode=a.mode,
        gp_corr_spacings=a.corr if a.mode == "nonlinear" else None,
        oracle_linear_share_median=float(np.median(list(lin_share.values()))) if lin_share else None,
    )
    info["grains"] = [
        dict(g, eps=eps[int(g["id"])].tolist(), strain_level=level[int(g["id"])],
             grad=grad[int(g["id"])].tolist(), centroid=cen[int(g["id"])].tolist(),
             radius=rad[int(g["id"])])
        for g in base["grains"]
    ]
    (out_dir / "sample.json").write_text(json.dumps(info) + "\n")
    np.save(out_dir / "grainmap.npy", gmap)
    np.save(out_dir / "strain_vox.npy", eps_vox)
    print(json.dumps(dict(tag=a.tag, n_voxels=n, n_grains=len(gids), sha256=sha[:12],
                          roundtrip=info["max_eps_roundtrip_error"],
                          vox_strain_rms_med=float(np.median(tensor_rms(eps_vox))),
                          oracle_linear_share_median=info["oracle_linear_share_median"])))


if __name__ == "__main__":
    main()
