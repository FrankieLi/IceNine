"""GO5: the ``grad`` ground truth plus a constant elastic strain per grain, as a rendered sample.

  uv run python scripts/grain_opt/make_strained_sample.py --tag strain

Reads the AR3 gradient sample (``gradient_truth.mic``, its grain map and ``sample.json``: the
orientations are unchanged) and gives every grain one crystal-frame symmetric strain tensor
eps (3x3). Levels: one draw per grain, uniformly from LEVELS (1e-4, 5e-4, 1e-3), generator seed
``--seed``. The level is the tensor RMS sqrt(trace(eps^2) / 3) (a pure dilation d I has level d);
the direction is a Gaussian symmetric tensor (off-diagonals with half the variance of the
diagonal) scaled to that level. Convention (``icenine.grain_fit.strain_to_m``): M = (I + eps)^-1,
g = O M g0. Writes ``benchmarks/grain_opt/samples/<tag>/``: ``truth.mic`` (19 columns, columns
13-18 = M00 M11 M22 M01 M12 M02 written with 9 digits after the point), ``grid.mic`` (the
gradient grid), ``sample.json`` (per grain: ``eps`` = [e11 e22 e33 e23 e13 e12], ``strain_level``,
the rotation fields of the grad sample), ``grainmap.npy``. Then render with the strain on
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
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar3"))
sys.path.insert(0, str(ICE / "scripts" / "autoresearch" / "ar1"))

import ar3_common as C  # noqa: E402,F401  (sets threads)
import gradient_common as G  # noqa: E402
import numpy as np  # noqa: E402

LEVELS = (1e-4, 5e-4, 1e-3)  # tensor RMS, sqrt(trace(eps^2) / 3)


def draw_strain(rng: np.random.Generator, level: float) -> np.ndarray:
    """Crystal-frame symmetric strain [e11 e22 e33 e23 e13 e12] of tensor RMS ``level``."""
    a = rng.standard_normal(6)
    a[3:] *= np.sqrt(0.5)
    e = np.array([[a[0], a[5], a[4]], [a[5], a[1], a[3]], [a[4], a[3], a[2]]])
    e *= level / np.sqrt(np.trace(e @ e) / 3.0)
    return np.array([e[0, 0], e[1, 1], e[2, 2], e[1, 2], e[0, 2], e[0, 1]])


def main() -> None:
    import misses as M
    import torch

    from icenine import grain_fit as GF
    from icenine.mic_file import MicFile

    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="strain")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    out_dir = ICE / "benchmarks" / "grain_opt" / "samples" / a.tag
    out_dir.mkdir(parents=True, exist_ok=True)
    mic = MicFile.read(str(G.TRUTH_MIC))
    gmap = np.load(M.GRAINMAP)
    n = len(mic.voxels)
    assert len(gmap) == n
    base = json.loads((G.SAMPLE / "sample.json").read_text())
    gids = sorted(int(g) for g in np.unique(gmap))
    assert gids == sorted(int(g["id"]) for g in base["grains"])
    rng = np.random.default_rng(a.seed)
    lv_idx = rng.integers(0, len(LEVELS), size=len(gids))
    eps = {g: draw_strain(rng, LEVELS[i]) for g, i in zip(gids, lv_idx)}
    level = {g: float(LEVELS[i]) for g, i in zip(gids, lv_idx)}
    eps_vox = np.stack([eps[int(g)] for g in gmap])
    Mv = GF.strain_to_m(torch.as_tensor(eps_vox, dtype=torch.float64)).numpy()
    for v, m in zip(mic.voxels, Mv):
        v.deformation = m  # float64: the writer rounds to 9 digits
    mic.write(str(out_dir / "truth.mic"), deformation_digits=9)
    shutil.copyfile(G.GRID_MIC, out_dir / "grid.mic")

    back = MicFile.read(str(out_dir / "truth.mic"))
    Mb = np.stack([v.deformation for v in back.voxels]).astype(np.float64)
    eps_back = GF.m_to_strain(Mb)
    gxy = np.array([v.position[:2] for v in MicFile.read(str(out_dir / "grid.mic")).voxels])
    assert np.abs(gxy - np.array([v.position[:2] for v in back.voxels])).max() < 1e-6
    sha = hashlib.sha256((out_dir / "truth.mic").read_bytes()).hexdigest()
    info = dict(base)
    info.update(
        tag=a.tag,
        strain_seed=a.seed,
        strain_levels=list(LEVELS),
        strain_level_meaning="tensor RMS sqrt(trace(eps^2)/3); eps = [e11 e22 e33 e23 e13 e12], "
        "crystal frame; M = (I + eps)^-1, g = O M g0",
        truth_sha256=sha,
        n_by_strain_level={str(L): int(sum(1 for g in gids if level[g] == L)) for L in LEVELS},
        max_eps_roundtrip_error=float(np.abs(eps_back - eps_vox).max()),
    )
    info["grains"] = [
        dict(g, eps=eps[int(g["id"])].tolist(), strain_level=level[int(g["id"])])
        for g in base["grains"]
    ]
    (out_dir / "sample.json").write_text(json.dumps(info) + "\n")
    np.save(out_dir / "grainmap.npy", gmap)
    print(
        json.dumps(
            dict(
                tag=a.tag,
                n_voxels=n,
                n_grains=len(gids),
                n_by_strain_level=info["n_by_strain_level"],
                max_eps_roundtrip_error=info["max_eps_roundtrip_error"],
                sha256=sha[:12],
            )
        )
    )


if __name__ == "__main__":
    main()
