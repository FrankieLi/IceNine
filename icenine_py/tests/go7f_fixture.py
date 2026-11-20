"""Fixed small render and fixed predictions of the GO7f bit-identical-defaults test.

``fingerprint()`` renders 3 voxels of Example2.ManyGrains through the batched path (recording each
triangle drawn) and predicts their spots with ``SpotModel``; ``tests/data/go7f_baseline.npz`` holds
its output from before the PSF and rocking options existed. Regenerate only for a deliberate
change of the default model: ``uv run python tests/go7f_fixture.py tests/data/go7f_baseline.npz``.
"""

import os
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
EX = HERE.parents[1] / "Examples" / "Example2.ManyGrains"


def fingerprint(n_vox: int = 3) -> Dict[str, np.ndarray]:
    from test_strain_render import _render, _setup

    import icenine.grain_fit as GF
    from icenine.reconstructor import _get_voxel_vertices
    from icenine.smooth_cost import SmoothCost

    cwd = os.getcwd()
    os.chdir(EX)
    try:
        sim, sample, dets, rmap, n_om = _setup(n_vox)
        log = _render(sim, sample, dets, rmap, n_om)
        out: Dict[str, np.ndarray] = {
            "tri_peak": np.array([[x[0][0], x[0][1], x[0][2], x[1]] for x in log], np.float64),
            "tri_xy": np.stack([x[2] for x in log]),
            "tri_int": np.array([x[3] for x in log], np.float64),
        }
        sc = SmoothCost(
            sim.simulator,
            dets,
            rmap,
            None,
            sample,
            sample.get_structure_list(),
            eta_limit=sim.exp_setup.get_eta_limit(),
        )
        mic = sample.get_mic()
        model = GF.SpotModel(sc, phase=int(mic.voxels[0].phase))
        Rs = torch.as_tensor(np.stack([v.orientation for v in mic.voxels]), dtype=torch.float64)
        Vs = torch.as_tensor(
            np.stack([_get_voxel_vertices(v).numpy() for v in mic.voxels]), dtype=torch.float64
        )
        rows = model.predict(Rs, Vs)
        for k in ("vox", "refl", "branch", "fpos", "frame", "col", "row", "hit"):
            out["rows_" + k] = getattr(rows, k).numpy()
        bary = torch.as_tensor(GF.subtriangle_bary(3))
        for det in range(model.n_det):
            pc, pr, m, ok = GF._rows_points(rows, bary, det)
            out[f"pts{det}_c"], out[f"pts{det}_r"] = pc.numpy(), pr.numpy()
            out[f"pts{det}_m"], out[f"pts{det}_ok"] = m.numpy(), ok.numpy()
        return out
    finally:
        os.chdir(cwd)


if __name__ == "__main__":
    import icenine

    print("icenine from", icenine.__file__)
    np.savez_compressed(sys.argv[1], **fingerprint())
