"""Task B step 1 (2-D slices): the cost and the number of frame changes on a plane through the
truth.

A plane in orientation space is spanned by two orthonormal rotation axes (random, seed 2000 +
voxel);
the orientation at (x, y) degrees is rotvec(x u + y v) applied to the truth. Grid: 61 x 61 points
over
+-0.05 degree and over +-0.5 degree, for the six interior full-sample voxels, complete cost (pixel
radius 0, all reflections) with the lenient, C++ and ray (2-pixel window) rules. Saves the cost, the
number of peaks evaluated, the number of peaks whose frame differs from the truth's, and the number
of
peaks with a different pixel signature from the truth's.

  uv run python scripts/eval_reduction/cost_edges_slices.py
"""

import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cost_edges as CE  # noqa: E402
import cost_edges_run as RUN  # noqa: E402

HALF = {"fine": 0.05, "wide": 0.5}
N = 61


def main() -> None:
    RUN._init("full")
    out = RUN.OUT / "slices"
    out.mkdir(parents=True, exist_ok=True)
    for v in RUN.cases("full"):
        base, Rt, verts, phase = RUN.case_data(v)
        rng = np.random.default_rng(2000 + v)
        u = rng.normal(size=3)
        u /= np.linalg.norm(u)
        w = rng.normal(size=3)
        w -= (w @ u) * u
        w /= np.linalg.norm(w)
        res = {"u": u, "v": w}
        for scale, half in HALF.items():
            xs = np.linspace(-half, half, N)
            X, Y = np.meshgrid(xs, xs)
            rv = np.radians(X)[..., None] * u + np.radians(Y)[..., None] * w
            Rs = Rotation.from_rotvec(rv.reshape(-1, 3)).as_matrix() @ Rt
            for rule in ("old", "cpp", "r2"):
                cf = CE.make_cf(base, "full", rule)
                t0 = CE.trace_eval(cf, Rt.astype(np.float32), verts, phase)
                cost = np.zeros(len(Rs))
                nfr = np.zeros(len(Rs), np.int16)
                npx = np.zeros(len(Rs), np.int16)
                npres = np.zeros(len(Rs), np.int16)
                for i, R in enumerate(Rs):
                    tr = CE.trace_eval(cf, R.astype(np.float32), verts, phase)
                    cost[i] = 1.0 - tr["agg"]["quality"]
                    both = (tr["status"] == CE.PRESENT) & (t0["status"] == CE.PRESENT)
                    nfr[i] = int((both & (tr["frame"] != t0["frame"])).sum())
                    npx[i] = int((both & (tr["sig"] != t0["sig"]).any(axis=(1, 2, 3))).sum())
                    npres[i] = int((tr["status"] == CE.PRESENT).sum())
                res[f"{scale}|{rule}|cost"] = cost.reshape(N, N)
                res[f"{scale}|{rule}|frame_changed"] = nfr.reshape(N, N)
                res[f"{scale}|{rule}|pixel_changed"] = npx.reshape(N, N)
                res[f"{scale}|{rule}|present"] = npres.reshape(N, N)
            res[f"{scale}|axis_deg"] = xs
        np.savez_compressed(out / f"v{v}.npz", **res)
        print("slices", v, flush=True)


if __name__ == "__main__":
    main()
